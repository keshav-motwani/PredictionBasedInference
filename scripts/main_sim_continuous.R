# library(furrr)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(reshape2)
library(parallel)
library(glue)
library(gam)
library(kableExtra)

postpi_sim = function(ss, n_sim, beta1, beta2, beta3, beta4){


  sim_dat = c()

  for(i in 1:n_sim){

    sim_dat =  rbind(sim_dat,data.frame(x1 = rnorm(sum(ss), mean = 1, sd = 1),

                                        x2 = rnorm(sum(ss), mean = 1, sd = 1),

                                        x3 = rnorm(sum(ss), mean = 1, sd = 1),

                                        x4 = rnorm(sum(ss), mean = 2, sd = 1),

                                        e_g = rnorm(sum(ss),mean = 0, sd=1),

                                        set = rep(c("training","testing","validation"),ss),

                                        sim = i))

  }## end for



  ## Set the ground truth model using x1, x2, x3

  g = function(beta1,beta2,beta3,beta4,x1,x2,x3,x4){

    return(beta1 * x1 + beta2 * x2 + beta3 * smooth(x3) + beta4 * smooth(x4))

  }


  ## simulate y values
  sim_dat = sim_dat %>% mutate(y = g(beta1,beta2,beta3,beta4,x1,x2,x3,x4) + e_g)


  ## use training set to train the model
  ## Get predicted values in the testing and validation set

  train_data = filter(sim_dat, set=="training")
  test_val_data = filter(sim_dat, set == "testing" | set == "validation")

  test_val_data$pred = rep(NA, nrow(test_val_data))

  for( i in 1:n_sim){

    trained_model = gam(y ~ s(x1)+s(x2)+s(x3) + s(x4),
                        data = filter(train_data, sim == i))

    pred_vals = predict(trained_model, newdata = filter(test_val_data, sim == i))

    test_val_data$pred[which(test_val_data$sim == i)] = pred_vals

  }

  train_data$pred = NA

  rbind(train_data, test_val_data)

}

truth_and_nc = function(sim_dat_tv){

  val_data = filter(sim_dat_tv, set == "validation")
  observed_data = filter(sim_dat_tv, set == "testing")

  df = c()

  for (i in 1:n_sim){

    truth_model = lm(y ~ x1, data = filter(val_data, sim == i))
    nc_model = lm(pred ~ x1, data = filter(val_data, sim == i))
    observed_model = lm(y ~ x1, data = filter(observed_data, sim == i))

    truth_beta = truth_model %>% tidy() %>% filter(term=="x1") %>% pull(estimate)
    truth_sd = truth_model %>% tidy() %>% filter(term=="x1") %>% pull(std.error)
    truth_t = truth_model %>% tidy() %>% filter(term=="x1") %>% pull(statistic)
    truth_p = truth_model %>% tidy() %>% filter(term=="x1") %>% pull(p.value)


    nc_beta = nc_model %>% tidy() %>% filter(term=="x1") %>% pull(estimate)
    nc_sd = nc_model %>% tidy() %>% filter(term=="x1") %>% pull(std.error)
    nc_t = nc_model %>% tidy() %>% filter(term=="x1") %>% pull(statistic)
    nc_p = nc_model %>% tidy() %>% filter(term=="x1") %>% pull(p.value)


    observed_beta = observed_model %>% tidy() %>% filter(term=="x1") %>% pull(estimate)
    observed_sd = observed_model %>% tidy() %>% filter(term=="x1") %>% pull(std.error)
    observed_t = observed_model %>% tidy() %>% filter(term=="x1") %>% pull(statistic)
    observed_p = observed_model %>% tidy() %>% filter(term=="x1") %>% pull(p.value)


    df = rbind(df,data.frame(sim = i,
                             truth_beta = truth_beta,
                             truth_sd = truth_sd,
                             truth_t = truth_t,
                             truth_p = truth_p,
                             nc_beta = nc_beta,
                             nc_sd = nc_sd,
                             nc_t = nc_t,
                             nc_p = nc_p,
                             observed_beta = observed_beta,
                             observed_sd = observed_sd,
                             observed_t = observed_t,
                             observed_p = observed_p))

  }

  df

}

predpowinf = function(sim_dat_tv){
  
  test_data = filter(sim_dat_tv, set == "testing")
  test_data$diff = test_data$pred - test_data$y
  val_data = filter(sim_dat_tv, set == "validation")
  
  df = c()
  
  for (i in 1:n_sim){
    
    val_data_i = filter(val_data, sim ==i)
    test_data_i = filter(test_data, sim == i)
    
    val_pred_x = lm(pred ~ x1, data = val_data_i)
    
    rectifier = lm(diff ~ x1, data = test_data_i)
    
    theta_hat_pp = coef(val_pred_x) - coef(rectifier)
    
    X_tilde = model.matrix(~val_data_i$x1)
    X = model.matrix(~test_data_i$x1)
    
    N = nrow(X_tilde)
    
    Sigma_tilde = 1 / N * crossprod(X_tilde)
    
    M_tilde = matrix(0, nrow = ncol(X_tilde), ncol = ncol(X_tilde))
    for (j in 1:N) {
      M_tilde = M_tilde + (val_data_i$pred[j] - crossprod(X_tilde[j, ], coef(val_pred_x))[1, 1])^2 * tcrossprod(X_tilde[j, ])
    }
    M_tilde = M_tilde / N
    
    Sigma_tilde_inv = solve(Sigma_tilde)
    V_tilde = Sigma_tilde_inv %*% M_tilde %*% Sigma_tilde_inv
    
    n = nrow(X)
    
    Sigma = 1 / n * crossprod(X)
    
    M = matrix(0, nrow = ncol(X), ncol = ncol(X))
    for (j in 1:n) {
      M = M + (test_data_i$pred[j] - test_data_i$y[j] - crossprod(X[j, ], coef(rectifier))[1, 1])^2 * tcrossprod(X[j, ])
    }
    M = M / n
    
    Sigma_inv = solve(Sigma)
    V = Sigma_inv %*% M %*% Sigma_inv
    
    se = sqrt((V[2, 2] / n) + (V_tilde[2, 2] / N))

    ## results
    df = rbind(df, data.frame(sim = i,
                              predpowinf_beta = theta_hat_pp[2],
                              predpowinf_se = se,
                              predpowinf_t = NA,
                              predpowinf_p = NA))
    
  }
  
  df
}

postpi_der = function(sim_dat_tv){



  test_data = filter(sim_dat_tv, set == "testing")
  val_data = filter(sim_dat_tv, set == "validation")

  df = c()

  for (i in 1:n_sim){


    val_pred_x = lm(pred ~ x1, data = filter(val_data, sim ==i))
    beta_val_pred = val_pred_x %>% tidy() %>% filter(term=="x1") %>% pull(estimate)


    ## correct variance using relationship moddel estimated in the testing set
    rel_model = lm(y ~ pred, data = filter(test_data, sim == i))

    gamma1 = rel_model %>% tidy() %>% filter(term=="pred") %>% pull(estimate)

    inf_factor = sigma(rel_model)^2 + gamma1^2 * (sigma(val_pred_x)^2)

    mod_matrix = model.matrix(~ x1,data=filter(val_data, sim ==i))

    der_se = sqrt(solve(t(mod_matrix) %*% mod_matrix)[2,2]*inf_factor)


    ## derived adjusted beta
    der_beta = gamma1 * beta_val_pred
    der_t = der_beta / der_se
    der_p = 2*pt(-abs(der_t), df = nrow(filter(val_data, sim ==i)) - 1 - 1)

    ## results
    df = rbind(df, data.frame(sim = i,
                              der_beta = der_beta,
                              der_se = der_se,
                              der_t = der_t,
                              der_p = der_p))

  }




  df
}

postpi_bs = function(sim_dat_tv){


  test_data = filter(sim_dat_tv, set == "testing")
  val_data = filter(sim_dat_tv, set == "validation")

  ## bootstrap bs times from the validation set
  bs = 10
  set.seed(2019)

  df = c()

  for (i in 1:n_sim){

    ## relationship model between continuous y and y_p
    ## we use gamma model for continuous outcomes gam(y ~ yp)
    rel_model = gam(y ~ s(pred), data = filter(test_data, sim == i))

    data = filter(val_data, sim == i)


    bs_step = map(1:bs, .f = function(i){

      bs_idx = sample(1:nrow(data),nrow(data),replace = TRUE)

      bs_data = data[bs_idx,]


      ## get simulated y using predicted y
      sim_y = rnorm(nrow(bs_data), mean=predict(rel_model, bs_data), sd=sigma(rel_model))
      bs_data$sim_y = sim_y

      inf_model = lm(sim_y ~ x1, data=bs_data)

      bs_beta = inf_model %>% tidy() %>% filter(term=="x1") %>% pull(estimate)
      model_se = inf_model %>% tidy() %>% filter(term=="x1") %>% pull(std.error)

      df = data.frame(bs_beta = bs_beta, model_se = model_se)

    }) %>% do.call(rbind,.)



    bs_beta = median(bs_step$bs_beta)

    ## parametric method
    model_se = median(bs_step$model_se)
    model_t = bs_beta / model_se
    model_p = 2*pt(-abs(model_t), df = nrow(data) - 1 - 1)


    ## non-parametric method
    sd_se = sd(bs_step$bs_beta)
    sd_t = bs_beta / sd_se
    sd_p = 2*pt(-abs(sd_t), df = nrow(data) - 1 - 1)

    df = rbind(df, data.frame(sim = i,
                              bs_beta = bs_beta,
                              model_se = model_se,
                              model_t = model_t,
                              model_p = model_p,
                              sd_se = sd_se,
                              sd_t = sd_t,
                              sd_p = sd_p))


  }


  df

}

plan(multicore, workers = 32)

n_sim = 1000

beta2 = 0.5
beta3 = 3
beta4 = 4

set.seed(2019)

n_trains = 300 # c(300, 3000, 30000)
n_vals = c(1000, 2000, 4000, 8000, 16000)
beta1s = c(0, 1) # c(0, 1) # , 3, 5)

methods = c("naive", "der-postpi", "bs-postpi-par", "bs-postpi-nonpar", "val*", "observed", "predpowinf")

coverage = function(true, est, se) {
  ((true >= est - 1.96 * se) & (true <= est + 1.96 * se))
}

for (k in 1:length(n_trains)) {
  
  n_train = n_trains[k]
  
  for (j in 1:length(beta1s)) {
    
    beta1 = beta1s[j]
    
    result = list()
      
      for (i in 1:length(n_vals)) {
        
        n_val = n_vals[i]
        n_test = n_val * 0.1
        
        print(beta1)
        
        sim_dat_tv = postpi_sim(c(n_train, n_test, n_val), n_sim, beta1, beta2, beta3, beta4)

        predpowinf_df = predpowinf(sim_dat_tv)
        truth_nc_df = truth_and_nc(sim_dat_tv)
        der_df = postpi_der(sim_dat_tv)
        bs_df = postpi_bs(sim_dat_tv)
        
        df = cbind(truth_nc_df,der_df,bs_df,predpowinf_df, beta1 = beta1)
        
        estimates = matrix(NA, nrow = n_sim, ncol = length(methods))
        colnames(estimates) = methods
        reported_ses = estimates
        
        estimates[, "naive"] = df$nc_beta
        estimates[, "der-postpi"] = df$der_beta
        estimates[, "bs-postpi-par"] = df$bs_beta
        estimates[, "bs-postpi-nonpar"] = df$bs_beta
        estimates[, "val*"] = df$truth_beta
        estimates[, "observed"] = df$observed_beta
        estimates[, "predpowinf"] = df$predpowinf_beta
        
        reported_ses[, "naive"] = df$nc_sd
        reported_ses[, "der-postpi"] = df$der_se
        reported_ses[, "bs-postpi-par"] = df$model_se
        reported_ses[, "bs-postpi-nonpar"] = df$sd_se
        reported_ses[, "val*"] = df$truth_sd
        reported_ses[, "observed"] = df$observed_sd
        reported_ses[, "predpowinf"] = df$predpowinf_se
        
        for (method in methods) {
          
          reported_var = reported_ses[, method]^2
          bias = estimates[, method] - beta1
          cov = coverage(beta1, estimates[, method], reported_ses[, method])
          z_stat = estimates[, method] / reported_ses[, method]
          p_value = 2 * (1 - pnorm(abs(z_stat)))
          
          result = c(result, list(data.frame(n_train = n_train, n_test = n_test, n_val = n_val, beta1 = beta1, method = method, 
                                             reported_var = reported_var, bias = bias, 
                                             coverage = cov, z_stat = z_stat, p_value = p_value)))
          
        }
        
      }
      
    dir.create("results")
    file = paste0("results/main_postpi_sim_results_beta1_", beta1s[j], "_ntrain_", n_trains[k], ".rds")
    saveRDS(do.call(rbind, result), file)
    
  }
  
}
