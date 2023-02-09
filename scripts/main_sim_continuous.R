library(furrr)
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

    print(theta_hat_pp)
    print(se)

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
  bs = 100
  set.seed(2019)

  df = c()

  for (i in 1:n_sim){

    ## relationship model between continuous y and y_p
    ## we use gamma model for continuous outcomes gam(y ~ yp)
    rel_model = gam(y ~ s(pred), data = filter(test_data, sim == i))

    data = filter(val_data, sim == i)


    bs_step = future_map(1:bs, .f = function(i){

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

postpi_classical_bs = function(sim_dat_tv){

  train_data = filter(sim_dat_tv, set == "training")
  val_data = filter(sim_dat_tv, set == "validation")

  ## bootstrap bs times from the validation set
  bs = 100
  set.seed(2019)

  df = c()

  for (i in 1:n_sim){
    print(i)


    tr_data = filter(train_data, sim == i)
    data = filter(val_data, sim == i)


    bs_step = future_map(1:bs, .f = function(b){tryCatch({

      bs_idx = sample(1:nrow(data),nrow(data),replace = TRUE)

      bs_data = data[bs_idx,]

      tr_bs_idx = sample(1:nrow(tr_data),nrow(tr_data),replace = TRUE)

      tr_bs_data = tr_data[tr_bs_idx,]

      trained_model = gam(y ~ s(x1)+s(x2)+s(x3) + s(x4),
                          data = tr_bs_data)

      ## get simulated y using predicted y
      pred_y = predict(trained_model, bs_data)
      if (any(abs(pred_y) > 1000)) stop(paste0("too large prediction - ", paste(tail(sort(abs(pred_y)), 10), collapse = ",")))
      bs_data$pred_y = pred_y
      inf_model = lm(pred_y ~ x1, data=bs_data)

      bs_beta = inf_model %>% tidy() %>% filter(term=="x1") %>% pull(estimate)
      model_se = inf_model %>% tidy() %>% filter(term=="x1") %>% pull(std.error)
      if(is.na(bs_beta) | is.na(model_se)) print(c(bs_beta, model_se, b))
      df = data.frame(bs_beta = bs_beta, model_se = model_se)

    }, error = function(e) {
      print(e)
      print("RETURNING NA")
      data.frame(bs_beta = NA, model_se = NA)
    })}) %>% do.call(rbind,.)

    bs_beta = median(bs_step$bs_beta, na.rm = TRUE)

    ## parametric method
    model_se = median(bs_step$model_se, na.rm = TRUE)
    model_t = bs_beta / model_se
    model_p = 2*pt(-abs(model_t), df = nrow(data) - 1 - 1)


    ## non-parametric method
    sd_se = sd(bs_step$bs_beta, na.rm = TRUE)
    sd_t = bs_beta / sd_se
    sd_p = 2*pt(-abs(sd_t), df = nrow(data) - 1 - 1)

    df = rbind(df, data.frame(sim = i,
                              classical_bs_beta = bs_beta,
                              classical_model_se = model_se,
                              classical_model_t = model_t,
                              classical_model_p = model_p,
                              classical_sd_se = sd_se,
                              classical_sd_t = sd_t,
                              classical_sd_p = sd_p))

  }

  df

}

plan(multicore, workers = 32)

n_sim = 100

beta2 = 0.5
beta3 = 3
beta4 = 4

set.seed(2019)

n_traintests = 300 # c(300, 600, 1200)
n_vals = c(150, 300, 600, 1200, 2400)
beta1s = c(0, 1) # , 3, 5)

methods = c("naive", "der-postpi", "bs-postpi-par", "bs-postpi-nonpar", "val*", "observed", "predpowinf") # , "bs-classical", "observed")

variance_results = as.list(beta1s)
names(variance_results) = paste0("beta1_", beta1s)
bias_results = mse_results = coverage_results = p_value_results = variance_results

coverage = function(true, est, se) {
  mean((true >= est - 1.96 * se) & (true <= est + 1.96 * se))
}

for (k in 1:length(n_traintests)) {

for (j in 1:length(beta1s)){

  var_result = matrix(NA, nrow = length(methods), ncol = length(n_vals))
  colnames(var_result) = as.character(n_vals)
  rownames(var_result) = methods

  bias_result = mse_result = coverage_result = var_result

  p_value_result = list()

  n_traintest = n_traintests[k]
  beta1 = beta1s[j]

  for (i in 1:length(n_vals)) {

    n_val = n_vals[i]

    print(beta1)

    sim_dat_tv = postpi_sim(c(n_traintest, n_traintest, n_val), n_sim, beta1, beta2, beta3,beta4)

    test_data = filter(sim_dat_tv, set == "testing")
    correlation = cor(test_data$pred, test_data$y)
    print(correlation)

    # classical_bs_df = postpi_classical_bs(sim_dat_tv)
    predpowinf_df = predpowinf(sim_dat_tv)
    truth_nc_df = truth_and_nc(sim_dat_tv)
    der_df = postpi_der(sim_dat_tv)
    bs_df = postpi_bs(sim_dat_tv)

    df = cbind(truth_nc_df,der_df,bs_df,predpowinf_df,correlation=correlation, beta1 = beta1) # ,classical_bs_df

var_result["naive", i] = paste0(round(mean(df$nc_sd), 3), "/",  round(sd(df$nc_beta), 3))
var_result["der-postpi", i] = paste0(round(mean(df$der_se), 3), "/",  round(sd(df$der_beta), 3))
var_result["bs-postpi-par", i] = paste0(round(mean(df$model_se), 3), "/",  round(sd(df$bs_beta), 3))
var_result["bs-postpi-nonpar", i] = paste0(round(mean(df$sd_se), 3), "/",  round(sd(df$bs_beta), 3))
# var_result["bs-classical", i] = paste0(round(mean(df$classical_sd_se), 3), "/",  round(sd(df$classical_bs_beta), 3))
var_result["val*", i] = paste0(round(mean(df$truth_sd), 3), "/",  round(sd(df$truth_beta), 3))
var_result["observed", i] = paste0(round(mean(df$observed_sd), 3), "/",  round(sd(df$observed_beta), 3))
var_result["predpowinf", i] = paste0(round(mean(df$predpowinf_se), 3), "/",  round(sd(df$predpowinf_beta), 3))

bias_result["naive", i] = round(mean(df$nc_beta) - beta1, 3)
bias_result["der-postpi", i] = round(mean(df$der_beta) - beta1, 3)
bias_result["bs-postpi-par", i] = round(mean(df$bs_beta) - beta1, 3)
bias_result["bs-postpi-nonpar", i] = round(mean(df$bs_beta) - beta1, 3)
# bias_result["bs-classical", i] = round(mean(df$classical_bs_beta) - beta1, 3)
bias_result["val*", i] = round(mean(df$truth_beta) - beta1, 3)
bias_result["observed", i] = round(mean(df$observed_beta) - beta1, 3)
bias_result["predpowinf", i] = round(mean(df$predpowinf_beta) - beta1, 3)

mse_result["naive", i] = round(mean((df$nc_beta - beta1)^2), 3)
mse_result["der-postpi", i] = round(mean((df$der_beta - beta1)^2), 3)
mse_result["bs-postpi-par", i] = round(mean((df$bs_beta - beta1)^2), 3)
mse_result["bs-postpi-nonpar", i] = round(mean((df$bs_beta - beta1)^2), 3)
# mse_result["bs-classical", i] = round(mean((df$classical_bs_beta - beta1)^2), 3)
mse_result["val*", i] = round(mean((df$truth_beta - beta1)^2), 3)
mse_result["observed", i] = round(mean((df$observed_beta - beta1)^2), 3)
mse_result["predpowinf", i] = round(mean((df$predpowinf_beta - beta1)^2), 3)

coverage_result["naive", i] = round(100 * coverage(beta1, df$nc_beta, df$nc_sd), 1)
coverage_result["der-postpi", i] = round(100 * coverage(beta1, df$der_beta, df$der_se), 1)
coverage_result["bs-postpi-par", i] = round(100 * coverage(beta1, df$bs_beta, df$model_se), 1)
coverage_result["bs-postpi-nonpar", i] = round(100 * coverage(beta1, df$bs_beta, df$sd_se), 1)
# coverage_result["bs-classical", i] = round(100 * coverage(beta1, df$classical_bs_beta, df$classical_sd_se), 1)
coverage_result["val*", i] = round(100 * coverage(beta1, df$truth_beta, df$truth_sd), 1)
coverage_result["observed", i] = round(100 * coverage(beta1, df$observed_beta, df$observed_sd), 1)
coverage_result["predpowinf", i] = round(100 * coverage(beta1, df$predpowinf_beta, df$predpowinf_se), 1)

z_stat = list()
z_stat[["naive"]] = df$nc_beta / df$nc_sd
z_stat[["der-postpi"]] = df$der_beta / df$der_se
z_stat[["bs-postpi-par"]] = df$bs_beta / df$model_se
z_stat[["bs-postpi-nonpar"]] = df$bs_beta / df$sd_se
# z_stat[["bs-classical"]] = df$classical_bs_beta / df$classical_sd_se
z_stat[["val*"]] = df$truth_beta / df$truth_sd
z_stat[["observed"]] = df$observed_beta / df$observed_sd
z_stat[["predpowinf"]] = df$predpowinf_beta / df$predpowinf_se

    p_value_result = c(p_value_result, list(do.call(rbind, lapply(1:length(z_stat), function(i) {
      data.frame(z = z_stat[[i]], method = names(z_stat)[i], p = 2 * (1 - pnorm(abs(z_stat[[i]]))), beta1 = beta1, n_val = n_val)
    }))))

    print(var_result)
    print(bias_result)
    print(mse_result)
    print(coverage_result)

  }

  variance_results[[j]] = var_result
  bias_results[[j]] = bias_result
  mse_results[[j]] = mse_result
  coverage_results[[j]] = coverage_result
  p_value_results[[j]] = do.call(rbind, p_value_result)

  file = paste0("results_observed_test_only/main_postpi_sim_results_beta1_", beta1s[j], "_ntraintest_", n_traintests[k], ".rds")
  saveRDS(list(var_result, bias_result, mse_result, coverage_result, do.call(rbind, p_value_result)), file)

}

}
