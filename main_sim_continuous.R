library(dplyr)
library(tidyr)
library(caret)
library(ggplot2)
library(broom)
library(reshape2)
library(parallel)
library(glue)
library(gam)

postpi_sim = function(ss, n_sim, beta1, beta2, beta3, beta4){
  
  
  sim_dat = c()
  
  for(i in 1:n_sim){
    
    sim_dat =  rbind(sim_dat,data.frame(x1 = rnorm(ss*3, mean = 1, sd = 1),
                                        
                                        x2 = rnorm(ss*3, mean = 1, sd = 1),
                                        
                                        x3 = rnorm(ss*3, mean = 1, sd = 1),
                                        
                                        x4 = rnorm(ss*3, mean = 2, sd = 1),
                                        
                                        e_g = rnorm(ss*3,mean = 0, sd=1),
                                        
                                        set = rep(c("training","testing","validation"),each=ss),
                                        
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
  
  
  test_val_data
  
}

truth_and_nc = function(sim_dat_tv){
  
  val_data = filter(sim_dat_tv, set == "validation")
  
  df = c()
  
  for (i in 1:n_sim){
    
    truth_model = lm(y ~ x1, data = filter(val_data, sim == i))
    nc_model = lm(pred ~ x1, data = filter(val_data, sim == i))
    
    truth_beta = truth_model %>% tidy() %>% filter(term=="x1") %>% pull(estimate)
    truth_sd = truth_model %>% tidy() %>% filter(term=="x1") %>% pull(std.error)
    truth_t = truth_model %>% tidy() %>% filter(term=="x1") %>% pull(statistic)
    truth_p = truth_model %>% tidy() %>% filter(term=="x1") %>% pull(p.value)
    
    
    nc_beta = nc_model %>% tidy() %>% filter(term=="x1") %>% pull(estimate)
    nc_sd = nc_model %>% tidy() %>% filter(term=="x1") %>% pull(std.error)
    nc_t = nc_model %>% tidy() %>% filter(term=="x1") %>% pull(statistic)
    nc_p = nc_model %>% tidy() %>% filter(term=="x1") %>% pull(p.value)
    
    
    df = rbind(df,data.frame(sim = i, 
                             truth_beta = truth_beta, 
                             truth_sd = truth_sd,
                             truth_t = truth_t,
                             truth_p = truth_p,
                             nc_beta = nc_beta, 
                             nc_sd = nc_sd, 
                             nc_t = nc_t,
                             nc_p = nc_p))
    
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
    
    
    bs_step = lapply(1:bs, function(i){
      
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


ss = 300
n_sim = 300

beta2 = 0.5
beta3 = 3
beta4 = 4


set.seed(2019)

for (beta1 in c(seq(-6,6,1))){
  
  print(beta1)
  
  sim_dat_tv = postpi_sim(ss, n_sim, beta1, beta2, beta3,beta4)
  
  test_data = filter(sim_dat_tv, set == "testing")
  correlation = cor(test_data$pred, test_data$y)
  print(correlation)
  
  truth_nc_df = truth_and_nc(sim_dat_tv)
  der_df = postpi_der(sim_dat_tv)
  bs_df = postpi_bs(sim_dat_tv)
  
  df = cbind(truth_nc_df,der_df,bs_df,correlation=correlation, beta1 = beta1)
  
  saveRDS(df, file = glue("./simulated_con_rda/v4_simcon_beta1_{beta1}.rds"))
  
}

