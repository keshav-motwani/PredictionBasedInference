sim_func = function(n_rep=500,n_train=100,n_test=100,n_valid=100){
  
  n_tot = n_train + n_test + n_valid
  
  
  beta_wang = beta_naive = beta_oracle = se_wang = se_naive = se_oracle = rep(NA,n_rep)
  
  for(k in 1:n_rep){
    ## Make some covariates
    x1 = rnorm(n_tot)
    x2 = rnorm(n_tot)
    x3 = rnorm(n_tot)
    err = rnorm(n_tot,sd=4)
    
    ## Make an outcome
    y = 2*x1 + 3*x2^2 + 0.5*x3^3 + err 
    
    ## Create an indicator of training,testing,or validation
    set_ind = c(rep("train",n_train),
                rep("test",n_test),
                rep("valid",n_valid))
    
    ## Build the data set
    dat = data.frame(y,x1,x2,x3,set_ind)
    
    
    ## Fit the model
    train_dat  = dat %>%
      filter(set_ind == "train") %>%
      select(y,x1,x2,x3)
    
    pred_mod = gam(y ~ s(x1) + s(x2) + s(x3))
    
    ## Add the predictions
    dat$pred = predict(pred_mod,dat)
    
    ## Filter to the test set
    test_dat = dat %>%
      filter(set_ind == "test")
    
    ## Fit the relationship model
    rel_mod = lm(y ~ ns(pred,df=3),data=test_dat)
    
    ## Get the residual variance for use in simulating
    sigma2 = sigma(rel_mod)^2
    
    
    ## Filter to the validation set
    valid_dat = dat %>%
      filter(set_ind == "valid")
    
    ## Get the simulated data
    sim_dat = valid_dat %>% sample_n(size=n_valid,replace=TRUE)
    
    ## Sample simulated y's
    y_sim = rnorm(n_valid,
                  mean=sim_dat$pred,
                  sd = sqrt(sigma2))
    sim_dat = cbind(sim_dat,y_sim)
    
    ## Fit the models
    
    lm_oracle = lm(y ~ x1,data=valid_dat)
    lm_naive = lm(pred ~ x1,data = valid_dat)
    lm_wang = lm(y_sim ~ x1,data=sim_dat)
    
    ## Get the coefficients
    beta_oracle[k] = lm_oracle$coefficients[2]
    beta_naive[k] = lm_naive$coefficients[2]
    beta_wang[k] = lm_wang$coefficients[2]
    
    ## Get the parametric ses
    se_oracle[k] = summary(lm_oracle)$coefficients[2,2]
    se_naive[k] = summary(lm_naive)$coefficients[2,2]
    se_wang[k] = summary(lm_wang)$coefficients[2,2]
    
  }
  
  results = data.frame(beta_oracle,se_oracle,beta_wang,se_wang,beta_naive,se_naive)
  
  results = results %>%
    mutate(up_oracle = beta_oracle + 1.96*se_oracle,
           low_oracle = beta_oracle - 1.96*se_oracle,
           cov_oracle = ((up_oracle > 2) & (low_oracle < 2)),
           ## Wang
           up_wang = beta_wang + 1.96*se_wang,
           low_wang = beta_wang - 1.96*se_wang,
           cov_wang = ((up_wang > 2) & (low_wang < 2)),
           
           ## Naive
           up_naive = beta_naive + 1.96*se_naive,
           low_naive = beta_naive - 1.96*se_naive,
           cov_naive = ((up_naive > 2) & (low_naive < 2)),
    ) 
  
  return(results)
  
  
}