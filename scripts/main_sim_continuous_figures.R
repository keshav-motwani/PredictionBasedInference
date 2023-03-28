library(ggplot2)

result_path = "results_fixed_train_100rep/"
n_traintests = n_traintest = 300

if (grepl("fixed", result_path)) {
  n_sim = 100
  var_label = expression(E(Var(hat(beta)[1]~"|"~hat(f))))
  bias_label = expression(E(E(hat(beta)[1] - beta[1] ~"|"~hat(f))))
  bias2_label = expression(E((E(hat(beta)[1] - beta[1] ~"|"~hat(f)))^2))
} else {
  n_sim = 1000
  var_label = expression(Var(hat(beta)[1]))
  bias_label = expression(E(hat(beta)[1] - beta[1]))
  bias2_label = expression((E(hat(beta)[1] - beta[1]))^2)
}

methods = c("naive", "der-postpi", "bs-postpi-par", "bs-postpi-nonpar", "predpowinf", "val*", "observed")
methods_new = c("Naive post-PI", "'Corrected' post-PI, analytical", "'Corrected' post-PI, 'parametric bootstrap'", "'Corrected' post-PI, 'nonparametric bootstrap'", "Prediction-Powered Inference", "Oracle, using validation y", "Use test data")
names(methods_new) = methods
beta1 = 0
pvals = readRDS(paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_ntraintest_", n_traintest, ".rds"))[[6]]
pvals = pvals[pvals$method %in% methods, ]
pvals$method = factor(methods_new[pvals$method], levels = methods_new)
# pvals$rep = rep(1:100, each = nrow(pvals) / 100)

pvals$theoretical = 0
for (n_val in unique(pvals$n_val)) {
  for (method in unique(pvals$method)) {
    pvals[pvals$n_val == n_val & pvals$method == method, "theoretical"] = qunif(ppoints(length(pvals[pvals$n_val == n_val & pvals$method == method, "theoretical"])))[rank(pvals[pvals$n_val == n_val & pvals$method == method, "p"])]
  }
}

# if (grepl(result_path, "fixed")) {
# 
#   indices = pvals$rep == 1
#   pvals_single = pvals[indices, ]
#   pvals_single$theoretical = 0
#   for (n_val in unique(pvals_single$n_val)) {
#     for (method in unique(pvals_single$method)) {
#       pvals_single[pvals_single$n_val == n_val & pvals_single$method == method, "theoretical"] = qunif(ppoints(length(pvals_single[pvals_single$n_val == n_val & pvals_single$method == method, "theoretical"])))[rank(pvals_single[pvals_single$n_val == n_val & pvals_single$method == method, "p"])]
#     }
#   }
#   pvals_single = pvals[sample(indices, length(indices)), ]
# ggplot(pvals_single, aes(x = theoretical, y = p, color = method)) +
#   geom_point(size = 0.5) +
#   xlab("Uniform(0, 1) Quantiles") +
#   ylab("Empirical p-value Quantiles") +
#   facet_wrap(~n_val, labeller = label_bquote(n[val]==.(n_val)), nrow = 1) +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   ggsci::scale_color_npg() +
#   coord_fixed() + xlim(0, 1) + ylim(0, 1) + geom_abline(slope=1, intercept=0, col="black", linetype = "dashed") +
#   guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size=2))) +
#   labs(color = "")
# file = paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_qqplot.pdf")
# ggsave(file, height = 3.3, width = 9)
# 
# }

pvals = pvals[sample(1:nrow(pvals), 10000), ]
ggplot(pvals, aes(x = theoretical, y = p, color = method)) +
  geom_point(size = 0.5) +
  xlab("Uniform(0, 1) Quantiles") +
  ylab("Empirical p-value Quantiles") +
  facet_wrap(~n_val, labeller = label_bquote(n[val]==.(n_val)), nrow = 1) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggsci::scale_color_npg() +
  coord_fixed() + xlim(0, 1) + ylim(0, 1) + geom_abline(slope=1, intercept=0, col="black", linetype = "dashed") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size=2))) +
  labs(color = "")
file = paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_qqplot.pdf")
ggsave(file, height = 3.3, width = 9)


beta1 = 1

extract = function(out, which) {
  
  out = as.data.frame(out)
  out$method = rownames(out)
  
  out = tidyr::pivot_longer(out, -c("method"), names_to = "n_val")
  out$n_val = as.numeric(out$n_val)
  out$type = which
  out$method = factor(methods_new[out$method], levels = methods_new)
  
  as.data.frame(out)
  
}

results = list()
for (n_traintest in n_traintests) {
  reported_var = readRDS(paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_ntraintest_", n_traintest, ".rds"))[[1]]
  reported_var = reported_var[intersect(rownames(reported_var), methods), ]
  true_var = readRDS(paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_ntraintest_", n_traintest, ".rds"))[[2]]
  true_var = true_var[intersect(rownames(true_var), methods), ]
  result = rbind(extract(reported_var, "Reported"), extract(true_var, "True"))
  result$beta_1 = beta1
  result$n_traintest = n_traintest
  results = c(results, list(result))
}
var_data = do.call(rbind, results)

newdata = list()
for (n_traintest in unique(var_data$n_traintest)) {
  for (method in setdiff(levels(var_data$method), c("Use test data", "Prediction-Powered Inference"))) {
    for (type in c("Reported")) {
      subsetted = var_data[var_data$n_traintest == n_traintest & var_data$method == method & var_data$type == type, ]
      subsetted$x = 1/subsetted$n_val
      fit = lm(value ~ 0 + x, data = subsetted)
      newdata = c(newdata, list(data.frame(n_val = seq(150, 2400, length.out = 1000), value = predict(fit, data.frame(x = 1/(seq(150, 2400, length.out = 1000)))), method = method, type = type, n_traintest = n_traintest)))
    }
  }
}
newdata = do.call(rbind, newdata)
newdata$method = factor(newdata$method, levels = methods_new)

library(ggplot2)
ggplot(var_data, aes(x = n_val, y = value, shape = type, color = type)) +
  geom_point(size = 2) +
  facet_wrap(~method, nrow = 2) + # labeller = label_bquote({n[train]==n[test]}==.(n_traintest))) +
  geom_line(data = newdata) +
  geom_vline(aes(xintercept = n_traintest), linetype = "dashed") +
  theme_bw() +
  labs(y = var_label, x = expression(n[val]), shape = "Type", color = "Type") +
  ggsci::scale_color_aaas() +
  ylim(0, NA) +
  theme(legend.position = "bottom")
file = paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_se_plot.pdf")
ggsave(file, height = 4.5, width = 9)

results = list()
for (n_traintest in n_traintests) {
  result = readRDS(paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_ntraintest_", n_traintest, ".rds"))[[3]][methods, ]
  result = extract(result, NA)
  result$beta_1 = beta1
  result$n_traintest = n_traintest
  results = c(results, list(result))
}
bias_data = do.call(rbind, results)
bias_data$n_traintest = factor(bias_data$n_traintest, levels = n_traintests)
# bias_data$se = se_data[se_data$type == "True", "value"] / sqrt(n_sim)
ggplot(bias_data, aes(x = n_val, y = value)) + # , color = n_traintest)) +
  geom_point(size = 1) +
  # geom_errorbar(aes(ymax = value + beta1 + se, ymin = value + beta1 - se)) +
  facet_wrap(~method, nrow = 2) +# , labeller = label_bquote({n[train]==n[test]}==.(n_traintest))) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  theme_bw() +
  labs(y = bias_label, x = expression(n[val]), color = expression(n[train]==n[test])) +
  # ggsci::scale_color_aaas() +
  theme(legend.position = "bottom")
file = paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_bias_plot.pdf")
ggsave(file, height = 4, width = 9)

results = list()
for (n_traintest in n_traintests) {
  result = readRDS(paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_ntraintest_", n_traintest, ".rds"))[[4]][methods, ] - readRDS(paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_ntraintest_", n_traintest, ".rds"))[[2]][methods, ] * (n_sim - 1) / n_sim
  result = extract(result, NA)
  result$beta_1 = beta1
  result$n_traintest = n_traintest
  results = c(results, list(result))
}
bias2_data = do.call(rbind, results)
bias2_data$n_traintest = factor(bias2_data$n_traintest, levels = n_traintests)
ggplot(bias2_data, aes(x = n_val, y = value)) + # , color = n_traintest)) +
  geom_point(size = 1) +
  facet_wrap(~method, nrow = 2) +# , labeller = label_bquote({n[train]==n[test]}==.(n_traintest))) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  theme_bw() +
  labs(y = bias2_label, x = expression(n[val]), color = expression(n[train]==n[test])) +
  # ggsci::scale_color_aaas() +
  theme(legend.position = "bottom")
file = paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_bias2_plot.pdf")
ggsave(file, height = 4, width = 9)


results = list()
for (n_traintest in n_traintests) {
  result = readRDS(paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_ntraintest_", n_traintest, ".rds"))[[5]][methods, ]
  result = extract(result, NA)
  result$beta_1 = beta1
  result$n_traintest = n_traintest
  results = c(results, list(result))
}
coverage_data = do.call(rbind, results)
coverage_data = coverage_data[order(coverage_data$n_traintest, decreasing = T), ]
coverage_data$n_traintest = factor(coverage_data$n_traintest, levels = n_traintests)
ggplot(coverage_data, aes(x = n_val, y = value)) + # , color = n_traintest)) +
  geom_point(size = 1) +
  facet_wrap(~method, nrow = 2) +# , labeller = label_bquote({n[train]==n[test]}==.(n_traintest))) +
  geom_hline(aes(yintercept = .95), linetype = "dashed", color = "black") +
  # geom_vline(aes(xintercept = as.numeric(as.character(n_traintest)), color = n_traintest), linetype = "dashed") +
  theme_bw() +
  labs(y = "95% Confidence Interval Coverage", x = expression(n[val]), color = expression(n[train]==n[test])) +
  # ggsci::scale_color_aaas() +
  theme(legend.position = "bottom") +
  ylim(0, 1)
file = paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_coverage_plot.pdf")
ggsave(file, height = 4, width = 9)

