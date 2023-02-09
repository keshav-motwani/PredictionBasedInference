library(ggplot2)

methods = c("naive", "der-postpi", "bs-postpi-par", "bs-postpi-nonpar", "predpowinf", "val*", "observed")
methods_new = c("Naive post-PI", "'Corrected' post-PI, analytical", "'Corrected' post-PI, 'parametric bootstrap'", "'Corrected' post-PI, 'nonparametric bootstrap'", "Prediction-Powered Inference", "Oracle, using validation y", "Use labeled data")
names(methods_new) = methods
beta1 = 0
n_traintest = 300
pvals = readRDS(paste0("results_observed_test_only/main_postpi_sim_results_beta1_", beta1, "_ntraintest_", n_traintest, ".rds"))[[5]]
pvals = pvals[pvals$method %in% methods, ]
pvals$method = factor(methods_new[pvals$method], levels = methods_new)

pvals$theoretical = 0
for (n_val in unique(pvals$n_val)) {
  for (method in unique(pvals$method)) {
    pvals[pvals$n_val == n_val & pvals$method == method, "theoretical"] = qunif(ppoints(length(pvals[pvals$n_val == n_val & pvals$method == method, "theoretical"])))[rank(pvals[pvals$n_val == n_val & pvals$method == method, "p"])]
  }
}

pvals = pvals[sample(1:nrow(pvals), nrow(pvals)), ]
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
file = paste0("results_observed_test_only/main_postpi_sim_results_beta1_", beta1, "_qqplot.pdf")
ggsave(file, height = 3.3, width = 9)


beta1 = 1
n_traintests = 300
n_sim = 1000

extract = function(result, which, se = FALSE) {

  if (se) {

    if (which == "Reported") i = 1
    else i = 2

    out = matrix(as.numeric(sapply(strsplit(result, "/"), `[`, i)), nrow = nrow(result))
    rownames(out) = rownames(result)
    colnames(out) = colnames(result)

  } else {

    out = result

  }

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
  result = readRDS(paste0("results_observed_test_only/main_postpi_sim_results_beta1_", beta1, "_ntraintest_", n_traintest, ".rds"))[[1]]
  result = result[intersect(rownames(result), methods), ]
  result = rbind(extract(result, "Reported", T), extract(result, "True", T))
  result$beta_1 = beta1
  result$n_traintest = n_traintest
  results = c(results, list(result))
}
se_data = do.call(rbind, results)

newdata = list()
for (n_traintest in unique(se_data$n_traintest)) {
  for (method in setdiff(levels(se_data$method), "Use labeled data")) {
    for (type in c("Reported")) {
      subsetted = se_data[se_data$n_traintest == n_traintest & se_data$method == method & se_data$type == type, ]
      subsetted$x = 1/sqrt(subsetted$n_val)
      fit = lm(value ~ 0 + x, data = subsetted)
      newdata = c(newdata, list(data.frame(n_val = seq(150, 2400, length.out = 1000), value = predict(fit, data.frame(x = 1/sqrt(seq(150, 2400, length.out = 1000)))), method = method, type = type, n_traintest = n_traintest)))
    }
  }
}
newdata = do.call(rbind, newdata)
newdata$method = factor(newdata$method, levels = methods_new)

library(ggplot2)
ggplot(se_data, aes(x = n_val, y = value, shape = type, color = type)) +
  geom_point(size = 2) +
  facet_wrap(~method, nrow = 2) + # labeller = label_bquote({n[train]==n[test]}==.(n_traintest))) +
  geom_line(data = newdata) +
  geom_vline(aes(xintercept = n_traintest), linetype = "dashed") +
  theme_bw() +
  labs(y = "Standard Error", x = expression(n[val]), shape = "Type", color = "Type") +
  ggsci::scale_color_aaas() +
  ylim(0, NA) +
  theme(legend.position = "bottom")
file = paste0("results_observed_test_only/main_postpi_sim_results_beta1_", beta1, "_se_plot.pdf")
ggsave(file, height = 4.5, width = 9)

results = list()
for (n_traintest in n_traintests) {
  result = readRDS(paste0("results_observed_test_only/main_postpi_sim_results_beta1_", beta1, "_ntraintest_", n_traintest, ".rds"))[[2]][methods, ]
  result = extract(result, NA, F)
  result$beta_1 = beta1
  result$n_traintest = n_traintest
  results = c(results, list(result))
}
expected_value_data = do.call(rbind, results)
expected_value_data$n_traintest = factor(expected_value_data$n_traintest, levels = n_traintests)
expected_value_data$se = se_data[se_data$type == "True", "value"] / sqrt(n_sim)
ggplot(expected_value_data, aes(x = n_val, y = value + beta1)) + # , color = n_traintest)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymax = value + beta1 + se, ymin = value + beta1 - se)) +
  facet_wrap(~method, nrow = 2) +# , labeller = label_bquote({n[train]==n[test]}==.(n_traintest))) +
  geom_hline(aes(yintercept = beta1), linetype = "dashed") +
  theme_bw() +
  labs(y = expression(E(hat(beta)[1])), x = expression(n[val]), color = expression(n[train]==n[test])) +
  # ggsci::scale_color_aaas() +
  theme(legend.position = "bottom") +
  ylim(beta1 - 0.1, beta1 + 0.1)
file = paste0("results_observed_test_only/main_postpi_sim_results_beta1_", beta1, "_expected_value_plot.pdf")
ggsave(file, height = 4, width = 9)

results = list()
for (n_traintest in n_traintests) {
  result = readRDS(paste0("results_observed_test_only/main_postpi_sim_results_beta1_", beta1, "_ntraintest_", n_traintest, ".rds"))[[4]][methods, ]
  result = extract(result, NA, F)
  result$beta_1 = beta1
  result$n_traintest = n_traintest
  results = c(results, list(result))
}
coverage_data = do.call(rbind, results)
coverage_data = coverage_data[order(coverage_data$n_traintest, decreasing = T), ]
coverage_data$n_traintest = factor(coverage_data$n_traintest, levels = n_traintests)
ggplot(coverage_data, aes(x = n_val, y = value)) + # , color = n_traintest)) +
  geom_point(size = 2) +
  facet_wrap(~method, nrow = 2) +# , labeller = label_bquote({n[train]==n[test]}==.(n_traintest))) +
  geom_hline(aes(yintercept = 95), linetype = "dashed", color = "black") +
  # geom_vline(aes(xintercept = as.numeric(as.character(n_traintest)), color = n_traintest), linetype = "dashed") +
  theme_bw() +
  labs(y = "95% Confidence Interval Coverage", x = expression(n[val]), color = expression(n[train]==n[test])) +
  # ggsci::scale_color_aaas() +
  theme(legend.position = "bottom") +
  ylim(0, 100)
file = paste0("results_observed_test_only/main_postpi_sim_results_beta1_", beta1, "_coverage_plot.pdf")
ggsave(file, height = 4, width = 9)

