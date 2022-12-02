methods = c("naive", "der-postpi", "bs-postpi-par", "bs-postpi-nonpar", "val*")
beta1 = 1
n_traintests = c(300, 600, 1200)
results = list()
for (n_traintest in n_traintests) {
  result = readRDS(paste0("results/main_postpi_sim_results_beta1_", beta1, "_ntraintest_", n_traintest, ".rds"))[[1]][methods, ]
  result = rbind(extract(result, "Reported"), extract(result, "True"))
  result$beta_1 = beta1
  result$n_traintest = n_traintest
  results = c(results, list(result))
}
data = do.call(rbind, results)

newdata = list()
for (n_traintest in unique(data$n_traintest)) {
  for (method in levels(data$method)) {
    for (type in c("Reported")) {
      subsetted = data[data$n_traintest == n_traintest & data$method == method & data$type == type, ]
      fit = nls(value ~ I(K * n_val^-0.5), data = subsetted, start = list(K = 1))
      newdata = c(newdata, list(data.frame(n_val = seq(150, 2400, length.out = 1000), value = predict(fit, data.frame(n_val = seq(150, 2400, length.out = 1000))), method = method, type = type, n_traintest = n_traintest)))
    }
  }
}
newdata = do.call(rbind, newdata)
newdata$method = factor(newdata$method, levels = methods)

library(ggplot2)
ggplot(data, aes(x = n_val, y = value, shape = type, color = type)) +
  geom_point(size = 2) +
  facet_grid(n_traintest~method, labeller = label_bquote({n[train]==n[test]}==.(n_traintest))) +
  geom_line(data = newdata) +
  geom_vline(aes(xintercept = n_traintest), linetype = "dashed") +
  theme_bw() +
  labs(y = "Standard Error", x = expression(n[val]), shape = "Type", color = "Type") +
  ggsci::scale_color_aaas() +
  ylim(0, NA) +
  theme(legend.position = "bottom")
file = paste0("results/main_postpi_sim_results_beta1_", beta1, "_se_plot.pdf")
ggsave(file, height = 6, width = 9.5)
