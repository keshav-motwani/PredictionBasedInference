library(ggplot2)

result_path = "results_fixed_train/"
n_train = 300

methods = c("naive", "der-postpi", "bs-postpi-par", "bs-postpi-nonpar", "predpowinf", "observed", "val*")
methods_new = c("Naive", "postPI, analytical", "postPI, 'parametric bootstrap'", "postPI, 'nonparametric bootstrap'", "Prediction-Powered Inference", "Classical, using labeled data", "Oracle, using unlabeled y")
names(methods_new) = methods

beta1 = 0
result = readRDS(paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_ntrain_", n_train, ".rds"))
result = result[result$method %in% methods, ]
result$method = factor(methods_new[result$method], levels = methods_new)
result$n_val = factor(result$n_val, levels = sort(unique(result$n_val)))

result$theoretical = 0

for (n_val in unique(result$n_val)) {
  for (method in unique(result$method)) {
    for (r in 1:max(result$rep)) {
      result[result$n_val == n_val & result$method == method & result$rep == r, "theoretical"] = qunif(ppoints(length(result[result$n_val == n_val & result$method == method & result$rep == r, "theoretical"])))[rank(result[result$n_val == n_val & result$method == method & result$rep == r, "p_value"])]
    }
  }
}

result$expr_n_val = factor(paste0("{n[lab]==", 0.1 * as.numeric(as.character(result$n_val)), "} / n[unlab]==", result$n_val), levels = paste0("{n[lab]==", 0.1 * as.numeric(levels(result$n_val)), "} / n[unlab]==", levels(result$n_val)))
ggplot(result, aes(x = theoretical, y = p_value, color = method)) +
  geom_point(size = 0.5) +
  xlab("Uniform(0, 1) Quantiles") +
  ylab("Empirical p-value Quantiles") +
  facet_grid(paste0("GAM ", rep)~expr_n_val, labeller = labeller(expr_n_val = label_parsed)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggsci::scale_color_npg() +
  coord_fixed() + xlim(0, 1) + ylim(0, 1) + geom_abline(slope=1, intercept=0, col="black", linetype = "dashed") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size=2))) +
  labs(color = "")
file = paste0(result_path, "/main_fixed_train_postpi_sim_results_beta1_", beta1, "_qqplot.pdf")
ggsave(file, height = 6.5, width = 9)

for (beta1 in c(0, 1)) {
  
  result = readRDS(paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_ntrain_", n_train, ".rds"))
  result = result[result$method %in% methods, ]
  result$method = factor(methods_new[result$method], levels = methods_new)
  result$n_val = factor(result$n_val, levels = sort(unique(result$n_val)))
  result$rep = paste0("GAM ", result$rep)
  
  if (beta1 == 1) {
  
    coverage = result %>%
      group_by(n_val, method, rep) %>%
      summarize(coverage = mean(coverage))
    
    ggplot(coverage, aes(x = as.numeric(as.character(n_val)), y = coverage, color = rep)) +
      geom_point(size = 1) +
      geom_line() +
      facet_wrap(~method, nrow = 2) +
      geom_hline(aes(yintercept = .95), linetype = "dashed", color = "black") +
      theme_bw() +
      labs(y = "Empirical Coverage", x = expression(n[unlab]), color = "Prediction Model") +
      ggsci::scale_color_aaas() +
      theme(legend.position = "bottom") +
      ylim(0, 1)
    file = paste0(result_path, "/main_fixed_train_postpi_sim_results_beta1_", beta1, "_coverage_plot.pdf")
    ggsave(file, height = 4.5, width = 9)
  
  }
  
  assumed = result %>% select(n_val, rep, method) %>% unique()
  assumed = do.call(rbind, replicate(100, assumed, simplify = FALSE))
  z = seq(-3, 3, length.out = 100)
  assumed$z = rep(z, each = nrow(assumed) / 100)
  
  plot = ggplot(result, aes(x = n_val, y = bias / sqrt(reported_var), color = rep)) +
    facet_wrap(~method, nrow = 2) +
    theme_bw() +
    ggsci::scale_color_aaas() +
    ggforce::geom_sina(scale = "area", size = 0.125) +
    geom_violin(data = assumed, aes(y = z), alpha = 1, fill = NA, show.legend = F) +
    geom_hline(yintercept = -1.96, linetype = "dashed") + 
    geom_hline(yintercept = 1.96, linetype = "dashed") +
    labs(color = "Prediction Model", y = expression(frac(hat(beta)[1]*' - '*beta[1]^'*', hat('SE')*'('*hat(beta[1])*')')), x = expression(n[unlab])) +
    theme(legend.position = "bottom") + 
    guides(colour = guide_legend(override.aes = list(size=1)))
  test = ggplot_build(plot)
  test$data[[2]] = left_join(test$data[[2]], test$data[[1]] %>% mutate(scale = round(sinawidth / density, 3)) %>% select(colour, PANEL, group, scale) %>% unique())
  test$data[[2]]$violinwidth = dnorm(test$data[[2]]$y) * tail(test$data[[2]] %>% filter(PANEL == 7) %>% pull(scale), 1) * 0.9
  test$data[[2]]$colour = "black"
  file = paste0(result_path, "/main_fixed_train_postpi_sim_results_beta1_", beta1, "_distribution_plot.pdf")
  pdf(file, height = 4.5, width = 9)
  plot(ggplot_gtable(test))
  dev.off()

}
