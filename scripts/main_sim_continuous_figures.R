library(ggplot2)

result_path = "results/"
n_train = 300

methods = c("naive", "der-postpi", "bs-postpi-par", "bs-postpi-nonpar", "predpowinf", "observed") # , "val*")
methods_new = c("Naive", "Wang et al., analytical", "Wang et al., 'parametric bootstrap'", "Wang et al., 'nonparametric bootstrap'", "Angelopoulos et al.", "Classical, using labeled data") # , "Oracle, using unlabeled y")
names(methods_new) = methods

beta1 = 0
result = readRDS(paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_ntrain_", n_train, ".rds"))
result = result[result$method %in% methods, ]
result$method = factor(methods_new[result$method], levels = methods_new)
result$n_val = factor(result$n_val, levels = sort(unique(result$n_val)))

result$theoretical = 0

for (n_val in unique(result$n_val)) {
  for (method in unique(result$method)) {
    result[result$n_val == n_val & result$method == method, "theoretical"] = qunif(ppoints(length(result[result$n_val == n_val & result$method == method, "theoretical"])))[rank(result[result$n_val == n_val & result$method == method, "p_value"])]
  }
}

result$expr_n_val = factor(paste0("{n[lab]==", 0.1 * as.numeric(as.character(result$n_val)), "} / n[unlab]==", result$n_val), levels = paste0("{n[lab]==", 0.1 * as.numeric(levels(result$n_val)), "} / n[unlab]==", levels(result$n_val)))
ggplot(result, aes(x = theoretical, y = p_value, color = method)) +
  geom_point(size = 0.5) +
  xlab("Uniform(0, 1) Quantiles") +
  ylab("Empirical p-value Quantiles") +
  facet_grid(~expr_n_val, labeller = labeller(expr_n_val = label_parsed)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ggsci::scale_color_npg() +
  coord_fixed() + xlim(0, 1) + ylim(0, 1) + geom_abline(slope=1, intercept=0, col="black", linetype = "dashed") +
  guides(color = guide_legend(nrow = 2, byrow = F, override.aes = list(size=2))) +
  labs(color = "")
file = paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_qqplot.pdf")
ggsave(file, height = 3.5, width = 9)


beta1 = 1
result = readRDS(paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_ntrain_", n_train, ".rds"))
result = result[result$method %in% methods, ]
result$method = factor(methods_new[result$method], levels = methods_new)
result$n_val = factor(result$n_val, levels = sort(unique(result$n_val)))

coverage = result %>%
  group_by(n_val, method) %>%
  summarize(coverage = mean(coverage))

ggplot(coverage, aes(x = as.numeric(as.character(n_val)), y = coverage, color = method)) +
  geom_point(size = 1) +
  geom_line(show.legend = FALSE) +
  geom_hline(aes(yintercept = .95), linetype = "dashed", color = "black") +
  theme_bw() +
  labs(y = "Empirical Coverage", x = expression(n[unlab]), color = "") +
  ggsci::scale_color_npg() +
  theme(legend.position = "right") +
  ylim(0, 1) +
  guides(color = guide_legend(ncol = 1, byrow = F, override.aes = list(size=2)))
file = paste0(result_path, "/main_postpi_sim_results_beta1_", beta1, "_coverage_plot.pdf")
ggsave(file, height = 3, width = 6)

