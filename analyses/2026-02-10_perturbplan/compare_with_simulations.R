#sce = readRDS("~/GitHub/broad-nnfc-hiPSC-to-EC/results/2026-02-10_simulation_power_analysis/perturb_sce.rds")

library(here)
sobj = readRDS(here("results/2026-04-22_sceptre_results/day0_grna20/sceptre_object.rds"))
sobj
sobj@response_precomputations
theta_vec <- vapply(
  sobj@response_precomputations,
  function(x) x$theta,
  numeric(1)
)
theta_vec

sceptre_object =sobj

#This reproduces the dispersion
gene_id <- "ZWILCH"
y <- sceptre::get_response_matrix(sceptre_object)[gene_id, ][sceptre_object@cells_in_use]
X <- sceptre_object@covariate_matrix[sceptre_object@cells_in_use, , drop = FALSE]

precomp_fun <- getFromNamespace("perform_response_precomputation", "sceptre")
out <- precomp_fun(expressions = as.numeric(y), covariate_matrix = X)

c(
  manual = out$theta,
  cached = sceptre_object@response_precomputations[[gene_id]]$theta
)

precomp <- sceptre_object@response_precomputations[[gene_id]]

X <- sceptre_object@covariate_matrix[
  sceptre_object@cells_in_use,
  ,
  drop = FALSE
]

mu <- exp(as.numeric(X %*% precomp$fitted_coefs))
names(mu) <- colnames(sceptre_object@response_matrix[[1]])[
  sceptre_object@cells_in_use
]

mu
y<- sceptre_object@response_matrix[[1]][gene_id, sceptre_object@cells_in_use]
plot(cbind(observed = as.numeric(y), fitted_mean = mu))

#

perturbplan = read.table(here("results/2026-04-24_perturbplan/day0/expression_stats_mean_theta.tsv"), header=TRUE, sep="\t")
plot(perturbplan$expression_theta, theta_vec[perturbplan$response_id])
# remove outlier from perturbplan$expression_theta
perturbplan_no_outlier = perturbplan[perturbplan$expression_theta < 10, ]
theta_vec_no_outlier = theta_vec[theta_vec < 50]
# remove na
theta_vec_no_outlier = theta_vec_no_outlier[!is.na(theta_vec_no_outlier)]
# keep common genes
common_genes = intersect(perturbplan_no_outlier$response_id, names(theta_vec_no_outlier))
# Plot
plot(perturbplan_no_outlier$expression_theta[perturbplan_no_outlier$response_id %in% common_genes], 
     theta_vec_no_outlier[common_genes], 
     xlab="perturbplan", 
     ylab="sceptre", 
     main="Per-gene computed theta")


gene_mean <- rowData(sce)$mean

length(gene_mean)

rowData(sce)

sce_df = rowData(sce)
sce_df$disp = unlist(sce_df$dispersion)


perturb_sim = readH5MU("results/2026-02-10_perturbplan/day2/perturbplan_es_15_pval_0.001_both_mean_dispersion_from_sim.h5mu")
perturbplan_sim_power = perturb_sim@metadata$power_results$individual_power
names(perturbplan_sim_power)
plot_sim_df = as.data.frame(perturbplan_sim_power)

mask_unique = !duplicated(plot_sim_df[, c("expression_mean","expression_dispersion", "response_id")])
plot_sim_df = plot_sim_df[mask_unique,c("response_id", "expression_mean","expression_dispersion")]
names(plot_sim_df)
head(plot_sim_df)

head(perturbplan_sim_power)
head(perturbplan_power)


plot(perturbplan_sim_power$power, perturbplan_power$power, xlab="perturbplan with mean/dispersion from simulation", ylab="perturbplan with mean/dispersion from SCEPTRE", main="Power scores for different perturbplan runs")

perturb_results = readH5MU("~/GitHub/broad-nnfc-hiPSC-to-EC/results/2026-02-10_perturbplan/day0/perturbplan_es_15_pval_0.001_both.h5mu")
perturbplan_power = perturb_results@metadata$power_results$individual_power
names(perturbplan_power)
plot_df = as.data.frame(perturbplan_power)
mask_unique = !duplicated(plot_df[, c("expression_mean","expression_dispersion", "response_id")])
plot_df = plot_df[mask_unique,c("response_id", "expression_mean","expression_dispersion")]
names(plot_df)
head(plot_df)

mask = intersect(plot_df$response_id, names(sce_df$dispersion))

plot_df = plot_df[plot_df$response_id %in% mask, ]
dim(plot_df)
plot_df$sim_mean = sce_df[plot_df$response_id,"mean"]
plot_df$sim_dispersion = unlist(sce_df$dispersion[plot_df$response_id])

plot(plot_df$sim_dispersion, 
     plot_df$expression_dispersion, 
     xlab="simulation", 
     ylab="perturbplan", 
     main="Per-gene computed dispersion(log scale)", log="xy"
     )
abline(0,1)

plot(plot_df$sim_mean,
     plot_df$expression_mean,
     xlab="simulation",
     ylab="perturbplan",
     main="Per-gene computed mean"
     )
abline(0,1)
