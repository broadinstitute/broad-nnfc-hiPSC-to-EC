function (discovery_pairs, cells_per_grna, baseline_expression_stats, 
    control_group, fold_change_mean, fold_change_sd, num_total_cells = NULL, 
    cutoff = NULL, n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L, 
    side = "both", multiple_testing_method = "BH", multiple_testing_alpha = 0.1) 
{
    input_check_posthoc(discovery_pairs = discovery_pairs, cells_per_grna = cells_per_grna, 
        baseline_expression_stats = baseline_expression_stats, 
        control_group = control_group, fold_change_mean = fold_change_mean, 
        fold_change_sd = fold_change_sd, num_total_cells = num_total_cells, 
        cutoff = cutoff, n_nonzero_trt_thresh = n_nonzero_trt_thresh, 
        n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh, side = side, 
        multiple_testing_method = multiple_testing_method, multiple_testing_alpha = multiple_testing_alpha)
    enhancer_gene <- dplyr::left_join(dplyr::left_join(discovery_pairs, 
        dplyr::ungroup(dplyr::summarize(dplyr::group_by(dplyr::filter(cells_per_grna, 
            grna_target != "non-targeting"), grna_target), num_trt_cells = sum(num_cells), 
            num_trt_cells_sq = sum(num_cells^2))), "grna_target", 
        relationship = "many-to-one"), baseline_expression_stats, 
        "response_id", relationship = "many-to-one")
    if (control_group == "nt_cells") {
        num_cntrl_cells <- dplyr::pull(dplyr::summarize(dplyr::filter(cells_per_grna, 
            grna_target == "non-targeting"), sum(num_cells)))
        enhancer_gene <- dplyr::mutate(enhancer_gene, num_cntrl_cells = num_cntrl_cells)
    }
    else {
        enhancer_gene <- dplyr::mutate(enhancer_gene, num_cntrl_cells = num_total_cells - 
            num_trt_cells)
    }
    if (is.numeric(fold_change_mean)) {
        enhancer_gene <- dplyr::mutate(enhancer_gene, fold_change_mean = fold_change_mean, 
            fold_change_sd = fold_change_sd)
    }
    else {
        enhancer_gene <- dplyr::left_join(dplyr::left_join(enhancer_gene, 
            fold_change_mean, c("grna_target", "response_id"), 
            relationship = "one-to-one"), fold_change_sd, c("grna_target", 
            "response_id"), relationship = "one-to-one")
    }
    enhancer_gene <- dplyr::ungroup(dplyr::select(dplyr::mutate(dplyr::group_by(enhancer_gene, 
        grna_target, response_id), test_stat_distribution = compute_distribution_teststat(num_trt_cells = num_trt_cells, 
        num_cntrl_cells = num_cntrl_cells, num_trt_cells_sq = num_trt_cells_sq, 
        expression_mean = expression_mean, expression_size = expression_size, 
        fold_change_mean = fold_change_mean, fold_change_sd = fold_change_sd), 
        mean_test_stat = unlist(test_stat_distribution)["mean"], 
        sd_test_stat = unlist(test_stat_distribution)["sd"], 
        QC_prob = compute_QC(fold_change_mean = fold_change_mean, 
            expression_mean = expression_mean, expression_size = expression_size, 
            num_cntrl_cells = num_cntrl_cells, num_trt_cells = num_trt_cells, 
            n_nonzero_trt_thresh = n_nonzero_trt_thresh, n_nonzero_cntrl_thresh = n_nonzero_cntrl_thresh)), 
        -test_stat_distribution))
    if (is.null(cutoff)) {
        cutoff <- dplyr::pull(dplyr::select(dplyr::summarize(enhancer_gene, 
            cutoff = adjusted_cutoff(mean_list = mean_test_stat, 
                sd_list = sd_test_stat, multiple_testing_alpha = multiple_testing_alpha, 
                multiple_testing_method = multiple_testing_method, 
                side = side, QC_prob = QC_prob)), cutoff))
    }
    enhancer_gene <- dplyr::mutate(enhancer_gene, cutoff = cutoff, 
        power = rejection_computation(mean_list = mean_test_stat, 
            sd_list = sd_test_stat, side = side, cutoff = cutoff) * 
            (1 - QC_prob))
    output <- list(individual_power = dplyr::select(enhancer_gene, 
        grna_target, response_id, power), expected_num_discoveries = sum(enhancer_gene$power))
    return(output)
}