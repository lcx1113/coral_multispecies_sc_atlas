# Shared configuration for Adig, Ahya, Amur, Fusp, Gfas, Hcoe, Pdam and Plut.
# Pcru differences are documented in ../README.md; no separate Pcru config is shipped.
qc_config <- list(
  random_seed = 100L,

  emptydrops = list(
    fdr_cutoff = 0.001,
    by_rank = 40000L,
    rank_multiplier = NA_real_,
    lower = 100L,
    retain = Inf,  # Descriptive only: stage function retains its source hard-coded Inf.
    niters = 100000L,
    test_ambient = TRUE
  ),

  soupx = list(
    min_cells = 5L,
    cluster_resolution = 0.6,
    temp_nfeatures = 2000L,
    temp_npcs = 20L,
    temp_neighbor_dims_max = 10L,
    round_to_int = TRUE,
    force_accept = "default",
    save_autoEstCont_plot = TRUE
  ),

  marker_filter = list(
    assay = "RNA",
    random_seed = 123L,
    sct_vst_flavor = "v2",  # Use SCT v2.
    cutoff = 15L,
    rna_resolution = 0.6,
    sct_leiden_resolution = 4,
    rna_nfeatures = 2000L,
    marker_padj_cutoff = 0.05,
    marker_log2fc_cutoff = 0.50,
    marker_pct_in_cutoff = 0.25,
    marker_delta_pct_cutoff = 0.15,
    max_cells_per_ident = 500L,
    heatmap_top_n = 10L
  ),

  top_mito = list(
    top_umi_quantile = 0.99,
    mito_cutoff = 1
  ),

  doubletfinder = list(
    doublet_rate = 0.01,
    pN = 0.25,
    dims_max = 20L,
    neighbor_dims_max = 10L,
    resolution = 0.6,
    use_homotypic_adjustment = TRUE
  ),

  final_recluster = list(
    workflow = "RNA",
    nfeatures = 2000L,
    npcs_max = 30L,
    dims_max = 20L,
    resolution = 0.6,
    marker_padj_cutoff = 0.05,
    marker_logfc_cutoff = 0.25,
    marker_min_pct = 0.10,
    heatmap_top_n = 10L
  ),

  output = list(
    save_soupx_object = TRUE,
    save_before_doubletfinder_object = TRUE,
    save_final_object = TRUE,
    save_raw_marker_gz = TRUE,
    overwrite_existing = FALSE,
    resume_from_existing = TRUE
  )
)
