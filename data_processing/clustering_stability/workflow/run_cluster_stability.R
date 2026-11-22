args <- commandArgs(trailingOnly = TRUE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
workflow_root <- dirname(script_path)

config <- file.path(workflow_root, "config", "config.yaml")
stage <- "all"
forward <- character()
for (arg in args) {
  if (!startsWith(arg, "--") && identical(config, file.path(workflow_root, "config", "config.yaml"))) {
    config <- normalizePath(arg, mustWork = TRUE)
  } else if (startsWith(arg, "--stage=")) {
    stage <- sub("^--stage=", "", arg)
  } else {
    forward <- c(forward, arg)
  }
}
valid <- c("scan", "pilot", "formal", "summarize", "marker-rediscovery", "all")
if (!stage %in% valid) stop("Unknown stage: ", stage, call. = FALSE)

run_script <- function(name, arguments = character()) {
  command <- file.path(workflow_root, "scripts", name)
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    args = c(command, paste0("--config=", config), arguments)
  )
  if (!identical(status, 0L)) stop("Stage failed: ", name, call. = FALSE)
}

if (stage == "scan") run_script("01_scan_resolution_grid.R", forward)
if (stage == "pilot") run_script("02_run_resampling.R", c("--mode=pilot", forward))
if (stage == "formal") run_script("02_run_resampling.R", c("--mode=formal", forward))
if (stage == "summarize") run_script("03_summarize_stability_range.R", forward)
if (stage == "marker-rediscovery") {
  run_script("02_run_resampling.R", c("--mode=formal", "--marker-mode=rediscovery", forward))
}
if (stage == "all") {
  run_script("01_scan_resolution_grid.R")
  run_script("02_run_resampling.R", c("--mode=pilot", forward))
  summary_forward <- forward[grepl("^--(?:count-mode|resolutions)=", forward)]
  run_script("03_summarize_stability_range.R", c("--mode=pilot", summary_forward))
}
