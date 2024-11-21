if (requireNamespace("INLA", quietly = TRUE)) {
  # Make sure to avoid deprecated mesh functions; use fmesher instead:
  INLA::inla.setOption(fmesher.evolution.warn = TRUE)
  INLA::inla.setOption(fmesher.evolution.verbosity = "stop")
}
