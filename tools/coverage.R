library(covr)

x <- covr::package_coverage(
  path = "~/syncomBuildR",
  type = c("examples", "tests", "vignettes"),
  combine_types = TRUE,
  quiet = TRUE,
  clean = TRUE,
  pre_clean = TRUE#,
#  line_exclusions = list("R/file.R")
)
x

report(x, file = "~/syncomBuildR/tools/scb-coverage-report.html")