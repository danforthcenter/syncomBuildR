library(lintr)
library(styler)
devtools::load_all("~/syncomBuildR")

x <- lintr::lint_package(path = "~/syncomBuildR/",
                    linters = linters_with_defaults(line_length_linter(length = 105L),
                                                    object_name_linter(styles = c("snake_case", "symbols",
                                                                                  "camelCase", "dotted.case",
                                                                                  "lowercase", "UPPERCASE")),
                                                    brace_linter(allow_single_line = TRUE),
                                                    return_linter = NULL
                    ))
length(x)
names(x)
x
#* dry run styling
style_pkg("~/syncomBuildR", dry = "off", scope = "line_breaks")
#* STOP HERE AND CHECK

style_pkg("~/syncomBuildR", dry = "off", scope = "tokens")

if(FALSE){
  c("R/plot_scbnet.R", "R/plot_thresh.R", "R/pullNode.R", "R/scbnet_class.R", 
    "R/thresh.R", "R/tresh_class.R")
  file = "vignettes/dada2.Rmd"
  styler::style_file(file, scope = "line_breaks")
  lintr::lint(file, linters = lintr::linters_with_defaults(lintr::line_length_linter(length = 105L),
                                                           lintr::object_name_linter(styles = c("snake_case", "symbols",
                                                                                                "camelCase", "dotted.case",
                                                                                                "lowercase", "UPPERCASE")),
                                                           lintr::brace_linter(allow_single_line = TRUE)
  ))
  #* or for tokens
  styler::style_file(file, dry = "off", scope = "tokens")
}
