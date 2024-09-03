wd <- getwd()
library(ggplot2)
setwd("~/syncomBuildR/R")
tdf <- do.call(rbind, lapply(dir("~/syncomBuildR/man", full.names = TRUE), function(doc) {
  message(paste0("Running ", gsub(".*/", "", doc), " examples"))
  x <- system.time({
    pkgload::run_example(doc, quiet = TRUE)
  })
  data.frame(time = x[["elapsed"]], fun = gsub(".*/", "", doc))
}))
setwd(wd)

ggplot(tdf, aes(x = 1, y = time, fill = fun)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d()+
  theme_light()
