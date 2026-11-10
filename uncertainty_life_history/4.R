
# Source all mean IPM scripts -------------------------------------------------

#"adler_2007_ks"
#"anderson_2016_az"
#"anderson_2016_mt"
#"christensen_2021_nm"
#"chu_2013_co"
#"moore_2021_az"
#"zachmann_2016_id"

dir_ipm <- "christensen_2021_nm"

ipm_files <- list.files(
  path = dir_ipm,
  pattern = "_ipm_mean\\.R$",
  recursive = TRUE, full.names = TRUE)

 for (file in ipm_files) {
  message("Sourcing: ", file)
  source(file)
}
