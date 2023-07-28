install.packages("data.table",repos="http://cloud.r-project.org")
install.packages("xlsx",repos="http://cloud.r-project.org")
install.packages("XLConnect",repos="http://cloud.r-project.org")
install.packages("jsonlite",repos = "http://cloud.r-project.org/")
install.packages("stringr",repos = "http://cloud.r-project.org/")
install.packages("stringdist",repos = "http://cloud.r-project.org/")
install.packages("dplyr",repos = "http://cloud.r-project.org/")
install.packages("RPostgreSQL",repos = "http://cloud.r-project.org/")

BiocManager::install()
BiocManager::install("ropls")
BiocManager::install("graph")

install.packages("NMF",repos = "http://cloud.r-project.org/")