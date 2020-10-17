#
setwd("C:/School/UVA/Research/Vignette")
library(rmarkdown)
library(tinytex)

#bioconductor style:
library(BiocStyle)
rmarkdown::render("adipogenesis.Rmd", BiocStyle::pdf_document(keep_tex=TRUE))
rmarkdown::render("adipogenesis.Rmd", BiocStyle::html_document())


#default style
#do NOT run library(BiocStyle) if you want to use default style
rmarkdown::render("adipogenesis.Rmd", pdf_document())
rmarkdown::render("adipogenesis.Rmd", html_document())
