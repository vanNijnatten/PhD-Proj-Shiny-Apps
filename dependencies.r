# Check installation, install and load dependencies
loadPackage <- function(lib, bioc = FALSE, gh = FALSE, gh.path = "") {
  if (!require(lib, character.only = TRUE, quietly = TRUE)) {
    if (!bioc & !gh) {
      do.call(install.packages, list(pkgs = c(lib), character.only = TRUE))
    } else if (bioc) {
      BiocManager::install(lib)
    } else if (gh) {
      remotes::install_github(gh.path)
    }
    do.call(install.packages, list(lib))
  }
}

#
loadPackage("BiocManager")
loadPackage("DT")
loadPackage("ggpubr")
loadPackage("ggrepel")
loadPackage("markdown")
loadPackage("patchwork")
loadPackage("pheatmap")
loadPackage("RColorBrewer")
loadPackage("shiny")
loadPackage("tidyverse")

loadPackage("biomaRt", bioc = TRUE)
loadPackage("edgeR", bioc = TRUE)
loadPackage("enrichplot", bioc = TRUE)
