# Print package version
pkg_info <- function(pkgs = pkgs) {
  x <- devtools::session_info(pkgs = pkgs)
  x <- as.data.frame(x$packages)
  x <- dplyr::filter(x, package %in% pkgs) %>% 
    dplyr::select("package", "loadedversion") %>% 
    dplyr::arrange(package)
  x
}