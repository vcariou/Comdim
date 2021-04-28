#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
.onUnload <- function (libpath) {
  library.dynam.unload("comdim", libpath)
}
NULL
