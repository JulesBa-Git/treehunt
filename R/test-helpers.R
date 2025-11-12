#' @keyword internal
#' @noRd
make_test_tree_ctr1 <- function(data, depth, upper_bound = NULL, name = NULL) {
  .test_create_tree_constructor1(data, depth, upper_bound, name)
}

#' @keyword internal
#' @noRd
make_test_tree_ctr2 <- function(depth) {
  .test_create_tree_constructor2(depth)
}

#' @keywords internal
#' @noRd
inspect_tree <- function(tree) {
  .test_tree_info(tree)
}