#' Not-in Operator (`%ni%`)
#'
#' A convenience operator to negate `%in%`, returning a logical vector indicating if there is **not** a match for each element of `x` in `vec`.
#'
#' @usage x \%ni\% vec
#'
#' @param x Vector or `NULL`; the values to be matched.
#' @param vec Vector or `NULL`; the values to be matched against.
#'
#' @return A logical vector of the same length as `x`, indicating which elements are **not** in `vec`.
#'
#' @examples
#' v1 = 1:3
#' v2 = 1:10
#' v1 %ni% v2
#' v2 %ni% v1
#'
#' @export
#' @name %ni%
#' @rdname ni_operator
`%ni%` = function(x, vec) {
  !(x %in% vec)
}
