#' Not in operator
#'
#' A convenience operator to negate %in%.
#' 
#' @usage x %ni% vector
#' @param x Values to be checked.
#' @param vector Values to be compared against.
#' @return Logical vector indicating which elements of x are not in the vector
#' @export
#' 
`%ni%` = base::Negate(`%in%`)