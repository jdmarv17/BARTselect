#' Function to get VIP from filtered tree data frame
#'
#' @param filtered.trees 
#'
#' @return A data frame of VIPs for each variable
#' 
#'
#' @examples
#' library(MASS)
#' library(dbarts)
#' library(dplyr)
#' library(janitor)
#' data(birthwt)
#' data = birthwt %>% dplyr::select(., -low)
#' formula = bwt ~ .
#' dbarts_data = dbartsData(formula, data)
#' dbarts_x = clean_names(dbarts_data@x)
#' dbarts_y = dbarts_data@y
#' vars = colnames(dbarts_x)
#' bart_run = run_bart(formula, data, vars = vars)
#' 
#' vip = get_vip(filtered_trees = bart_run[[4]])
#' 
get_vip = function(filtered_trees) {
  
  # get vips
  num.splits = nrow(filtered_trees)
  vip = 
    as.data.frame(filtered_trees) %>% 
    group_by(var) %>% 
    summarise(prop.split = n() / num.splits)
  
  return(vip)
}