#' Function to return a vector of indicators
#'
#' @param var_indices A vector of integers specifying position of 1's in the output vector- largest value no greater than p.
#' @param p Length of output indicator vector
#'
#' @return A vector of 1's and 0's. 1's are in positions specified by var_indices.
#' 
#' @keywords internal
#'
#' @examples
#' var_indices = c(1, 4, 5, 10)
#' p = 10
#' indicators = make_indicators(var_indices, p)
#' 
make_indicators = function(var_indices, p) {

  # get unique, sorted list and initialize the vector that will be returned
  ordered = unique(sort(var_indices))
  ind_vec = c()
  for (i in 1:length(ordered)) {
    if (i == 1) { # if its the first in the list
      next_vec = c(rep(0, times = (ordered[i] - 1)), 1)
      # if var_1 = 5, this makes a vector with four 0's and one 1
      
    } else { # otherwise need to repeat based on difference between last var
      next_vec = c(rep(0, times = (ordered[i] - ordered[i - 1] - 1)), 1)
      # this makes a vector with a series of 0's (possibly none), followed by a 1
    }
    
    ind_vec = c(ind_vec, next_vec) # bind vectors together to keep running list
    # this will have all the 1's we need after the loop is done
  }
  
  # now check if we need to append 0's to the end
  # if last var = p then we are done, otherwise append 0's
  if (ordered[length(ordered)] == p) {
    return(ind_vec)
    
  } else {
    to_append = p - ordered[length(ordered)]
    ind_vec = c(ind_vec, rep(0, times = to_append))
    return(ind_vec)
  }
  
}