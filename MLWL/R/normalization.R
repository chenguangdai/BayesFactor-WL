#' Normalizing two log constants to sum 1
#'
#' @param target_logZ        Unnormalized log normalizing constants of the target distribution
#' @param surrogate_logZ     Unnormalized log normalizing constants of the surrogate distribution
#' @return Two normalized log constants
#' @export
normalization <- function(target_logZ, surrogate_logZ){
  maxlogZ <- max(target_logZ, surrogate_logZ)
  logsumZ <- log(sum(exp(target_logZ - maxlogZ) + exp(surrogate_logZ - maxlogZ)))
  target_logZ <- target_logZ - maxlogZ - logsumZ
  surrogate_logZ <- surrogate_logZ - maxlogZ - logsumZ
  return(list(target_logZ = target_logZ, surrogate_logZ = surrogate_logZ))
}