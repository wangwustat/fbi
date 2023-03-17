#' @title Double-Project Imputation of Missing Value in Panel Data
#'
#' @description
#' \code{tp_apc} imputes the missing values in a given panel data using the
#' method of "Double-Project".
#'
#' @import stats
#' @export
#'
#' @param X a matrix of size T by N with missing values.
#' @param kmax integer, indicating the maximum number of factors.
#' @param center logical, indicating whether or not X should be demeaned
#' @param standardize logical, indicating whether or not X should be scaled.
#' @param re_estimate logical, indicating whether or not output factors,
#' `Fhat`, `Lamhat`, `Dhat`, and `Chat`, should be re-estimated from the imputed data.
#'
#' @return a list of elements:
#' \item{Fhat}{estimated F}
#' \item{Lamhat}{estimated Lambda}
#' \item{Dhat}{estimated D}
#' \item{Chat}{euqals Fhat x Lamhat'}
#' \item{ehat}{equals Xhat - Chat}
#' \item{data}{X with missing data imputed}
#' \item{X}{the original data with missing values}
#' \item{kmax}{the maximum number of factors}
#' \item{center}{logical, indicating whether or not X was demeaned in the algorithm}
#' \item{standardize}{logical, indicating whether or not X was scaled in the algorithm}
#' \item{re_estimate}{logical, indicating whether or not output factors,
#' `Fhat`, `Lamhat`, `Dhat`, and `Chat`, were re-estimated from the imputed data}

dp_apc <- function(X, kmax, center = FALSE, standardize = FALSE,
                   re_estimate = TRUE) {


  missing <- is.na(X) # missing patters
  goodT <- rowSums(is.na(X)) == 0
  goodN <- colSums(is.na(X)) == 0
  T1 <- sum(goodT)
  N1 <- sum(goodN)

  # The first step tp (Tall-projection)
  res_tp = tp_apc(X, kmax, center , standardize , re_estimate = FALSE)

  # The second step: Wide-projection
  res_wp = tp_apc(t(X), kmax, center , standardize , re_estimate = FALSE)

  # combine the two "by averaging"
  # There might be more efficient ways for combing
  Chat = (res_tp$Chat + t(res_wp$Chat))/2

  Xhat <- X   # estimated data
  Xhat[missing] <- Chat[missing]

  out = list()
  out$Chat = Chat
  if (re_estimate){
    reest <- fbi::apc(Xhat, kmax)
    out$Fhat <- reest$Fhat
    out$Lamhat <- reest$Lamhat
    out$Dhat <- reest$Dhat
    out$Chat <- reest$Chat
    data <- X
    data[missing] <- out$Chat[missing]
    #out$data <- data

  }
  return(out)


}





