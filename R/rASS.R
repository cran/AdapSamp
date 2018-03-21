#' Adaptive Slice Sampling Algorithm With Stepping-Out Procedures
#'
#' rASS generates a sequence of random numbers by the adaptive slice sampling algorithm with stepping-out procedures.
#' @param n Desired sample size;
#' @param x0 Initial value;
#' @param formula Target density function p(x);
#' @param w Length of the coverage interval.
#' @references Neal R M. Slice sampling - Rejoinder[J]. Annals of Statistics, 2003, 31(3):758-767.
#' @author Dong Zhang <\url{dzhang0716@126.com}>
#'
#' @export
#'
#' @examples
#'
#' # Example 1: Sampling from exponential distribution with bounded domain
#' x<-rASS(100,-1,"1.114283*exp(-(4-x^2)^2)",3)
#' plot(density(x))
#'
rASS <- function(n,x0=0,formula,w=3){
  f <- function(x){eval(parse(text=formula))}
  x_final=NULL
  Slice=NULL
  x_final[1]=x0
  for (i in 1:n){
    Slice[i]=runif(1,0,f(x_final[i]))
    left=x_final[i]-runif(1,0,w)
    right=left+w
    while(!((f(left)<Slice[i])&(f(right)<Slice[i]))){
      left=left-w
      right=right+w
    }
    x=runif(1,left,right)
    while(f(x)<Slice[i]){
      if(x>x_final[i]) {right=x} else {left=x}
      x=runif(1,left,right)
    }
    if(f(x)>Slice[i]){x_final[i+1]=x}
  }
  return(x_final)
}
