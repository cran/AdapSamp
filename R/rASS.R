#'Adaptive Slice Sampling Algorithm With Stepping-Out Procedures
#' rASS generates a sequence by Adaptive Slice Sampling Algorithm With Stepping-Out Procedures.
#' @param n Desired sample size
#' @param x0 Initial value
#' @param formula Target density function p(x)
#' @param w Length of the coverage interval
#' @references Neal R M. Slice sampling - Rejoinder[J]. Annals of Statistics, 2003, 31(3):758-767.
#' @examples
#' #Expotential distribution
#' x<-rASS(500,-1,"1.114283*exp(-(4-x^2)^2)",3)
#' f <- function(x){1.114283*exp(-(4-x^2)^2)}
#' plot(seq(-3,3,0.01),f(seq(-3,3,0.01)),lwd=2,lty=3,col="blue",type="l")
#' lines(density(x,bw=0.05),lwd=2,lty=2,col="red")
#'
#' #Mixed normal distribution
#' x <- rASS(500,2,"0.2/sqrt(2*pi)*exp(-x^2/2)+0.8/sqrt(2*pi*9)*exp(-(x-3)^2/2/9)",0.2)
#' f <- function(x){0.2/sqrt(2*pi)*exp(-x^2/2)+0.8/sqrt(2*pi*9)*exp(-(x-3)^2/2/9)}
#' plot(seq(-10,15,0.01),f(seq(-10,15,0.01)),lwd=2,lty=3,col="blue",type="l")
#' lines(density(x,bw=0.5),lwd=2,lty=2,col="red")
#'
#' #Mixed gamma distribution
#' x <- rASS(500,6,"0.3*2^8/gamma(8)*x^7*exp(-2*x)+0.7*5^4/gamma(4)*x^3*exp(-5*x)",0.3)
#' f <- function(x){0.3*2^8/gamma(8)*x^7*exp(-2*x)+0.7*5^4/gamma(4)*x^3*exp(-5*x)}
#' plot(seq(0,8,0.01),f(seq(0,8,0.01)),lwd=2,lty=3,col="blue",type="l")
#' lines(density(x,bw=0.2),lwd=2,lty=2,col="red")
#' @export
rASS <- function(n,x0=0,formula,w=3){
  f <- function(x){eval(parse(text=formula))}
  x_final=NULL
  Slice=NULL
  x_final[1]=x0
  for (i in 1:n){
    Slice[i]=stats::runif(1,0,f(x_final[i]))
    left=x_final[i]-stats::runif(1,0,w)
    right=left+w

    while(!((f(left)<Slice[i])&(f(right)<Slice[i]))){
      left=left-w
      right=right+w
    }

    x=stats::runif(1,left,right)
    while(f(x)<Slice[i]){
      if(x>x_final[i]) {right=x} else {left=x}
      x=stats::runif(1,left,right)
    }
    if(f(x)>Slice[i]){x_final[i+1]=x}
  }
  return(x_final)
}
