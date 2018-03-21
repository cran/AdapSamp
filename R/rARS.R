#' Adaptive Rejection Sampling Algorithm
#'
#' rARS generates a sequence of random numbers using the adaptive rejection sampling algorithm.
#'
#' @param n Desired sample size;
#' @param formula Kernal of the target density;
#' @param min,max Domain including positive and negative infinity of the target distribution;
#' @param sp Supporting set.
#' @export
#' @author Dong Zhang <\url{dzhang0716@126.com}>
#'
#' @examples
#'
#' # Example 1: Standard normal distribution
#' x1 <- rARS(100,"exp(-x^2/2)",-Inf,Inf,c(-2,2))
#'
#' # Example 2: Truncated normal distribution
#' x2 <- rARS(100,"exp(-x^2/2)",-2.1,2.1,c(-2,2))
#'
#' # Example 3: Normal distribution with mean=2 and sd=2
#' x3 <- rARS(100,"exp(-(x-2)^2/(2*4))",-Inf,Inf,c(-3,3))
#'
#' # Example 4: Exponential distribution with rate=3
#' x4 <- rARS(100,"exp(-3*x)",0,Inf,c(2,3,100))
#'
#' # Example 5: Beta distribution with alpha=3 and beta=4
#' x5 <- rARS(100,"x^2*(1-x)^3",0,1,c(0.4,0.6))
#'
#' # Example 6: Gamma distribution with alpha=5 and lambda=2
#' x6 <- rARS(100,"x^(5-1)*exp(-2*x)",0,Inf,c(1,10))
#'
#' # Example 7: Student distribution with df=10
#' x7 <- rARS(100,"(1+x^2/10)^(-(10+1)/2)",-Inf,Inf,c(-10,2))
#'
#' # Example 8: F distribution with m=10 and n=5
#' x8 <- rARS(100,"x^(10/2-1)/(1+10/5*x)^(15/2)",0,Inf,c(3,10))
#'
#' # Example 9:Cauchy distribution
#' x9 <- rARS(100,"1/(1+(x-1)^2)",-Inf,Inf,c(-2,2,10))
#'
#' # Example 10:Rayleigh distribution with lambda=1
#' x10 <- rARS(100,"2*x*exp(-x^2)",0,Inf,c(0.01,10))
#'
rARS <-
  function(n,formula,min=-Inf,max=Inf,sp){
  sp <- sort(sp)
  if(!is.character(formula)) stop("Unsuitable density function.")
  if (n<=0) stop("Unsuitable sample size.")
  if(min >= max) stop("Unsuitable domain.")
  p <- function(x){eval(parse(text=formula))}
  V <- function(x){-log(p(x))}
  x_final <- numeric(n)
  for(j in 1:n){
    Support <- sp
    if (!identical(Support,sort(Support))) stop("Put the supporting points in ascending order.")
    u=0
    compareprop=-1
    while(u>compareprop){
      tangent <- fderiv(V,Support,1)
      crosspoint=numeric(length(Support)+1)
      crosspoint[1]=min
      crosspoint[length(crosspoint)]=max
      crossvalue=numeric(length(Support)-1)
      for( i in 1:(length(Support)-1)){
        A=matrix(c(tangent[i],-1,tangent[i+1],-1),nrow=2,byrow=T)
        b=c(tangent[i]*Support[i]-V(Support)[i],tangent[i+1]*Support[i+1]-V(Support)[i+1])
        solve(A,b)
        crosspoint[i+1]=solve(A,b)[1]
        crossvalue[i]=solve(A,b)[2]
      }
       IntSum <- numeric(length(Support))
      for (i in 1:length(IntSum)){
        expfun=function(x){
          exp(-tangent[i]*(x-Support[i])-V(Support)[i])
        }
        IntSum[i]= integrate(expfun,crosspoint[i],crosspoint[i+1])[[1]]
      }
      rdm <-  runif(1)
      cum=c(0, cumsum(IntSum/sum(IntSum)))
      idx <- which(rdm<cumsum(IntSum/sum(IntSum)))[1]
      x_star <-
        log((rdm-cum[idx]+exp(tangent[idx]*Support[idx]-V(Support)[idx])*
               exp(-tangent[idx]*crosspoint[idx])/sum(IntSum)/(-tangent[idx]))*
              sum(IntSum)*(-tangent[idx])/exp(tangent[idx]*Support[idx]-V(Support)[idx]))/(-tangent[idx])
      u <-  runif(1)
      compareprop <- p(x_star)/exp(-tangent[idx]*(x_star-Support[idx])-V(Support)[idx])
      Support <- sort(c(Support,x_star))
    }
    x_final[j]=x_star
  }
  x_final
}
