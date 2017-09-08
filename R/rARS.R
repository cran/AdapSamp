#' Adaptive Rejection Sampling Algorithm
#' rARS generates a sequence of random variables using ARS algorithm.
#' @param n Desired sample size
#' @param formula Kernal density of target log target density
#' @param min,max Domain including positive and negative infinity
#' @param sp Supporting set.
#' @author Zhangdong<\url{dzhang0716@126.com}>
#' @examples
#' ###Running the following codes my take you a few minutes!
#' 
#' #Standard normal distribution
#' x<-rARS(100,"exp(-x^2/2)",-Inf,Inf,c(-2,2))
#'
#' #Truncated normal distribution
#' rARS(100,"exp(-x^2/2)",-2.1,2.1,c(-2,2))
#'
#' #Normal distribution with mean=2 and sd=2
#' rARS(100,"exp(-(x-2)^2/(2*4))",-Inf,Inf,c(-3,3))
#'
#' #Exponential distribution with rate=3
#' rARS(100,"exp(-3*x)",0,Inf,c(2,3,100))
#'
#' #Beta distribution
#' rARS(100,"x^2*(1-x)^3",0,1,c(0.4,0.6))
#'
#' #Gamma distribution
#' rARS(100,"x^(5-1)*exp(-2*x)",0,Inf,c(1,10))
#'
#' #Student distribution
#' rARS(100,"(1+x^2/10)^(-(10+1)/2)",-Inf,Inf,c(-10,2))
#'
#' #F distribution
#' rARS(100,"x^(10/2-1)/(1+10/5*x)^(15/2)",0,Inf,c(3,10))
#'
#' #Cauchy distribution
#' rARS(100,"1/(1+(x-1)^2)",-Inf,Inf,c(-2,2,10))
#'
#' #Rayleigh distribution with lambda=1
#' rARS(100,"2*x*exp(-x^2)",0,Inf,c(0.01,10))
#' @export
rARS <-
  function(n,formula,min=-Inf,max=Inf,sp){
  sp <- sort(sp)
  if(!is.character(formula)) stop("Density function is inappropriate, please look up examples for help")
  if (n<=0) stop("Length of sequence shouble be large than 0")
  if(min >= max) stop("Minimum of domain shouble be larger than maximum")
  p <- function(x){eval(parse(text=formula))}
  V <- function(x){-log(p(x))}
  x_final <- numeric(n)
  for(j in 1:n){
    Support <- sp
    if (!identical(Support,sort(Support))) stop("Support points should be in ascending order")
    u=0
    compareprop=-1
    while(u>compareprop){
      tangent <- pracma::fderiv(V,Support,1)
      crosspoint=numeric(length(Support)+1)
      crosspoint[1]=min
      crosspoint[length(crosspoint)]=max
      crossvalue=numeric(length(Support)-1)
      for( i in 1:(length(Support)-1)){
        A=matrix(c(tangent[i],-1,tangent[i+1],-1),nrow=2,byrow=TRUE)
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
        IntSum[i]=stats::integrate(expfun,crosspoint[i],crosspoint[i+1])[[1]]
      }
      rdm <- stats::runif(1)
      cum=c(0, cumsum(IntSum/sum(IntSum)))
      idx <- which(rdm<cumsum(IntSum/sum(IntSum)))[1]
      x_star <-
        log((rdm-cum[idx]+exp(tangent[idx]*Support[idx]-V(Support)[idx])*
               exp(-tangent[idx]*crosspoint[idx])/sum(IntSum)/(-tangent[idx]))*
              sum(IntSum)*(-tangent[idx])/exp(tangent[idx]*Support[idx]-V(Support)[idx]))/(-tangent[idx])
      u <- stats::runif(1)
      compareprop <- p(x_star)/exp(-tangent[idx]*(x_star-Support[idx])-V(Support)[idx])
      Support <- sort(c(Support,x_star))
    }
    x_final[j]=x_star
  }
  x_final
}

