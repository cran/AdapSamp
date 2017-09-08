#' Concave-Convex Adaptive Rejection Sampling Algorithm
#' rCCARS generates a sequence by Concave- Convex Adaptive Rejection Sampling algorithm from target densities with bounded domain. For unbounded domain's rCCARS can also be used for sampling approximately, especially for densities whose probability is too small in tails.
#' @param n Desired sample size
#' @param cvformula,ccformula Convex and concave decompositions for -ln(p(x)) where p(x) is the kernal of target density
#' @param min,max Domain except positive and negative infinity
#' @param sp Supporting set
#' @details Strictly speaking,Concave-Convex Adaptive Rejection Sampling Algorithm can generate samples from target densities who have bounded domains. However, if the target distribution has a tiny probability in the tails, we can also use rCCARS for sampling approximately. For example, if we want to get a sequence drawn from N(0,1), and we know that
#' @references Teh Y W. Concave-Convex Adaptive Rejection Sampling[J]. Journal of Computational & Graphical Statistics, 2011, 20(3):670-691.
#' @examples
#' #Generalized inverse bounded gaussian distribution with lambda=-1 and a=b=2
#' x<-rCCARS(100,"x+x^-1","2*log(x)",0.001,100,1)
#' hist(x,breaks=20,probability = TRUE);lines(density(x,bw=0.1),col="red",lwd=2,lty=2)
#' f <- function(x) {x^(-2)*exp(-x-x^(-1))/0.2797318}
#' lines(seq(0,5,0.01),f(seq(0,5,0.01)),lwd=2,lty=3,col="blue")
#'
#' #Expontional bounded distribution
#' x<-rCCARS(100,"x^4","-8*x^2+16",-3,4,c(-2,1))
#' hist(x,breaks=30,probability=TRUE);lines(density(x,bw=0.05),col="blue",lwd=2,lty=2)
#' f <- function(x) exp(-(x^2-4)^2)/ 0.8974381
#' lines(seq(-3,4,0.01),f(seq(-3,4,0.01)),col="red",lty=3,lwd=2)
#'
#' #Makeham bounded distribution
#' x<-rCCARS(100,"x+1/log(2)*(2^x-1)","-log(1+2^x)",0,5,c(1,2,3))
#' hist(x,breaks=30,probability=TRUE);lines(density(x,bw=0.05),col="blue",lwd=2,lty=2)
#' f <- function(x){(1+2^x)*exp(-x-1/log(2)*(2^x-1))}
#' plot(seq(0,5,0.01),f(seq(0,5,0.01)),col="red",lty=3,lwd=2,type="l")
#' @export
rCCARS <- function(n,cvformula,ccformula,min,max,sp){
p <- function(x){eval(parse(text=paste("exp(-(",cvformula,")-(",ccformula,"))",sp="")))}
x_final <- numeric(n)
for( k in 1:n){
  support <- sp
  xrange <- c(min,max)
  u=0
  prop=-1
  while(u>prop){
    allpt <- sort(c(xrange,support))
    convex <- function(x){eval(parse(text=cvformula))}
    drv1or<- function(x){eval(stats::D(parse(text=cvformula),"x"))}
    der <- drv1or(allpt)
    crossx <- c()
    crossy <- c()
    for(i in 1:(length(allpt)-1)){
      A <- matrix(c(der[i],-1,der[i+1],-1),nrow=2,byrow=1)
      b <- c(der[i]*allpt[i]-convex(allpt)[i],der[i+1]*allpt[i+1]-convex(allpt)[i+1])
      crossx[i] <- solve(A,b)[1]
      crossy[i] <- solve(A,b)[2]
    }
    rubbish1 <- data.frame(X=c(crossx,allpt),Y=c(crossy,convex(allpt)))
    xconvex <- c(rubbish1[order(rubbish1$X),][1])$X
    yconvex <- c(rubbish1[order(rubbish1$X),][2])$Y
 tan1<- numeric(length(xconvex)-1)
    int1<- numeric(length(xconvex)-1)
    for (i in 1:length(tan1)){
      tan1[i] <- (yconvex[i+1]-yconvex[i])/(xconvex[i+1]-xconvex[i])
      int1[i] <- (yconvex[i+1]-yconvex[i])/(xconvex[i+1]-xconvex[i])*(-xconvex[i])+yconvex[i]
    }
    concave <- function(x){eval(parse(text=ccformula))}
    tan2 <- numeric(length(allpt)-1)
    for(i in 1:length(tan2)){
      tan2[i] <- (concave(allpt[i+1])-concave(allpt[i]))/(allpt[i+1]-allpt[i])
    }
    int2<- numeric(length(allpt)-1)
    for(i in 1:length(tan2)){
      int2[i] <- -tan2[i]*allpt[i]+concave(allpt[i])
    }
    xconcave <- rep(allpt[1:length(allpt)-1],rep(2,length(allpt)-1))
    yconcave <- concave(xconcave)
    tan2 <- rep(tan2,rep(2,length(tan2)))
    int2 <- rep(int2,rep(2,length(int2)))
    IntSum <- numeric(length(tan2))
    for(i in 1:length(tan2)){
      fun <- function(x){
        exp(-(tan1[i]+tan2[i])*x-int1[i]-int2[i])
      }
      IntSum[i] <- stats::integrate(fun,xconvex[i],xconvex[i+1])[[1]]
    }
    cum=c(0,cumsum(IntSum/sum(IntSum)))
    rdm <- stats::runif(1)
    idx <- which(rdm<cumsum(IntSum/sum(IntSum)))[1]
    x_star <- log(((rdm-cum[idx])*(-tan1[idx]-tan2[idx])*sum(IntSum))*exp(int1[idx]+int2[idx])+exp(-(tan1[idx]+tan2[idx])*xconvex[idx]))/(-tan2[idx]-tan1[idx])
    u <- stats::runif(1)
    prop=p(x_star)/exp(-(tan1[idx]+tan2[idx])*x_star-int1[idx]-int2[idx])
    support <- sort(c(x_star,support))
  }
  x_final[k]=x_star
}
x_final
}

