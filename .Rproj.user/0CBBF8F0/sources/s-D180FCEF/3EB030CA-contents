#' Apply Simple Linear Regression
#'
#' @param x x-coordinate
#' @param y y-coordinate
#'
#' @return beta0,beta1,etc
#' @export
#'
#' @examples
#' library(slr)
#' slr(iris$Sepal.Length,iris$Sepal.Width)
slr<-function(x,y){
    n<-length(x)
    split.screen(c(1,3))
    screen(1)
    plot(x,y,pch=16)
    title(main="data scatterplot")

    meanx<-mean(x)
    meany<-mean(y)
    lxx<-sum((x-meanx)^2)
    lyy<-sum((y-meany)^2)
    lxy<-sum((x-meanx)*(y-meany))

    beta_1<-lxy/lxx#beta_1
    beta_0<-meany-beta_1*meanx#beta_0

    screen(2)
    plot(x,y,pch=16)
    points(x,beta_0+beta_1*x,type="l")
    title(main="regression plot")

    y_hat<-beta_0+beta_1*x #y hat
    sse<-sum((y_hat-y)^2)
    ssr<-sum((y_hat-meany)^2)
    sst<-ssr+sse
    sigma_hat<-sqrt(1/(n-2)*sse) #sigma hat

    alpha<-0.05

    sd.beta_0<-sqrt((1/n+(meanx^2)/lxx))*sigma_hat
    sd.beta_1<-sqrt(sigma_hat^2/lxx)
    beta_1l<-beta_1-qt(1-alpha/2,n-2)*sd.beta_1
    beta_1u<-beta_1+qt(1-alpha/2,n-2)*sd.beta_1
    beta_0l<-beta_0-qt(1-alpha/2,n-2)*sd.beta_0
    beta_0u<-beta_0+qt(1-alpha/2,n-2)*sd.beta_0

    R<-ssr/sst

    f<-(ssr/1)/(sse/(n-2))#f test
    p1<-pf(f,1,n-2)

    t1<-beta_1/sd.beta_1#t test
    p2<-pt(t1,n-2)

    r<-lxy/(sqrt(lxx*lyy))
    t2<-sqrt(n-2)*r/sqrt(1-r^2);
    p3<-pt(t2,n-2,lower.tail=TRUE)#p value

    screen(3)
    e<-y_hat-y
    n<-length(e)
    sigma_u<-seq(2*sigma_hat,2*sigma_hat,length.out=n)
    sigma_l<-seq(-2*sigma_hat,-2*sigma_hat,length.out=n)
    plot(x,e,pch=16,ylim=c(5,-5))
    points(x,sigma_u,type="l")
    points(x,sigma_l,type="l")
    title(main="residual plot")

    cat("Beta 0: ",beta_0,
        "\nBeta 1: ",beta_1,
        "\nr square:",R,
        "\nF-statistic: ",f,"p-value: ",1-p1,
        "\nt-statistic: ",t1,"p-value: ",2*p2)
}
