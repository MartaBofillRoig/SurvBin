
# install.packages("boot")
library(boot) 


head(iris)


foo <- function(data, indices, cor.type){
  dt<-data[indices,]
  c(
    cor(dt[,1], dt[,2], method=cor.type),
    median(dt[,1]),
    median(dt[,2])
  )
}


set.seed(12345)
myBootstrap <- boot(iris, foo, R=1000, cor.type='s')
myBootstrap

head(myBootstrap)