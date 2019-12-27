
##################################################################################
# FUNCTION bintest
# Computes the proportions' test
# Marta Bofill Roig
# 
# First version 2019-05-09
##################################################################################

##################################################################################
# ARGUMENTS:
#   * binary outcome 
#   * group indicator

bintest <- function(x1,x2,x3='Unpooled'){
  # x1=data$binary; x2=data$treat;
  db=cbind.data.frame(binary=x1, treat=x2)
  
  db1=subset(db,db$treat==1)
  db0=subset(db,db$treat==0)
  
  n1=dim(db1)[1]
  n0=dim(db0)[1]
  n=n0+n1 
  
  # BINARY TEST
  ###################################### 
  phat_group0 = sum(subset(db, db$treat==0)$binary)/n0
  phat_group1 = sum(subset(db, db$treat==1)$binary)/n1 
  phat_pooled = (phat_group0*n0 + phat_group1*n1)/n  
  if(x3=='Unpooled'){
    var.bin= phat_group0*(1-phat_group0)/n0+phat_group1*(1-phat_group1)/n1
  }else{
    # phat_pooled = (phat_group0*n0 + phat_group1*n1)/n  
    var.bin = (n/(n0*n1))*phat_pooled*(1-phat_pooled)
  }
  # test risk difference with pooled variance /unpooled v
  # Zb = sqrt(ss/2)*(phat_group1-phat_group0)/sqrt(phat_pooled*(1-phat_pooled))  
  Zb = (phat_group1-phat_group0)/sqrt(var.bin) 
  
  
  # return(list=c(Test=Zb,Ub=sqrt((n0*n1)/n)*(phat_group1-phat_group0),sd=sqrt(var.bin) ))
  return(list=c(Test=Zb,Ub=sqrt((n0*n1)/n)*(phat_group1-phat_group0),sd=sqrt(phat_pooled*(1-phat_pooled))  ))
  # return(list=c(Test=Zb,d=phat_group1-phat_group0,sd=sqrt(var.bin) ))
  # return(Zb)
} 


##################################################################################
# Example
# bintest(data$binary, data$treat)





