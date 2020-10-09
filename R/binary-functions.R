#' Tests for binary outcomes
#'
#' @description performs a test for testing the null hypothesis that the proportions (probabilities of success) in two groups are the same.
#'
#' @param binary 0=nonresponse, 1=response.
#' @param treat The treatment-group indicator, normally 0=control, 1=intervention.
#' @param var_est variance estimator: Unpooled/Pooled
#'
#' @export
#'
#' @return List: standardized statistic, statistic and variance.
#' @author Marta Bofill Roig

bintest <- function(binary,treat,var_est='Unpooled'){

  db=cbind.data.frame(binary=binary, treat=treat)
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
  if(var_est=='Unpooled'){
    var.bin= phat_group0*(1-phat_group0)/n0+phat_group1*(1-phat_group1)/n1
  }else{
    var.bin = (n/(n0*n1))*phat_pooled*(1-phat_pooled)
  }
  Zb = (phat_group1-phat_group0)/sqrt(var.bin)

  # return(list=c(Test=Zb,Ub=sqrt((n0*n1)/n)*(phat_group1-phat_group0),sd=sqrt(phat_pooled*(1-phat_pooled))  ))
  return(list=c(Test=Zb,Ub=sqrt((n0*n1)/n)*(phat_group1-phat_group0),sd=sqrt(var.bin/(n/(n0*n1)))  ))

}



