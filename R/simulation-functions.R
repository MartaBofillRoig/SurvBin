#' Simulation binary and time-to-event data
#'
#' @description simulates two-sample binary and time-to-event dataset.
#'
#' @param a.shape shape parameter Weibull
#' @param b.scale scale parameter Weibull
#' @param rate.param Censoring distribution parameter (Rate for the exponential, Max for Uniform(0,max))
#' @param prob0 probability binary outcome
#' @param ass.par Association between binary and time-to-event outcome according to a Frank Copula
#' @param ss Sample size per arm
#' @param censoring Censoring distribution. Options: "Exp": Exponential; "Unif": Uniform
#'
#' @export
#'
#' @return Binary and time-to-event data
#' @references Trivedi et al, Appendix A
#' @author Marta Bofill Roig
#'

##################################################################################
simsurvbin <- function(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring="Exp"){

  # CENSORING
  ######################################
  if(censoring=="Exp"){
    TC1 = rexp(n= ss, rate = rate.param)
    TC0 = rexp(n= ss, rate = rate.param)
  }
  if(censoring=="Unif"){
    TC1 = runif(n= ss, min=0,max=rate.param)
    TC0 = runif(n= ss, min=0,max=rate.param)
  }


  # TREATMENT GROUP
  ######################################
  v1 <- runif(n=ss)
  v2 <- runif(n=ss)
  u1 <- v1
  u2 <- cop2(v1,v2,ass.par)

  v = cbind(u1,u2)

  # without latent variable
  BE1 = ifelse(v[,2]<prob0, 1, 0)

  # time-to-event (survival)
  TE1 = b.scale*(-log(1-v[,1]))^(1/a.shape)
  time1= ifelse(TE1<=TC1, TE1,TC1)
  status1 = ifelse(TE1<=TC1,1,0)
  treat1 = rep(1,ss)

  # CONTROL GROUP
  ######################################
  v1 <- runif(n=ss)
  v2 <- runif(n=ss)
  u1 <- v1
  u2 <- cop2(v1,v2,ass.par)

  v = cbind(u1,u2)

  # without latent variable
  BE0 = ifelse(v[,2]<prob0, 1, 0)

  # time-to-event (survival)
  TE0 = b.scale*(-log(1-v[,1]))^(1/a.shape)
  time0= ifelse(TE0<=TC0, TE0,TC0)
  status0 = ifelse(TE0<=TC0,1,0)
  treat0 = rep(0,ss)

  # TWO-SAMPLE db
  ######################################
  treat0=as.vector(treat0)
  treat1=as.vector(treat1)
  db = data.frame(binary=c(BE1, BE0), time=c(time1,time0), status=c(status1,status0),treat=c(treat1,treat0))


  return(db)
}

##################################################################################

# Auxiliar function
# cop2 : conditional transformation for copula FRANK
# Table A1 conditional transformation for copula FRANK
cop2 <- function(v1,v2,theta){
  u2 = -(1/theta)*log(1+(v2*(1-exp(-theta)))/(v2*(exp(-theta*v1)-1)-exp(-theta*v1)))
}
