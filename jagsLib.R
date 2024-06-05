require(R2jags)
require(infinitefactor)
require(abind)

condMan.mod="
model{
  for (j in 1:J){
    pdel[j] ~ dgamma(.5,.5*rDel^2)
    mu[j] ~ dnorm(a,1/b^2)
    for (d in 1:D){
      lambda[j,d] ~ dnorm(0,1/(rLam^2))}}
  for (i in 1:I){
    for (d in 1:D){
      phi[i,d]~dnorm(0,1)}}
  for (i in 1:I){
    for (j in 1:J){
      y[i,j] ~ dnorm(mu[j]+sum(phi[i,]*lambda[j,]),pdel[j])}}
}
"

runCondMan=function(y,D=2,prior,M=2000){
  I=dim(y)[1]
  J=dim(y)[2]
  setup=list("I"=I,"J"=J,"D"=D,"y"=y)
  out=jags(data=c(setup,prior), 
           parameters=c("lambda","phi","pdel","mu"), 
           model.file = textConnection(condMan.mod), 
           n.chains=1,n.iter=M,n.burnin=M/10,n.thin=1)
  return(out)
}
#runCondMan(y)


condManL.mod="
model{
  for (j in 1:J){
    pdel[j] ~ dgamma(.5,.5*rDel^2)
    mu[j] ~ dnorm(a,1/b^2)
    for (d in 1:D){
      lambda[j,d] ~ dnorm(0, ifelse(d>j,sv^-2,rLam^-2))T(ifelse(d==j,0,-100),)}}
  for (i in 1:I){
    for (d in 1:D){
      phi[i,d]~dnorm(0,1)}}
  for (i in 1:I){
    for (j in 1:J){
      y[i,j] ~ dnorm(mu[j]+sum(phi[i,]*lambda[j,]),pdel[j])}}
}
"

runCondManL=function(y,D=2,prior,M=2000){
  I=dim(y)[1]
  J=dim(y)[2]
  setup=list("I"=I,"J"=J,"D"=D,"y"=y,"sv"=.00000001)
  out=jags(data=c(setup,prior), 
           parameters=c("lambda","phi","pdel","mu"), 
           model.file = textConnection(condManL.mod), 
           n.chains=1,n.iter=M,n.burnin=M/10,n.thin=1)
  return(out)
}
