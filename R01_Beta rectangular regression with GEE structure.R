#=================================================================================================
# Code: Beta rectangular regression models based in GEE (Ribeiro, 2020)
# Author: Vinicius Osterne (vinicius@osterne -- www.osterne.com)
#=================================================================================================




#=================================================================================================
# Clearing memory
#=================================================================================================
rm(list=ls())






#=================================================================================================
# Required packages
#=================================================================================================
library("betareg")
library("matrixcalc")
library("MASS")
library("psych")
require("gridExtra")
library("betareg")
library("matrixcalc")
library("CVTuningCov")
library("KernSmooth")
library("tidyverse")






#=================================================================================================
# Functions of distribution  (dbetar, pbetar, rbetar, qbetar)
#=================================================================================================


#Beta density - Based in of Cribari-Neto and Ferrari (2004) --------------------------------------
dbetanova = function(x,mu,phi){
  a = mu * phi
  b = (1 - mu) * phi
  betanova = dbeta(x, shape1 = a, shape2 = b)
  return(betanova)
} 


#Density of rectangular beta ---------------------------------------------------------------------
dbetar = function(x,gama,phi,alpha){
  BR = c()
  for(i in 1:length(x)){
    nucleo = sqrt(1-4*alpha*gama*(1-gama))
    mu = (gama-0.5+0.5*nucleo)/(nucleo)
    BR[i] = 1 - nucleo + nucleo*dbetanova(x[i],mu,phi)   
  }
  return(BR)
} 



#Generating random numbers -----------------------------------------------------------------------
rbetar = function(n,gama,phi,alpha){
  nucleo = sqrt(1-4*alpha*gama*(1-gama))
  mu = (gama - (1/2) + (1/2)*(nucleo))/nucleo
  theta = 1 - nucleo  
  uni = runif(n)
  bet = rbetanova(n,mu,phi) #ajustar isso!
  unif = runif(n)
  y = (theta > unif)*uni + (theta < unif)*bet
  print(y)
}

#Rectangular beta cumulative distribution function ------------------------------------------------
pbetar = function(x,gama,phi,alpha){
  betar = function(t) dbetar(t,gama,phi,alpha)
  acum = integrate(betar,lower=10^(-100), upper = x)$value
  return(acum)
}




#Quantile function ---------------------------------------------------------------------------------
qbetar = function(q, gama, phi, alpha){
  library("GoFKernel")
  f = function(x) pbetar(x, gama, phi, alpha)
  f.inv = inverse(f,lower=0,upper=1)
  xq=c()
  aux=c()
  for (i in 1:length(q)) {
    xq[i] = f.inv(q[i])
    aux[i] = xq[i]
    if(aux[i]==1) {return(0.999)} else {
      if(aux[i]==0){return(0.001)} else{return(aux[i])}}
  }
  #return(xq)
}

















data = read.csv(choose.files(), header=T, sep=";", dec=",")
data = data.frame(cbind(data[5], log(data[2]), (log(data[2]))^2, data[3], data[4]))
head(data)

ind = data[,1]
y = data[,5]
X = data[,c(2,3,4)]
link = "cloglog"
#link = "logit"
corstructure = "AR1"
#corstructure = "Ind"






#=================================================================================================
#Estimating by beta rectangular (Ribeiro, 2020)
#=================================================================================================

betarreg = function(ind, y, X, data, link, corstructure){


#Previous informations -----------------------------------------------------
N = dim(data)[1]
X = data.matrix(X)
y = data.matrix(y)
ind = as.vector(ind)

#Initial values for model-----------------------------------------------------
model.reg = betareg(y ~ X, link = link)
betas = matrix(c(model.reg$coef$mean), nrow = dim(X)[2]+1, ncol=1)
uns = rep(1, N)
X = cbind(uns, X)
alpha = 2*sum(y^2)/(sum(y^2)+N)
phi = model.reg$coef$precision


#Estimation process -----------------------------------------------------
dif = 1

while (dif > 0.01){

eta = X%*%betas

#Link function
if(link == "cloglog"){
  gamaa = as.vector(1-exp(-exp(eta))) 	   
  g = exp(eta)*exp(-exp(eta)) 
}

if(link == "logit"){
  gamaa = as.vector(exp(eta)/(1+exp(eta)))  
  g = exp(eta)*exp(-exp(eta)) 
}



#Variance of Y
nucleo = sqrt(1-4*alpha*gamaa*(1-gamaa))
theta = 1 - nucleo
mu = (gamaa - 1/2 + 1/2*nucleo)/nucleo
var.y = (1/(4*nucleo*(1+phi)))*(4*gamaa*(1-gamaa)+nucleo^2-1)*(2-nucleo+phi*(1-nucleo)) + (1-nucleo)*(1+3*nucleo)/12 #By Vinicius
a = var.y


#Matrices A and G
A = diag(as.vector(a))
G = diag(as.vector(g))


#SI for QIC measure
SI = t(X)%*%G%*%A%*%G%*%X
    

#Vector with mean zero
b = y - gamaa

freq1 = unique(ind)
freq2 = data.frame(table(ind))[,2] #quantas vezes cada individuo aparece
freq3 = cumsum(freq2) #numero de vezes acumuladas
freq4 = c(0,freq3) #adapatação para incluir o zero

numInd = length(freq1) #completo é 31, pois é o numero do ultimo individuo
numTime = max(freq2) #numero maximo de vezes em que um individuo se repete

B = list()
for (i in 2:numInd+1) {
  interval = seq(freq4[i-1]+1, freq4[i])
  auxB = c(b[interval])
  B[[i]] = c(auxB, rep(0, numTime - length(auxB)))
}

aux.b = do.call("rbind", B)


    
# Correlation structure:
if(corstructure == "AR1"){
  m = numTime
  num = sum(apply(aux.b[,-m]*aux.b[,-1],1,sum))
  den = sqrt(sum(apply(aux.b[,-m]^2,1,sum))*sum(apply(aux.b[,-1]^2,1,sum)))
  rho  = num/den
  
  library(dae)
  R.rho = list()
  for (i in 1:numInd) {
    R.rho[[i]] = AR1(freq2[i],rho)
  }
  R.rho = mat.dirsum(R.rho)
}


if(corstructure == "Ind"){
  R.rho = diag(rep(1,N))
  rho = 0
}







#Others matrices
Lamb = G
Omega = A^(1/2)%*%R.rho%*%A^(1/2)
W = Lamb%*%solve(Omega)%*%Lamb
z = eta + solve(Lamb)%*%b
matriz1 = t(X)%*%W%*%X
matriz2 = t(X)%*%W%*%z
    
#Calculating betas
beta.atual = solve(matriz1)%*%matriz2
    
#Selection criterion
dif = norm(beta.atual-betas)
betas = beta.atual
    
    
  
#EM algorithm
dif2 = 1
  
while (dif2 > 0.01){
      
par.est.anterior = matrix(c(phi, alpha))
      
#Dado o chute, definimos u:
theta.u = 1 - sqrt(1 - 4*par.est.anterior[2,1]*(unlist(gamaa)*(1-unlist(gamaa))))
mu.u = (unlist(gamaa) - 1/2 + 1/2*sqrt(1 - 4*par.est.anterior[2,1]*(unlist(gamaa)*(1-unlist(gamaa)))))/
sqrt(1 - 4*par.est.anterior[2,1]*(unlist(gamaa)*(1-unlist(gamaa))))
      
u = theta.u/(theta.u + (1-theta.u)*dbetanova(y,mu.u,par.est.anterior[1,1]))
    
#Dado u, estimamos alpha e phi:
log.vero = function(par){
    phi = par[1]
    alpha = par[2]
        
    theta.aux = 1 - sqrt(1 - 4*alpha*(unlist(gamaa)*(1-unlist(gamaa))))
    mu.aux = (unlist(gamaa) - 1/2 + 1/2*sqrt(1 - 4*alpha*(unlist(gamaa)*(1-unlist(gamaa)))))/
      sqrt(1 - 4*alpha*(unlist(gamaa)*(1-unlist(gamaa))))
        
    log.vero.betar = sum(u*log(theta.aux) + (1-u)*log(1-theta.aux) + (1-u)*log(dbetanova(y,mu.aux,phi)))
    return(-log.vero.betar)
    }
    
    #est1 = optim(par.est.anterior, log.vero, control=list(ndeps=rep(1e-8,4)))
    #estimativas = est1$par
    est2 = optim(par.est.anterior, log.vero, gr=NULL, method='BFGS', hessian=TRUE)
    estimativas = est2$par
    #est2 = optim(par.est.anterior, log.vero, gr=NULL, method='BFGS', hessian=TRUE, control=list(ndeps=rep(1e-8,4)))
    #estimativas = est2$par
    #method = c("nlm", "BFGS", "Rcgmin", "nlminb")
    
    par.est.atual = matrix(c(estimativas[1], estimativas[2]))
      
dif2 = norm(par.est.atual - par.est.anterior)
par.est.anterior = par.est.atual
phi = par.est.anterior[1,1]
alpha = par.est.anterior[2,1]

}


#Output
mylist = list("beta" = betas, "rho" = rho, "phi" = phi, "alpha" = alpha)
#mylist = list("beta" = betas, "rho" = rho, "phi" = phi, "alpha" = alpha, "W" = W, "X" = X, "Lamb" = Lamb, "Omega" = Omega, "b" = b)
return(mylist)

}

}

betas
rho
phi
alpha


#For phi and alpha
est2
SE = sqrt(diag(solve(est2$hessian)))
cbind(estimate=c(phi=est2$par[1],alpha=est2$par[2]),SE)











#=================================================================================================
# Outputs for function
#=================================================================================================
data = read.csv(choose.files(), header=T, sep=";", dec=",")
data = data.frame(cbind(data[5], log(data[2]), (log(data[2]))^2, data[3], data[4]))
head(data)

ind = data[,1]
y = data[,5]
X = data[,c(2,3,4)]
link = "cloglog"
#link = "logit"
corstructure = "AR1"
#corstructure = "Ind"

model = betarreg(ind, y, X, data, link, corstructure)
model$beta
model$phi
model$rho
model$alpha














#=================================================================================================
# Jacknife
#=================================================================================================
data = read.csv(choose.files(), header=T, sep=";", dec=",")
data = data.frame(cbind(data[5], log(data[2]), (log(data[2]))^2, data[3], data[4]))
head(data)


phi.est = c()
rho.est = c()
alpha.est = c()

for (i in 1:31) {
  
  data1 = subset(data, ind != i)
  ind = data1[,1]
  y = data1[,5]
  X = data1[,c(2,3,4)]
  link = "cloglog"
  corstructure = "AR1"
  
  phi.est[i] = betarreg(ind, y, X, data1, link, corstructure)$phi
  rho.est[i] = betarreg(ind, y, X, data1, link, corstructure)$rho
  alpha.est[i] = betarreg(ind, y, X, data1, link, corstructure)$alpha
}


phi1 = phi.est
phi1 = phi1[-c(11,12,15,24,28)]
phiJacknife = mean(phi.est)
phiEpJacknife = ((30/31)*sum((phi1 - phiJacknife)^2))^(1/2);phiEpJacknife



rho1 = rho.est
rho2 = mean(rho.est)
rho.ep = sqrt((30/31)*sum((rho1 - rho2)^2))

alpha1 = alpha.est
alpha2 = mean(alpha.est)
alpha.ep = sqrt((30/31)*sum((alpha1 - alpha2)^2))

phi.ep
rho.ep
alpha.ep










#=================================================================================================
# Bootstrap não paramétrico
#=================================================================================================
data = read.csv(choose.files(), header=T, sep=";", dec=",")
data = data.frame(cbind(data[5], log(data[2]), (log(data[2]))^2, data[3], data[4]))
head(data)

phi.est = c()
rho.est = c()
alpha.est = c()
B = 1000

for (i in 1:B) {
  
  indBoot = sample(1:31, size = 31, replace = TRUE)

  dataBootaux = list()
  for (j in 1:31) {
    dataBootaux[[j]] = cbind(j, subset(data, ind == indBoot[j]))
  }
  dataBoot = do.call("rbind", dataBootaux)

  ind = dataBoot[,1]
  y = dataBoot[,6]
  X = dataBoot[,c(3,4,5)]
  link = "cloglog"
  corstructure = "AR1"
  
  phi.est[i] = betarreg(ind, y, X, dataBoot, link, corstructure)$phi
  rho.est[i] = betarreg(ind, y, X, dataBoot, link, corstructure)$rho
  alpha.est[i] = betarreg(ind, y, X, dataBoot, link, corstructure)$alpha
}


phi1 = phi.est
phiEstBoot = mean(phi.est)
phiEp = ((1/(B-1))*sum((phi1 - phiEstBoot)^2))^(1/2)
phiEp

phiBias = phiEstBoot - 3.2646
phiBias
phiCorrected = 3.2646 - phiBias
phiCorrected


phiEp = ((1/(B-1))*sum((phi1 - phiCorrected)^2))^(1/2)
phiEp






rho1 = rho.est
rho2 = mean(rho.est)
rhoBoot = ((1/(B-1))*sum((rho1 - rho2)^2))^(1/2)

alpha1 = alpha.est
alpha2 = mean(alpha.est)
alphaBoot = ((1/(B-1))*sum((alpha1 - alpha2)^2))^(1/2)

phiBoot
rhoBoot
alphaBoot








#=================================================================================================
# Bootstrap paramétrico
#=================================================================================================
data = read.csv(choose.files(), header=T, sep=";", dec=",")
data = data.frame(cbind(data[5], log(data[2]), (log(data[2]))^2, data[3], data[4]))
head(data)















#=================================================================================================
#Graph of fitted values
#=================================================================================================

beta.est = c(model$beta[1], model$beta[2], model$beta[3], model$beta[4])
eta.fit = beta.est[1] + beta.est[2]*data[,2] + beta.est[3]*data[,3] + beta.est[4]*data[,4]
fit = 1 - exp(-exp(eta.fit))  #cloglog


dadosaux = data
aux.con = c(rep("15%",34), rep("20%",74), rep("25%",73))
with.fit = data.frame(cbind(tempo, concentração, fit))

ggplot(dadosaux, aes(x = dadosaux[,2], y = dadosaux[,4], shape = factor(aux.con))) +
  geom_point() +
  scale_x_continuous(name="Time (days)", breaks=c(0, 20, 40, 60 , 80)) +
  scale_y_continuous(name="Gas concentration (%)", breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))+
  scale_shape_manual("Gas", values = c(19, 3, 1))+
  theme(axis.text.y = element_text(angle = 90))+
  geom_line(aes(x = with.fit[,1], y = with.fit[,3], linetype = factor(with.fit[,2])))+
  scale_linetype_manual("Gas", values=c("solid", "dashed", "dotted"), labels = c("15%", "20%", "25%"))+
  theme(panel.background = element_rect(fill = 'white', colour = 'black'))+
  theme(legend.position = c(0.9, 0.8))
  #scale_fill_discrete(name = "Gas", labels = c("15%", "20%", "25%"))
  
  
  
  
  

  
  
  
  
  
#=================================================================================================
#Fit values
#=================================================================================================
beta.est = c(betas[1], betas[2], betas[3], betas[4])
eta.fit = beta.est[1] + beta.est[2]*data[,2] + beta.est[3]*data[,3] + beta.est[4]*data[,4]
fit = 1 - exp(-exp(eta.fit))  #cloglog
  

  









#=================================================================================================
#Standard error and p-value
#=================================================================================================
X = model$X
W = model$W
Lamb = model$Lamb
b = model$b
Omega = model$Omega

S = -t(X)%*%W%*%X
V = t(X)%*%Lamb%*%solve(Omega)%*%b%*%t(b)%*%solve(Omega)%*%Lamb%*%X

Jrobusto = round(solve(t(S))%*%V%*%solve(S), 4);Jrobusto
Jnaive = round(-solve(S), 4);Jnaive

dpJrobusto = round(sqrt(diag(Jrobusto)), 4);dpJrobusto
dpJnaive = round(sqrt(diag(Jnaive)), 4);dpJnaive




#p-valor:
t0 = betas[1]/(dpJrobusto[1])
p0 = 2*(1-pnorm(abs(t0),0,1))

t1 = betas[2]/(dpJrobusto[2])
p1 = 2*(1-pnorm(abs(t1),0,1))

t2 = betas[3]/(dpJrobusto[3])
p2 = 2*(1-pnorm(abs(t2),0,1))

t3 = betas[4]/(dpJrobusto[4])
p3 = 2*(1-pnorm(abs(t3),0,1))

valor.p = c(p0, p1, p2, p3)

tabela = cbind(betas, dpJrobusto, valor.p);tabela



#For phi and alpha
est
SE = sqrt(diag(solve(est$hessian)))
cbind(estimate=c(phi=est$par[1],alpha=est$par[2]),SE)







#=================================================================================================
#Measure QIC
#=================================================================================================
betaaa = (gamma(phi)/gamma(mu*phi)*gamma((1-mu)*phi))*(y^(mu*phi))*((1-y)^((1-mu)*phi - 1))
LogVer = sum(log(theta + (1-theta)*betaaa))
#SI = t(X)%*%G%*%A%*%G%*%X #already calculated
VarBeta = solve(S)%*%V%*%solve(S)
QIC = -2*LogVer+2*sum(diag(SI%*%VarBeta))
QIC


















#=================================================================================================
# Diagnostics measures - Resíduos, alavanca e distância de Cook
#=================================================================================================
dados = read.csv(choose.files(),header=T, sep=";", dec=",")
attach(dados)
head(dados)
indice = seq(1,181,1)
head(dados)

nome = paste("(", paciente, "," , tempo, ")")




#Resíduo ordinário:
ro = chol(W)%*%(z-eta)
plot(indice, ro, xlab="Índice", ylab="Resíduo",pch=16)
#identify(indice,ro, labels = nome, plot=TRUE,n = 1, cex=0.7)



#Matriz hat:
problema = solve(t(X)%*%W%*%X)
h = diag(chol(W)%*%X%*%problema%*%t(X)%*%chol(W))
plot(indice, h, xlab="Índice", ylab="Resíduo",pch=16)



#Resíduo padronizado:
rp = ro/(sqrt(1-h))
plot(indice,rp, xlab="Paciente", ylab="Resíduo padronizado",pch=16)
#identify(indice,rp, labels = nome, plot=TRUE,n = 4, cex=0.7)



#Distância de Cook:
dc = (1/4)*rp^2*(h/(1-h))
plot(indice,dc, xlab="Paciente", ylab="Distância de Cook",pch=16)
#identify(indice, dc, labels = nome, plot=TRUE,n = 2, cex=0.7)


















#=================================================================================================
#Local Influence 
#=================================================================================================

#Information for local Influence
S = -t(X)%*%W%*%X #Sensibility Matrix
b = y - gamaa #vector b = y - gamma
N = 181 #size of sample


#Local Influence ---> Case-Weights
Delta1 = t(X)%*%W%*%solve(Lamb)%*%diag(as.vector(b),N)
B1 = -t(Delta1)%*%solve(S)%*%Delta1
lmax1 = abs(eigen(B1, symmetric=T)$vectors[,1])
plot(fit, lmax1, xlab="Predicted value", ylab="Case weighting pertubation",pch=16)
identify(fit,lmax1, labels = nome, plot=TRUE,n = 4, cex=0.7)




#Local Influence ---> Response Pertubation
BB = diag(sqrt(var.y))
Delta2 = t(X)%*%W%*%solve(Lamb)%*%BB
B2 = -t(Delta2)%*%solve(S)%*%Delta2
lmax2 = abs(eigen(B2, symmetric=T)$vectors[,1])
plot(fit,lmax2, xlab="Predicted value", ylab="Response perturbation",pch=16)
identify(fit, lmax2, labels = nome, plot=TRUE,n = 5, cex=0.7)




#Local Influence ---> Covariate Pertubation Xk
k = 2 #position of parameter pertubation
p = 4 #number col of X
sxk =  sd(X[,k])
dX = matrix(0,p,N)
dX[k,] = sxk
db = -G*betas[k]*sxk
dG = diag(as.vector((1-exp(eta))*exp(eta)*exp(-exp(eta))),N)	 #derivate of G when link mu is cloglog
dLambda = dG*betas[k]*sxk


aux1 = (4*alpha*gamaa^2 - 4*alpha*gamaa + 1)^(3/2)
aux2 = (-1)*((((alpha-1)*(2*gamaa-1)*(-2*alpha*gamaa^2*phi + phi*aux1 - 4*alpha*gamaa^2 + aux1 + 2*alpha*gamaa^2*phi+ 4*alpha*gamaa - phi - 2)))/((phi+1)*aux1)
             + alpha - 2*alpha*gamaa + (2*alpha*gamaa)/(3*nucleo) - (alpha)/(3*nucleo))
EE = diag(aux2*diag(betas[k]*sxk*dG))#derivate of A
AA = solve(sqrt(A))%*%EE
dOmega = (1/2)*betas[k]*sxk*(sqrt(A)%*%R.rho%*%AA + AA%*%R.rho%*%sqrt(A)) 
dinvOmega = -solve(Omega)%*%dOmega%*%solve(Omega)
Delta3 = t(X)%*%Lamb%*%(solve(Omega)%*%db+dinvOmega%*%diag(as.vector(b),N))+(t(X)%*%dLambda+dX%*%Lamb)%*%solve(Omega)%*%diag(as.vector(b),N)
B3 = -t(Delta3)%*%solve(S)%*%Delta3
lmax3 = abs(eigen(B3, symmetric=T)$vectors[,1])
plot(fit,lmax3, xlab="Predicted value", ylab="Covariate perturbation",pch=16)
identify(fit, lmax3, labels = nome, plot=TRUE,n = 3, cex=0.7)


#





























#=================================================================================================
#Simulated envelope
#=================================================================================================
#5x5
#alpha = 0.01
N=181
R = R.rho

#Function that generates beta rectangular dependent random variables
va_br_dep = function(size_sample){
  z = mvrnorm(size_sample, mu = rep(0, size_sample), Sigma = R, empirical = T)
  u = pnorm(z)
  #Simulating the data from the model parameters:
  resp = list()
  for (i in 1:size_sample) {
    resp[i] = qbetar(u[i,1], fit[i], phi, alpha)
  }
  y_sim = unlist(resp) #This vector contains the 181 simulated values
  return(y_sim)
}
va_br_dep(N)





#Construindo o envelope simulado com resíduo padronizado
new.res.geral = list()
res.ordenados = list()
for (i in 1:25) {
  #Gerando novas amostras correlacionadas
  yij = va_br_dep(N)
  #Calculando o resíduo:
  new.res.geral = (yij - fit)/(sqrt(1-h))                 #padronizado
  #new.res.geral = (yij - fit)/sqrt(phi*mu*(1-mu)*(1-h))         #pearson
  res.ordenados[[i]] = sort(abs(new.res.geral))
}

erro = (yij - fit)/(sqrt(1-h))

#Estatísticas de ordem
valor = matrix(nrow = 25, ncol=N)
for (k in 1:N) {
  for (j in 1:25) {
    valor[j,k] = res.ordenados[[j]][k]
  }
}

minimos = 0
medianas = 0
maximos = 0
for (k in 1:N){
  minimos[k] = min(valor[,k])
  medianas[k] = median(valor[,k])
  maximos[k] = max(valor[,k])
}

min_med_max = matrix(c(minimos, medianas, maximos), nrow = N, ncol = 3)

Z=0
for(i in 1:N)   {Z[i]<-qnorm((i+N-1/8)/(2*N+1/2))}

linhas = cbind(Z,min_med_max,sort(abs(erro)))
faixa = range(linhas[,5],linhas[,2],linhas[,4])

#Graph
plot(linhas[,1],linhas[,5],xlab="Half-normal quantiles",
     ylab="Standardized residual", ylim=faixa, pch=20, col =  1)
par(new=TRUE)
lines(linhas[,1],linhas[,2])
lines(linhas[,1],linhas[,3],lty=2)
lines(linhas[,1],linhas[,4])

#




















































#Extra
################################################################################
#Construindo o envelope simulado com resíduo de pearson padronizado
new.res.geral = list()
res.ordenados = list()
for (i in 1:25) {
  #Gerando novas amostras correlacionadas
  yij = va_br_dep(N)
  #Calculando o resíduo padronizado:
  new.res.geral[[i]] = (yij - fit)/sqrt(phi*mu*(1-mu)*(1-h))         #pearson
  #new.res.geral[[i]] = (yij - ajustado.vini)/(sqrt(1-h.geral))                 #padronizado
  res.ordenados[[i]] = sort(abs(new.res.geral[[i]]))
}

erro = (yij - fit)/sqrt(phi*mu*(1-mu)*(1-h))

#Estatísticas de ordem
valor = matrix(nrow = 25, ncol=N)
for (k in 1:N) {
  for (j in 1:25) {
    valor[j,k] = res.ordenados[[j]][k]
  }
}

minimos = 0
medianas = 0
maximos = 0
for (k in 1:N){
  minimos[k] = min(valor[,k])
  medianas[k] = median(valor[,k])
  maximos[k] = max(valor[,k])
}

min_med_max = matrix(c(minimos, medianas, maximos), nrow = N, ncol = 3)


Z=0
for(i in 1:N)   {Z[i]<-qnorm((i+N-1/8)/(2*N+1/2))}


linhas = cbind(Z,min_med_max,sort(abs(erro)))
faixa = range(linhas[,5],linhas[,2],linhas[,4])

#Construindo o gráfico
plot(linhas[,1],linhas[,5],xlab="Valor esperado da estatística de ordem meio-normal",
     ylab="Valor absoluto ordenado do resíduo padronizado", ylim=faixa, pch=16)
par(new=TRUE)
lines(linhas[,1],linhas[,2])
lines(linhas[,1],linhas[,3],lty=2)
lines(linhas[,1],linhas[,4])









#----------------------------------------------------------------------
### Density
#----------------------------------------------------------------------


#Antes
a.b0 = 0.6629
a.b1 = 0.1879
a.b2 = -0.1843
a.b3 = 0.2503
a.phi = 3.2646
a.alpha = 0.1572

b.b0 = 0.6790
b.b1 = 0.1679
b.b2 = -0.1804
b.b3 = 0.2515
b.phi = 5.7884
b.alpha = 0.5575


x = seq(0, 1, 0.01)
par(mfrow= c(1,2))
plot(dbetar(x, a.b0, a.phi, a.alpha))
plot(dbetar(x, b.b0, b.phi, b.alpha))

require(e1071)
kurtosis(dbetar(x, a.b0, a.phi, a.alpha))
kurtosis(dbetar(x, b.b0, b.phi, b.alpha))





