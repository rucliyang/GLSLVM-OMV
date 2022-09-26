library(ordinal)
library(MASS)
source("ordinal function.R")
source('estima_prop.R')
source('nipals.R')

dat = read.csv('sd3.csv',sep='\t')
dat = dat[,1:27]

dat[dat==0] = NA
datt = na.omit(dat)

set.seed(123)
n1 = sample(1:nrow(datt),600,replace=F)
dat1 = datt[n1,]

#dat1 = datt

lev = c()
for(i in 1:27){
 lev = rbind(lev, table(dat1[,i]))
               };lev

fit_score = stage1(dat1, list(1:9,10:18,19:27))
ksi1 = fit_score$ksi
eta1 = fit_score$eta

fit_other = stage2(dat1, list(1:9,10:18,19:27), ksi1, eta1)

fit_other$lambda
fit_other$lambda0
fit_other$delta
fit_other$delta0

#bootstrap
set.seed(123)
B = 500
lambda_B = matrix(nrow=27,ncol=B)
delta_B = matrix(nrow=3,ncol=B)
lambda_B0 = matrix(nrow=27,ncol=B)
delta_B0 = matrix(nrow=3,ncol=B)
ksi_B = array(0, dim = c(600,3,B))
eta_B = matrix(nrow=600, ncol=B)

for(b in 1:B){
  print(b)
  n_b = sample(1:600,replace=T)
  dat_b = dat1[n_b,]
  score_b = stage1(dat_b, list(1:9,10:18,19:27))
  ksi_B[,,b] = score_b$ksi
  eta_B[,b] = score_b$eta
  other_b = stage2(dat_b, list(1:9,10:18,19:27), score_b$ksi, score_b$eta)
  lambda_B[,b] = other_b$lambda
  delta_B[,b] = other_b$delta
  lambda_B0[,b] = other_b$lambda0
  delta_B0[,b] = other_b$delta0
              }

out_index = function(x)
{
  xout = boxplot.stats(x)$out
  idx = which(x %in% xout)
  return(idx) 
                       }

indok = c(1:200)[-Reduce(union, apply(lambda_B, 1, out_index))]


seB = function(x){
  n = nrow(x);p=ncol(x)
  bmean = apply(x,1,mean)
  bse = sqrt(apply((x - bmean)^2,1,sum)/p)
  return(list(bmean=bmean,bse=bse))
                 }

lab_B_out = seB(lambda_B[,indok])
del_B_out = seB(delta_B[,indok])

lab_B0_out = seB(lambda_B0[,indok])
del_B0_out = seB(delta_B0[,indok])

