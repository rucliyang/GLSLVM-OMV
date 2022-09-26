
ds_or=function(x,jj,zh,df){
n=nrow(x)
p=ncol(x)
m=length(table(x[,1]))
#m=max(x)-min(x)+1 
########
y=list()
for(j in 1:p){
  xx=matrix(0,n,m)
  for(i in 1:n){
    xx[i,x[,j][i]]=1
             }
  y[[j]]=xx
}
########
z=array(0,dim=c(p,m-1,n))
for(j in 1:p){
  for(k in 1:(m-1)){
    for(i in 1:n){
    z[j,k,i]=exp(jj[j,k]+zh[j,1]*df[i,1])
                 }
                }
             }
z[z>1e30]=1e30  ##

#################
ff=matrix(0,n,p)
f1=matrix(0,n,1)
c1=array(0,dim=c(n,p,m))

#################
fff=matrix(0,n,p)
f2=matrix(0,n,1)#二阶导
c2=array(0,dim=c(n,p,m))

###############
for(i in 1:n){
  for(j in 1:p){ 
    for(k in 1:m){
      if(k==1) {
       c1[i,j,k]=1/(1+z[j,k,i])
              }
      else if(k==m){ 
       c1[i,j,k]=-z[j,k-1,i]/(1+z[j,k-1,i])
                  }
      else {
       c1[i,j,k]=(1-z[j,k,i]*z[j,k-1,i])/((1+z[j,k,i])*(1+z[j,k-1,i]))
            }
                 }
    ff[i,j]=zh[j,1]*(y[[j]][i,]%*%c1[i,j,])
                }
  f1[i,1]=sum(ff[i,])
              }

##########################
for(i in 1:n){
  for(j in 1:p){
    for(k in 1:m){
      if(k==1){ 
       c2[i,j,k]=z[j,k,i]/((1+z[j,k,i])^2)
              }
      else if(k==m){ 
       c2[i,j,k]=z[j,k-1,i]/((1+z[j,k-1,i])^2)
                   }
      else {
       c2[i,j,k]=(4*z[j,k,i]*z[j,k-1,i]+(z[j,k,i]+z[j,k-1,i])*(1+z[j,k,i]*z[j,k-1,i]))/(((1+z[j,k,i])*(1+z[j,k-1,i]))^2)
           }
                 }
    fff[i,j]=(zh[j,1]^2)*y[[j]][i,]%*%c2[i,j,]
               }
   f2[i,1]=-sum(fff[i,])
             }
out=list(f1,f2)
return(out)
}

########################
zh_ordinal=function(x,jj,zh,df){
n=nrow(x)
p=ncol(x)
m=length(table(x[,1]))
#m=max(x)-min(x)+1
##########哑变量
y=list()
for(j in 1:p){
  xx=matrix(0,n,m)
  for(i in 1:n){
    xx[i,x[,j][i]]=1
               }
  y[[j]]=xx
             }
##################
z=array(0,dim=c(p,m-1,n))
for(j in 1:p){
  for(k in 1:(m-1)){
    for(i in 1:n){
      z[j,k,i]=exp(jj[j,k]+zh[j,1]*df[i,1])
                 }
                }
             }
z[z>1e30]=1e30  
#################
gg=matrix(0,n,p)
g1=matrix(0,p,1)
d1=array(0,dim=c(n,p,m))
#################
ggg=matrix(0,n,p)
g2=matrix(0,p,1)
d2=array(0,dim=c(n,p,m))

############################
for(j in 1:p){
  for(i in 1:n){    
    for(k in 1:m){
      if(k==1) {
       d1[i,j,k]=1/(1+z[j,k,i])
                }
      else if(k==m) {
       d1[i,j,k]=-z[j,k-1,i]/(1+z[j,k-1,i])
                    }
      else {
       d1[i,j,k]=(1-z[j,k,i]*z[j,k-1,i])/((1+z[j,k,i])*(1+z[j,k-1,i]))
            }
                 }
    gg[i,j]=df[i,1]*(y[[j]][i,]%*%d1[i,j,])
                }
   g1[j,1]=sum(gg[,j])
              }
#######################################
for(j in 1:p){
  for(i in 1:n){
    for(k in 1:m){
      if(k==1) {
       d2[i,j,k]=z[j,k,i]/((1+z[j,k,i])^2)
               }
      else if(k==m) {
       d2[i,j,k]=z[j,k-1,i]/((1+z[j,k-1,i])^2)
                    }
      else {
       d2[i,j,k]=(4*z[j,k,i]*z[j,k-1,i]+(z[j,k,i]+z[j,k-1,i])*(1+z[j,k,i]*z[j,k-1,i]))/(((1+z[j,k,i])*(1+z[j,k-1,i]))^2)
            }
                 }
    ggg[i,j]=(df[i,1]^2)*y[[j]][i,]%*%d2[i,j,]
               }
   g2[j,1]=-sum(ggg[,j])
             }
###################################
for(j in 1:p){   
  zh[j,1]=zh[j,1]-g1[j,]/g2[j,1]
             }
return(zh)
}

