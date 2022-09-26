#ordinal;MASS;"ordinal function.R";"upls.R";'nipals.R'

stage1 = function(x, block){
  #x should be matrix
  #block is list

  n=nrow(x)
  categ=length(table(x[,1]))
  jj = matrix(1,nrow=ncol(x),ncol=categ-1)

  b = matrix(1,nrow=ncol(x),ncol=1)
  
  Q = length(block)
  esp = rep(1,Q+1)
  ksi = matrix(1,nrow=length(x[,1]),ncol=Q)
  eta = c()
  w = rep(1, Q)

  w=rep(1,Q)
  ee=0.1
  s=1
  pin=apply(x,2,table)/n
  ppin=apply(pin,2,cumsum)

  while(ee>0.000001 & s<300){  
    
    for(q in 1:Q){
      #q=1
      jj2 = jj[block[[q]],]

      kk1 = ksi[,q]
      G = ds_or(x[,block[[q]]],jj2,as.matrix(b[block[[q]],]),as.matrix(ksi[,q]))
      ggg1=G[[1]]
      ggg2=G[[2]]

      ksi[,q] = ksi[,q] - ggg1/ggg2
      ksi[,q] = scale(ksi[,q])*sqrt((n-1)/n)

      esp[q] = max(abs((kk1-ksi[,q])/kk1))

      B = zh_ordinal(x[,block[[q]]],jj2,as.matrix(b[block[[q]],]),as.matrix(ksi[,q]))
      B = -abs(B)
      b[block[[q]],] = B/c(sqrt(t(B)%*%B))

      jieju=matrix(nrow=ncol(x[,block[[q]]]),ncol=categ-1)
        for(j in 1:(categ-1)){
          jieju[,j] = quantile(ksi[,q],ppin[,block[[q]]][j,])
                             }
      jj[block[[q]],] = jieju
                   }
      
      eta = nipals(ksi)$scores[,1]     
      w0 = w
      w = cor(ksi,eta)
      esp[Q+1] = max(abs((w0-w)/w0))

      ee = max(esp)
      s = s+1
      
                           }
      #end while
      return(list(ksi=ksi,eta=eta))
               }
      #end stage1 function

scale0 = function(x){
  x = x/x[1]
  return(x)         
                    }

stage2 = function(x, block, ksi, eta){
  #x = dat1;block=list(1:9,10:18,19:27);ksi=ksi1;eta=eta1
  lambda = lambda0 = c();delta = delta0 = c()
  for(q in 1:length(block)){
    #q=3
    for(i in block[[q]]){
    lambda[i] = c(clm(as.factor(x[,i]) ~ ksi[,q])$beta)
                        }
  
  lambda0[block[[q]]] = scale0(lambda[block[[q]]])
  #delta[q]=lm(ksi[,q] ~ eta-1)$coefficients
  delta[q] = lm(nipals(x[,block[[q]]])$scores[,1] ~ eta-1)$coefficients
                           }
  delta0 = scale0(delta)
  return(list(lambda=lambda, lambda0=lambda0, delta=delta, delta0=delta0))
                                     }  

 
 
  


