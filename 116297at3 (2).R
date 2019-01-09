foodCalibrate <- read.csv("foodCalibrate.csv")
foodPredict <- read.csv("foodTrain.csv")

#função de predicao
pred.cluster <- function(x, centers) {
  tmp = matrix(NA, nrow = nrow(x), ncol = nrow(centers))
  for(i in seq_len(nrow(x))){
    for (j in seq_len(nrow(centers))){
      tmp[i,j]= sum((x[i,]- centers[j,])^2)
    }
  } 
  max.col(-(tmp))  
}

#calculando BW em dados novos
BW= function(xpred,centers){
  k= nrow(centers)
  n= ncol(xpred)-1  
  #Matriz B
  xbarra= apply(xpred[2:(n+1)],2,mean)
  dB= apply(centers, 1, function(i) i - xbarra)
  B= matrix(0,nrow = n,ncol = n)
  for(l in 1:k){
    nk= length(which(xpred$pred == l))
    dB2= dB[,l] %*% t(dB[,l])
    B= B+ nk*dB2
  }
  
  #Matriz W
  W= matrix(0,ncol = n,nrow = n)
  vazio= 0
  for (i in 1:k){
    idx = which(xpred$pred == i) 
    xk= as.matrix(xpred[idx,2:(n+1)])
    if(nrow(xk)>0){
    for(j in 1: nrow(xk)){
      
      xk[j,]= xk[j,] - centers[i,]
    }
    dW= matrix(0,ncol = n,nrow = n)
    for(j in 1:nrow(xk)){
      dW2= xk[j,] %*% t(xk[j,])
      dW= dW + dW2
    }
    W= W + dW
    }
  else vazio= vazio + 1
  }
  propVar= det(B)/det(B+W)
  saida= t(c(propVar,vazio))
  return(saida)
}

x= foodCalibrate[3:9]
x.scale= as.data.frame(scale(x))
SVD= svd(x.scale)
U = SVD$u
D = diag(SVD$d)
V = SVD$v
scores <- U %*% D
xtrain= as.data.frame(scores[,1:2])

mu= apply(x,2,mean)
sig= apply(x, 2, sd)
set.seed(13)
#idx= sample(nrow(foodPredict),1000)
#xtune= foodPredict[idx,3:9]        
xtune= foodPredict[,3:9]
xtune.std= t(apply(xtune, 1, function(x) (x - mu)/sig))
xtune.pca = as.data.frame(xtune.std %*% V[,1:2])
#2 componetes principais explicam cerca de 50% da variabilidade dos dados de treino.
#Quanto menor o numero de componentes, mais rapido k converge para uma região de corte, 
#segundo a regra do cotovelo, e maior a proporcao da variabilidade explicada nos dados de
#teste. 

k= 9 
erro= matrix(NA,nrow=length(k),ncol = 2)

for(i in 1:length(k)){
  fit= kmeans(xtrain,k[i],iter.max = 100)
  centers= fit$centers
  pred= pred.cluster(xtune.pca,centers)
  xpred = cbind(pred,xtune.pca)
  erro[i,]=BW(xpred,centers)
  #print(i)
}

prop= erro[1]

write(pred, file= "116297label.txt",sep= " ")
write(prop, file= "116297prop.txt")
