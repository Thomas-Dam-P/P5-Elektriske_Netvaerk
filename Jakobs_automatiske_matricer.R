# order of data is voltages first, then currents
# voltages ordered according to vertex number
# currents ordered according to from-vertex first, then to-vertex

### Data
res = read.table("R_small_grid.txt")
curr = c(-2,-160,-3320,-400,1320,-80,-400,-4160,8960,-280)
nedg = sum(res!=0); nver = dim(res)[1]
n = nver+nedg

### matrix containing numbering of edge from vertex i to vertex j
ver2edge = res
k = 1
for (i in 1:nver){
  for (j in 1:nver){
    if (res[i,j]!=0){
      ver2edge[i,j] = k; k = k+1
    }
  }
}

### C matrix
C = matrix(nrow=0,ncol=n)
for (i in 1:nver){   # Ohm's law
  for (j in 1:nver){
    if (res[i,j]!=0){
      temp = rep(0,n); temp[i] = 1; temp[j] = -1; temp[nver+ver2edge[i,j]] = res[i,j]
      C = rbind(C,temp,deparse.level = 0)
    }
  }
}
for (i in 1:nver){   # Kirchhoff's law
  if (sum(res[,i])>0 & sum(res[i,])>0){
    ingoing = ver2edge[,i]; ingoing = ingoing[ingoing>0]
    outgoing = ver2edge[i,]; outgoing = outgoing[outgoing>0]
    temp = rep(0,n); temp[nver+ingoing] = -1; temp[nver+outgoing] = rep(1,length(outgoing))
    C = rbind(C,temp,deparse.level = 0)
  }
}

### D matrix 
D = matrix(nrow=0,ncol=n)
measp = c(3,4,7,9,11,13,15,17,18,19)  # ID of measured powers, prosumers
#for (i in 1:length(measv)){
#  temp = rep(0, n); temp[measv[i]] = 1
#  D = rbind(D,temp,deparse.level = 0) 
#}
temp = c(1,rep(0,n-1))   # voltage at substation
D = rbind(D,temp,deparse.level = 0)
for (i in 1:length(measp)){  # currents at prosumers
  ingoing = ver2edge[,measp[i]]; ingoing = ingoing[ingoing>0]
  temp = rep(0,n); temp[nver+ingoing] = 1   # add nver to get ID of measured currents
  D = rbind(D,temp,deparse.level = 0)
}

### d vector
d = c(400, curr)  # 400V at substation, currents calculated from measured powers

### Sigma
Sigma = diag((0.01*d)^2)    ### std dev of 1% of measurements (should be actual values)

Dt = solve(Sigma)^(1/2) %*% D
dt = solve(Sigma)^(1/2) %*% d
A = rbind(cbind(t(Dt)%*%Dt,t(C)),
          cbind(C,matrix(0,dim(C)[1],dim(C)[1])))
b = rbind(t(Dt)%*%dt,matrix(0,dim(C)[1],1))
xt = solve(A,b); round(xt,2)
xhat = as.matrix(xt[1:dim(C)[2]]); round(xhat,2)
dhat = D%*%xhat;dhat