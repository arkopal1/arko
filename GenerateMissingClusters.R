set.seed(1)
library("clusterGeneration")
mu = matrix(NaN,nrow = 5,ncol = 4)
sigma = array(NaN,c(5,5,4))

# Generating mean vectors and covariance matrices
for (j in 1:4)
{
  mu[,j] = runif(5)
  sigma[,,j] = genPositiveDefMat (5,covMethod = "unifcorrmat")$Sigma
}

# Generating the complete data
library("MASS")
dataFull = array(NaN,c(5,100,4))
for (j in 1:4)
{
  dataFull[,,j] = t(mvrnorm(100,mu[,j],sigma[,,j]))
}

# Missing at Random Structure
d = array(1,c(5,100,4))

# Parameters of the logistic distribution for generating missing observations
pmiss = array(NaN,c(4,2,4))
pmiss[,1,1] = c(0.5,1,1.5,2)
pmiss[,2,1] = c(0.6,1,1.4,1.8)
pmiss[,1,2] = c(0.5,1.2,1.9,2.6)
pmiss[,2,2] = c(0.6,1.2,1.8,2.4)
pmiss[,1,3] = c(0.5,1.5,2.5,3.5)
pmiss[,2,3] = c(0.6,1.5,2.4,3.3)
pmiss[,1,4] = c(0.5,2,3.5,5)
pmiss[,2,4] = c(0.6,2,3.4,4.8)
pmiss = pmiss

# matrix of bernoullis d to create observed matrix with missing variables
library("LaplacesDemon")
for (i in 1:100)
{
  for (j in 1:4)
  {
    prob4 = pmiss[1,1,j]+pmiss[2,1,j]*dataFull[1,i,j]+pmiss[3,1,j]*dataFull[2,i,j]+pmiss[4,1,j]*dataFull[3,i,j]
    d[4,i,j] = rbern(1,(exp(prob4)/(1+exp(prob4))))
    prob5 = pmiss[1,2,j]+pmiss[2,2,j]*dataFull[1,i,j]+pmiss[3,2,j]*dataFull[2,i,j]+pmiss[4,2,j]*dataFull[3,i,j]
    d[5,i,j] = rbern(1,(exp(prob5)/(1+exp(prob5))))
  }
}
dataNMiss = d*dataFull
dataNMiss[dataNMiss==0] = NA
library("mice")
print(md.pattern(t(as.data.frame(dataNMiss[,,1]))))
print(md.pattern(t(as.data.frame(dataNMiss[,,2]))))
print(md.pattern(t(as.data.frame(dataNMiss[,,3]))))
print(md.pattern(t(as.data.frame(dataNMiss[,,4]))))
print(md.pattern(t(as.data.frame(dataNMiss))))

# Starting our analysis of the missing data
p1 = rowMeans(d[,,1]) # Proportion of missingness in each group
p2 = rowMeans(d[,,2])
p3 = rowMeans(d[,,3])
p4 = rowMeans(d[,,4])

# Means of the missing data where each column corresponds to a group
means = cbind(rowMeans(dataNMiss[,,1],na.rm = TRUE),rowMeans(dataNMiss[,,2],na.rm = TRUE),
              rowMeans(dataNMiss[,,3],na.rm = TRUE),rowMeans(dataNMiss[,,4],na.rm = TRUE))
meansAct = cbind(rowMeans(dataFull[,,1],na.rm = TRUE),rowMeans(dataFull[,,2],na.rm = TRUE),
                 rowMeans(dataFull[,,3],na.rm = TRUE),rowMeans(dataFull[,,4],na.rm = TRUE))

# Imputing the observations of the data matrix
distances = matrix(0,4)
for (l in 1:4)
{
  g1 = dataNMiss[,,l]
  d1 = d[,,l]
  g1Full = dataFull[,,l]

  # Reorganizing the missing values
  k = 1
  for (i in 1:100)
  {
    for (j in 1:5)
      if (is.na(g1[j,i])==TRUE)
      {
        k = i
        while(k<100)
        {
          k=k+1
          nope = 0
          for (j in 1:5)
          {
            if (is.na(g1[j,k])==TRUE)
             nope = 1
          }
          if (nope == 0)
          {
            swap = g1[,i]
            g1[,i] = g1[,k]
            g1[,k] = swap
            swap = d1[,i]
            d1[,i] = d1[,k]
            d1[,k] = swap
            swap = g1Full[,i]
            g1Full[,i] = g1Full[,k]
            g1Full[,k] = swap
            break
          }
        }
        break
      }
    if (k==100)
      break
  }

  # Counting the number of complete observations
  k = NA
  for (i in 1:100)
  {
    for (j in 1:5)
    {
      if (is.na(g1[j,i])==TRUE)
      {
        k = i-1
        break
      }
    }
    if (is.na(k)==FALSE)
      break
  }
  
  # Imputation using 3 nearest neighbours
  for (i in (k+1):100)
  {
    dist1 = Inf
    dist2 = Inf
    dist3 = Inf
    min1 = NA
    min2 = NA
    min3 = NA
    for (j in 1:k)
    {
      distance = sqrt(sum((g1[,i]-g1[,j])^2,na.rm=TRUE))
      if (distance < dist1)
      {
        min1 = j
        dist1 = distance
      }
    }
    for (j in 1:k)
      if (j!=min1)
      {
        distance = sqrt(sum((g1[,i]-g1[,j])^2,na.rm=TRUE))
        if (distance < dist2)
        {
          min2 = j
          dist2 = distance
        }
      }
    for (j in 1:k)
      if (j!=min1 && j!= min2)
      {
        distance = sqrt(sum((g1[,i]-g1[,j])^2,na.rm=TRUE))
        if (distance < dist3)
        {
          min3 = j
          dist3 = distance
        }
      }
    for (j in 4:5)
      if (is.na(g1[j,i])==TRUE)
      {
        p = runif(1)
        if (p < (1/dist1)/(1/dist1+1/dist2+1/dist3))
          g1[j,i] = g1[j,min1]
        else if (p < 1-(1/dist3)/(1/dist1+1/dist2+1/dist3))
          g1[j,i] = g1[j,min2]
        else
          g1[j,i] = g1[j,min3]
        #g1[j,i] = (g1[j,min1]/dist1+g1[j,min2]/dist2+g1[j,min3]/dist3)*(1/(1/dist1+1/dist2+1/dist3))
        #g1[j,i] = (g1[j,min1]/dist1+g1[j,min2]/dist2)*(1/(1/dist1+1/dist2))
        #g1[j,i] = g1[j,min1]
      }
  }
  #distances[l] = colDist(g1,g1Full)
  print(cbind(rowMeans(g1),rowMeans(g1Full)))
  print(ks.test(g1[5,c((k+1):100)],g1Full[5,c((k+1):100)],alternative = "two.sided",exact = NULL))
}
#print(distances)

# creating a test matrix
n = 100
class = sample(1:4,n,replace=T)
test = matrix(NaN,5,n)
td = matrix(1,5,n)
for (i in 1:n)
{
  test[,i] = mvrnorm(1,mu[,class[i]],sigma[,,class[i]])
  
  # Introducing Missingness
  prob4 = pmiss[1,1,class[i]]+pmiss[2,1,class[i]]*test[1,i]+pmiss[3,1,class[i]]*test[2,i]+pmiss[4,1,class[i]]*test[3,i]
  td[4,class[i]] = rbern(1,(exp(prob4)/(1+exp(prob4))))
  prob5 = pmiss[1,2,class[i]]+pmiss[2,2,class[i]]*test[1,i]+pmiss[3,2,class[i]]*test[2,i]+pmiss[4,2,class[i]]*test[3,i]
  td[5,class[i]] = rbern(1,(exp(prob5)/(1+exp(prob5))))
}
testNMiss = td*test
testNMiss[testNMiss==0] = NA
test2 = rbind(class,test)

tdist1 = matrix(Inf,1,n)
tclass = matrix(NA,1,n)
for (i in 1:n)
{
  for (j in 1:4)
  {
    for (k in 1:100)
    {
      distance = sqrt(sum((test[,i]-dataNMiss[,k,j])^2,na.rm=TRUE))
      if (distance < tdist1[i])
      {
        tdist1[i] = distance
        tclass[i] = j
      }
    }
  }
}  
print(prop.table(table(class!=tclass)))
