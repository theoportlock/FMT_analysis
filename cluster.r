require(cluster)
require(clusterSim)
data=read.table("genus.txt", header=T, row.names=1, dec=".", sep="\t")
data=data[-1,]

JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
KLD <- function(x,y) sum(x * log(x/y))

dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
       KLD <- function(x,y) sum(x *log(x/y))
       JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
       matrixColSize <- length(colnames(inMatrix))
       matrixRowSize <- length(rownames(inMatrix))
       colnames <- colnames(inMatrix)
       resultsMatrix <- matrix(0, matrixColSize, matrixColSize)

inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))

       for(i in 1:matrixColSize) {
       	for(j in 1:matrixColSize) {
       		resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
       		as.vector(inMatrix[,j]))
       	}
       }
       colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
       as.dist(resultsMatrix)->resultsMatrix
       attr(resultsMatrix, "method") <- "dist"
       return(resultsMatrix)
}

data.dist=dist.JSD(data)

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
                        require(cluster)
                        cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
                        return(cluster)
                       }

data.cluster=pam.clustering(data.dist, k=3)

nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids")

nclusters=NULL

for (k in 1:20) { 
	if (k==1) {
		nclusters[k]=NA 
	} else {
		data.cluster_temp=pam.clustering(data.dist, k)
		nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
		centrotypes = "medoids")
	}
}

plot(nclusters, type="h", xlab="k clusters", ylab="CH index")

data.cluster=pam.clustering(data.dist, k=3)

obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])

noise.removal <- function(dataframe, percent=0.01, top=NULL){
	dataframe->Matrix
	bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent
	Matrix_1 <- Matrix[bigones,]
	print(percent)
	return(Matrix_1)
 }

data.denoized=noise.removal(data, percent=0.01)

require(ade4)
obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1)
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F)

s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(4,2,3))
