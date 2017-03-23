t_start<-Sys.time()
.libPaths()

####Read in image data####
library(R.matlab)
b_readdata<-Sys.time()

#dvfx is a 4D array. We need to add a singleton as the 5th dimension
dvfx<-readMat('/home/kangk4/dissertation/data/dvfx.mat')[[1]]
dim1<-length(dvfx[,1,1,1])
dim2<-length(dvfx[1,,1,1])
dim3<-length(dvfx[1,1,,1])
dim4<-length(dvfx[1,1,1,])
dvf<-array(0,c(dim1,dim2,dim3,dim4,1))
dvf[,,,,1]<-dvfx
dvf<-list(dvf)

#dvf<-readMat('/home/kangk4/dissertation/data/dvf.mat')
#dvf1<-readMat('/home/kangk4/dissertation/data/dvf1.mat')
#dvf2<-readMat('/home/kangk4/dissertation/data/dvf2.mat')
#dvf3<-readMat('/home/kangk4/dissertation/data/dvf3.mat')
#If you have the DVF saved as an R object, then use the following will be faster
#load("/home/kangk4/dissertation/data/dvfx.rda")
e_readdata<-Sys.time()
e_readdata-b_readdata
#make a list containing all subjects
s<-list(dvf)
str(s)

#The following algorithm will vectorize a matrix by stacking its columns
mtov<-function(m){
	m    <- data.matrix(m)
	dim1 <- dim(m)[1]
	dim2 <- dim(m)[2]
	v    <- rep(0,dim1*dim2)
	for (i in 1:dim2){
		v[(dim1*(i-1)+1):(dim1*i)] <- m[,i]
	}
	return(v)
}

#The following algorithm will convert a two-level ([[sub]][[phase]]) list of 2D arrays into a matrix by (1) converting 2D arrays into a vector for each level-one+level-two list element and (2) cbind-ing each vector into a longer vector
ltov<-function(l){
	sub    <- length(l)
	phase  <- length(l[[1]])
	dim1   <- length(l[[1]][[1]][,1])
	dim2   <- length(l[[1]][[1]][1,])
	dim    <- dim1*dim2
	sdim   <- dim1*dim2*phase
	v      <- rep(0,dim*sub*phase)
	m      <- matrix(0,nrow=dim1,ncol=dim2)
	for (i in 1:sub){
		for (j in 1:phase){
			m <- l[[i]][[j]]
			v[(sdim*(i-1)+dim*(j-1)+1):(sdim*(i-1)+dim*j)]<-mtov(m)
		}
	}
	return(v)
}

#The following algorithm reaaranges the image structure from a list of 5D arrays into a 2-level ([[sub]][[phase]]) list of 2D arrays [dim1*dim3*3(DVF contains x,y,z cor),dim2] aka for each subject at each phase the algorithm is stacking the transverse slices for each cor in DVF to make a long matrix.
#input needed:
#s: the original DVF image data in stored as a list [[sub]] of 5D arrays [dim1,dim2,dim3,Phase,DVF(1:3)]
restr<-function(s){
	
I<-length(s)
P<-length(s[[1]][[1]][1,1,1,,1])
X<-length(s[[1]][[1]][,1,1,1,1])
Y<-length(s[[1]][[1]][1,,1,1,1])
Z<-length(s[[1]][[1]][1,1,,1,1])
COR<-length(s[[1]][[1]][1,1,1,1,])
rs <- vector('list',I)

	for (i in 1:I)
	{
	for (p in 1:P)
	{
		rs[[i]][[p]] <- matrix(0,nrow=X*Z*COR, ncol=Y)
		for (cor in 1:COR){
			for (z in 1:Z){
				rs[[i]][[p]][(X*Z*(cor-1)+X*(z-1)+1):(X*Z*(cor-1)+X*z),] <- s[[i]][[1]][,,z,p,cor]
			}
		}
	}
	}
	return(rs)
}
rs<-restr(s)
#save(rs,file="/home/kangk4/dissertation/data/pvd_rs.rda")

begin_t_pvd<-Sys.time()
b_t_svd<-Sys.time()
	I<-length(rs)
	P<-length(rs[[1]])
	X<-dim(rs[[1]][[1]])[1]
	Y<-dim(rs[[1]][[1]])[2]

	svdrs <- vector('list',I)
	
	#perform SVD (Y=UDV^T)
	u <- vector('list',I)
	v <- vector('list',I)
	
	sigma <- vector('list',I)
	
	listofu <- list()
	listofv <- list()
	
	for (i in 1:I)
	{
	for (p in 1:P)
	{svdrs[[i]][[p]] <- svd(rs[[i]][[p]],LINPACK=FALSE)}
	}
e_t_svd<-Sys.time()
e_t_svd-b_t_svd

b_t_p1<-Sys.time()	
	#approximate SVD based on variations explianed p1 threshold
	p1<-.9
	for (i in 1:I)
	{
	for (p in 1:P)
	{
		var <- 0
		#for (k in 1:4) #for a fixed number of Vi for each subject set to. It was set to 4 because the max number for Vi is 4 in the prcedure where Vi is not fixed to achieve a certain percentage of variations.
		for (k in 1:length(svdrs[[i]][[p]]$d)) # for difference size of Vi for each subject
		{
			var <- var+svdrs[[i]][[p]]$d[k]^2/sum((svdrs[[i]][[p]]$d)^2)
			
			#if the variantions explianed is above a certain percentage then stop
			if (var > p1) break #for different size of Vi for each subject
		}	
		
		u[[i]][[p]] <- svdrs[[i]][[p]]$u[,1:k]
		v[[i]][[p]] <- svdrs[[i]][[p]]$v[,1:k]

		sigma[[i]][[p]]<- c(svdrs[[i]][[p]]$d[1:k])
	}
	}
	
	#combine each subject's each vist's U and V
	for (i in 1:I)
	{
	for (p in 1:P)
	{
		listofu[[2*i+p-2]] <- u[[i]][[p]]
		listofv[[2*i+p-2]] <- v[[i]][[p]]
	}
	}
	U <- do.call(cbind,listofu)
	V <- do.call(cbind,listofv)
	dim(U)
	dim(V)
e_t_p1<-Sys.time()
e_t_p1-b_t_p1
	
#perform the second level reduction on combined U and V
b_t_ssvd<-Sys.time()
p2<-.9
svdu <- svd(U)
object.size(svdu)
varu <- 0
listofp <- list()
for (x in 1:X)
{
	varu <- varu+svdu$d[x]^2/sum((svdu$d)^2)
	listofp[[x]] <- svdu$u[,x]
	if (varu > p2) break
}

svdv <- svd(V)
object.size(svdv)
varv <- 0
listofd <- list()
for (y in 1:Y)
{
	varv <- varv+svdv$d[y]^2/sum((svdv$d)^2)
	listofd[[y]] <- svdv$u[,y]
	if (varv > p2) break
}

#be careful P here is the matrix in PVD 
#before it is the index for phase
P <- do.call(cbind,listofp)
D <- do.call(cbind,listofd)
dim(P)
dim(D)
e_t_ssvd<-Sys.time()
e_t_ssvd-b_t_ssvd
end_t_pvd<-Sys.time()
end_t_pvd-begin_t_pvd

b_t_recon_v<-Sys.time()
#construct Vi
value <- vector('list',I)
for (i in 1:I)
{
	for (p in 1:length(rs[[1]]))
	{
		if (length(sigma[[i]][[p]])<2){
			#here block matrix multiplication is not necesary b/c matrix products in () have small dim
			value[[i]][[p]] <- (t(P)%*%u[[i]][[p]])%*%(sigma[[i]][[p]]%*%diag(1))%*%(t(v[[i]][[p]])%*%D)
		}
		else
		{
			#here block matrix multiplication is not necesary b/c matrix products in () have small dim
			value[[i]][[p]] <- (t(P)%*%u[[i]][[p]])%*%diag(sigma[[i]][[p]])%*%(t(v[[i]][[p]])%*%D)
		}
	}
}
e_t_recon_v<-Sys.time()
e_t_recon_v-b_t_recon_v

dim(u[[1]][[1]])
dim(u[[1]][[2]])
length(sigma[[1]][[1]])
length(sigma[[1]][[2]])
dim(v[[1]][[1]])
dim(v[[1]][[2]])
b_t_recon_img<-Sys.time()
#Reconstruct all images
apprs <- vector('list',I)

for (i in 1:I)
{
		for (p in 1:length(rs[[1]]))
		{
			apprs[[i]][[p]] <- matrix(0,nrow=X, ncol=length(rs[[1]]))
			if (length(sigma[[i]][[p]])<2){
				#here block matrix multiplication can be used if dim(p)[1] is too big
				apprs[[i]][[p]] <- P%*%((t(P)%*%u[[i]][[p]])%*%(sigma[[i]][[p]]%*%diag(1))%*%(t(v[[i]][[p]])%*%D)%*%t(D))
		}
		else
		{
			#here block matrix multiplication can be used if dim(p)[1] is too big
			apprs[[i]][[p]] <- P%*%((t(P)%*%u[[i]][[p]])%*%diag(sigma[[i]][[p]])%*%(t(v[[i]][[p]])%*%D)%*%t(D))
		}
}
}
e_t_recon_img<-Sys.time()
e_t_recon_img-b_t_recon_img
object.size(apprs)
#save(apprs,file="/home/kangk4/dissertation/data/pvd_apprs.rda")

b_t_dif<-Sys.time()
vrs<-ltov(rs)
vapprs<-ltov(apprs)
difvec<-abs(vapprs-vrs)
e_t_dif<-Sys.time()
e_t_dif-b_t_dif
#save(difvec,file="/home/kangk4/dissertation/data/pvd_difvec.rda")
save(difvec,file="/home/kangk4/dissertation/data/pvd_difvec_x.rda")
length(difvec)
sum(difvec)/length(difvec)
tailpvd0.2<- difvec[difvec>0.2]
length(tailpvd0.2)
tailpvd0.15<- difvec[difvec>0.15]
length(tailpvd0.15)

b_t_mse<-Sys.time()
mse_pvd<-sum(difvec^2)/length(difvec)
e_t_mse<-Sys.time()
e_t_mse-b_t_mse
mse_pvd
max(difvec)
quantile(difvec,c(0.75,0.80,0.85,0.90,0.95,0.99,0.9999))

t_end<-Sys.time()
t_end-t_start
