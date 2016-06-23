#prefix<-"~/Desktop/"

#setwd(paste(prefix, "Dropbox/1-work/old/K profile/code/code", sep=""))
setwd("/Users/liuhaodong1/simu")
source("clusGap.R")
source("DCOL_Kmeans.R")
source("GDHC.R")

#setwd(paste(prefix, "Dropbox/1-work/old/K profile/simu", sep=""))

#library(doSNOW)
#n.node=16
#cl<-makeSOCKcluster(rep("localhost",n.node))
#registerDoSNOW(cl)
#clusterEvalQ(cl, prefix<-"C:/Users/tyu8/")
#clusterEvalQ(cl, prefix<-"~/Desktop/")
#clusterEvalQ(cl,source(paste(prefix, "Dropbox/1-work/old/K profile/code/code/clusGap.R", sep="")))
#clusterEvalQ(cl,source(paste(prefix, "Dropbox/1-work/old/K profile/code/code/DCOL_Kmeans.R", sep="")))
#clusterEvalQ(cl,source(paste(prefix, "Dropbox/1-work/old/K profile/code/code/GDHC.R", sep="")))


############# the function to conduct simulation for one scenario
# epsilon: noise level
# n.depend: data generation dependence structure. can be 0, 1, 2
# n.samples: the number of columns of the matrix
# n.grps: the number of hidden clusters
# aver.grp.size: averge number of genes in a cluster
# n.fun.types: number of function types to use
# total.n: number of simulations
# min.size.ratio: to be passed to the dynamic tree cutting

one.simu<-function(epsilon, n.depend, n.samples=100, n.grps=10, aver.grp.size=100, n.fun.types=4, total.n=20, min.size.ratio)
{
  r<-rep(0,4)
	for(N in 1:total.n)
	{
	    n.genes<- (1+n.grps) * aver.grp.size
	    dataset=data.gen(n.genes, n.samples, n.grps, aver.grp.size, n.fun.types, epsilon, n.depend)
			if (0 %in% unique(dataset$grps))
      {
		      ngrp = length(unique(dataset$grps))-1
		  }
      else
      {
		      ngrp = length(unique(dataset$grps))
		  }
	 results<-data.cluster.2(dataset$data, ngrp, 100)$cluster
   kresult= kmeans(dataset$data, ngrp)    
   print("kresult")
   print(kresult)
    processable<-"haha"
    hc= nlhc(dataset$data)
    hgrps<-cutreeDynamic(hc, minClusterSize= aver.grp.size*min.size.ratio)
    processable<-try(hgrps2<-cutree(hc, h=-2))

    if(substr(processable,1,5)=="Error")
    {
        r<-rbind(r,c(
        adjustedRandIndex(dataset$grps +1 , results+1),
        adjustedRandIndex(dataset$grps + 1, hgrps),
        NA,
        adjustedRandIndex(dataset$grps +1 , kresult$cluster)))
	}else{
    	r<-rbind(r,c(
        adjustedRandIndex(dataset$grps +1 , results+1),
        adjustedRandIndex(dataset$grps + 1, hgrps),
        adjustedRandIndex(dataset$grps + 1, hgrps2),
        adjustedRandIndex(dataset$grps +1 , kresult$cluster)))
	}
    }
  print("end")
	r<-r[-1,]
	r
}


############# actual simulation

all.epsilon<-c(1.2, 1.6, 2)
all.n.depend<-c(0,1,2)
all.grp.size<-c(100)
min.size.ratio<-c(0.5, 0.7, 0.9, 1)
n.grps<-c(10, 20)
combos<-expand.grid(all.epsilon, all.n.depend,all.grp.size,min.size.ratio, n.grps)
combos
#foreach(m=nrow(combos):1) %dopar%
for(m in 1:nrow(combos))
{
    rec<-one.simu(epsilon=combos[m,1], n.depend=combos[m,2], aver.grp.size=combos[m,3], min.size.ratio=combos[m,4], n.grps=combos[m,5], n.fun.types=4, total.n=10)
    print(rec)
    save(rec, file=paste("simu epsilon", combos[m,1], "n.depend", combos[m,2],"grp size",combos[m,3], "min size ratio", combos[m,4], "n grps", combos[m,5], ".bin"))
}
    #}


############# plotting results
par(mfcol=c(2,3))
for(grp.size in c(100))
{
    for(n.depend in c(0,1,2))
    {
	  for(n.grps in c(10, 20))
	  {
        r<-rep(NA,8)
        for(epsilon in c(0.2, 0.6, 0.8, 1.2, 1.6))
        {
		  for(min.size.ratio in c(0.5, 0.7, 0.9, 1))
		  {
            this.name<-paste("linear simu epsilon", epsilon, "n.depend", n.depend,"grp size",grp.size, "min size ratio", min.size.ratio, "n grps", n.grps, ".bin")
            #load(this.name)
            load("simu epsilon 1.2 n.depend 2 grp size 10 min size ratio 1 n grps 20.bin")
			if(min.size.ratio == 0.5){
				this<-apply(rec,2,mean)[c(1,4,2)]
			}else{
				this<-c(this, apply(rec,2,mean)[2])
			}
		  }
		   r<-rbind(r, c(epsilon, this, apply(rec,2,mean)[3]))

        }	
		r<-r[-1,]
		plot(r[,1], r[,2], xlab="noise SD/signal SD", ylab="adjusted Rand Index", type="l", col="red",ylim=c(0,1),lwd=2, main=paste("# clusters", n.grps,", n.depend", n.depend))
		lines(r[,1], r[,3],col="green",lwd=2, lty=1)

		lines(r[,1], r[,4],col="blue",lwd=2, lty=1)
		lines(r[,1], r[,5],col="blue",lwd=2, lty=2)
		lines(r[,1], r[,6],col="blue",lwd=2, lty=3)
		lines(r[,1], r[,7],col="blue",lwd=2, lty=4)
		#lines(r[,1], r[,8],col="cyan",lwd=2, lty=1)
	  }
    }
}

