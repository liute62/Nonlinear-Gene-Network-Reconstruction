#implementation of DCOL-based kmeans clustering
#contact: kai.wang.magic@gmail.com  Kai

library(TSP)
library(doSNOW)
#library(clues)
library(dynamicTreeCut)
library(fossil)
library(mclust)
#library(phyclust, quiet = TRUE)
# to use function adjustedRand() for quality assessment

gene.specific.null<-function(array, B=200)
{
	null.mat<-matrix(0, nrow=nrow(array), ncol=B)
	l<-ncol(array)
	d.array<-array[,1:(l-1)]
	for(i in 1:B)
	{
		this.order<-sample(l, l, replace=FALSE)
		for(j in 1:(l-1)) d.array[,j]<-abs(array[,this.order[j+1]]-array[,this.order[j]])
		null.mat[,i]<-apply(d.array, 1, sum)
	}
	r<-cbind(apply(null.mat, 1, mean), apply(null.mat, 1, sd))
	return(r)
}

# data.gen<-function(n.genes, n.samples, n.grps, aver.grp.size, n.fun.types, epsilon, n.depend)
# {
#     set.seed(Sys.time())
# 	link<-function(x, type)
# 	{
# 		x<-(x-mean(x))/sd(x)
# 		if(type == 1) return(x)
# 		if(type == 2) return(sin(2*x))
# 		if(type == 3) return(x^2)
# 		if(type == 4) return(abs(x))
# 		if(type == 5) return(x^3)
# 		if(type == 6) return(atan(4*x))
# 	}
	
# 	a<-matrix(rnorm(n.genes*n.samples),ncol=n.samples)
# 	curr.count<-0
# 	g<-new("list")
# 	for(i in 1:n.grps)
# 	{
# #		this.size<-rpois(1, aver.grp.size)
# 		this.size<-aver.grp.size
# 		if(this.size < 2) this.size<-2
		
# 		this.mat<-matrix(0, nrow=this.size, ncol=n.samples)
# 		this.mat[1,]<-rnorm(n.samples)
# 		for(j in 2:this.size)
# 		{
# 			if(n.depend==0)
# 			{
# 				this.basis<-c(1, rep(0,j-2))
# 			}else{
# #				this.basis<-sample(c(1,0), j-1, replace=T, prob=c(min(1, n.depend/(j-1)), 1-min(1, n.depend/(j-1))))
# 				if(j-1 <= n.depend) 
# 				{
# 					this.basis<-rep(1, j-1)
# 				}else{
# 					this.basis<-sample(c(rep(1, n.depend), rep(0,j-1-n.depend)), j-1, replace=F)
# 				}
				
# 			}
# 			if(sum(this.basis) > 0)
# 			{
# 				x<-rep(0,n.samples)
# 				for(k in which(this.basis == 1))
# 				{
# 					x<-x+link(this.mat[k,], sample(n.fun.types,1))*runif(1,min=-1,max=1)
# 				}
# #				x[x>quantile(x, 0.95)]<-quantile(x, 0.95)
# #				x[x<quantile(x, 0.05)]<-quantile(x, 0.05)
# 				this.mat[j,]<-x
# 				this.mat[j,]<-(this.mat[j,]-mean(this.mat[j,]))/sd(this.mat[j,])
# 			}else{
# 				this.mat[j,]<-rnorm(n.samples)
# 			}
# 		}
# 		if(n.depend == 0)
# 		{
# 			this.mat[1,]<-link(this.mat[1,], sample(n.fun.types,1))
# 			this.mat[1,]<-(this.mat[1,]-mean(this.mat[1,]))/sd(this.mat[1,])
# 		}
		
# 		if(curr.count+this.size <= n.genes)
# 		{
# 			a[(curr.count+1):(curr.count+this.size),]<-this.mat
# 			g[[length(g)+1]]<-(curr.count+1):(curr.count+this.size)
# 		}
# 		curr.count<-curr.count+this.size		
# 	}
# 	a<-a+matrix(rnorm(n.genes*n.samples, sd=epsilon),ncol=n.samples)
	
# 	g2<-rep(0, nrow(a))
# 	for(i in 1:length(g)) g2[g[[i]]]<-i
	
# 	r<-new("list")
# 	r$data<-a
# 	r$grps<-g2
# 	return(r)
# }

data.gen.2<-function(n.genes, n.samples, n.grps, aver.grp.size, n.fun.types, epsilon, n.depend)
# no flipping signs
{
    set.seed(Sys.time())
    link<-function(x, type)
    {
        x<-(x-mean(x))/sd(x)
        if(type == 1) return(x)
        if(type == 2) return(sin(2*x))
        if(type == 3) return(x^2)
        if(type == 4) return(abs(x))
        if(type == 5) return(x^3)
        if(type == 6) return(atan(4*x))
    }
    
    a<-matrix(rnorm(n.genes*n.samples),ncol=n.samples)
    curr.count<-0
    g<-new("list")
    for(i in 1:n.grps)
    {
        #		this.size<-rpois(1, aver.grp.size)
        this.size<-aver.grp.size
        if(this.size < 2) this.size<-2
        
        this.mat<-matrix(0, nrow=this.size, ncol=n.samples)
        this.mat[1,]<-rnorm(n.samples)
        for(j in 2:this.size)
        {
            if(n.depend==0)
            {
                this.basis<-c(1, rep(0,j-2))
            }else{
                #				this.basis<-sample(c(1,0), j-1, replace=T, prob=c(min(1, n.depend/(j-1)), 1-min(1, n.depend/(j-1))))
                if(j-1 <= n.depend)
                {
                    this.basis<-rep(1, j-1)
                }else{
                    this.basis<-sample(c(rep(1, n.depend), rep(0,j-1-n.depend)), j-1, replace=F)
                }
                
            }
            if(sum(this.basis) > 0)
            {
                x<-rep(0,n.samples)
                for(k in which(this.basis == 1))
                {
                    x<-x+link(this.mat[k,], sample(n.fun.types,1))
                }
                #				x[x>quantile(x, 0.95)]<-quantile(x, 0.95)
                #				x[x<quantile(x, 0.05)]<-quantile(x, 0.05)
                this.mat[j,]<-x
                this.mat[j,]<-(this.mat[j,]-mean(this.mat[j,]))/sd(this.mat[j,])
            }else{
                this.mat[j,]<-rnorm(n.samples)
            }
        }
        if(n.depend == 0)
        {
            this.mat[1,]<-link(this.mat[1,], sample(n.fun.types,1))
            this.mat[1,]<-(this.mat[1,]-mean(this.mat[1,]))/sd(this.mat[1,])
        }
        
        if(curr.count+this.size <= n.genes)
        {
            a[(curr.count+1):(curr.count+this.size),]<-this.mat
            g[[length(g)+1]]<-(curr.count+1):(curr.count+this.size)
        }
        curr.count<-curr.count+this.size		
    }
    for(i in 1:nrow(a)) a[i,]<-a[i,]*runif(1, min=0.7, max=1.3)
    a<-a+matrix(rnorm(n.genes*n.samples, sd=epsilon),ncol=n.samples)
    
    g2<-rep(0, nrow(a))
    for(i in 1:length(g)) g2[g[[i]]]<-i
    
    r<-new("list")
    r$data<-a
    r$grps<-g2
    return(r)
}

calDistance<- function(path, gene){
    path= strtoi(path)
    len= length(path)
    #print(len)
    if(len<2){
        print('number of samples < 2')
        q()           
    }
    d= 0
    for(i in 2:len){
        d = d+ abs(gene[path[i]]-gene[path[i-1]])            
    }
    return( d)
}


outputCluster<- function(cluster, path, indCluster){
   
    f= file('clusters.txt')
    for(i in 1:length(cluster)){
        sprintf('cluster %d contains %d members\n', i, length(cluster[[i]]) )
        #cat(path[[i]],'\n')
        #write('%d\n',i)
        
        nrow= dim(cluster[[i]])[1]
        print()          
    }            
}
          
data.cluster.2<-function(dataset, nCluster, maxIter =100, p.max=0.2, p.min=0.05){
    
    null.mat<- gene.specific.null(dataset)  #mean and sd for each gene
    
    dims= dim(dataset)
    cat('dim:',dims,'\n')
    ngene= dims[1]
    nsample= dims[2]
    
    if (nCluster > ngene){
        sprintf('datasize %d smaller than the amount of clusters %d\n',nCluster, size)
        q()
    }
    
    indCluster<- list()
    # initialization randomly
    indCluster= sample(1:nCluster, ngene, replace=TRUE)
    path= list()
    
    for(i in 1:nCluster){
            d= dist(t(dataset[indCluster==i,]))
            tsp= TSP(d)
			set.seed(999999)
            tour= solve_TSP(tsp, method="nn")
            path[[i]]= as.integer(tour)
    }
    
    p_max= p.max           # initial p-value threshold
    p_min= p.min           # more strict p-value
    p_val= p_max
    
    iloop= 0
    unclustered=c()
    while (TRUE){
        preindCluster=indCluster
        isEmpty= rep(TRUE, nCluster)
        unclustered=c()
        
        iloop= iloop+1
        if (iloop%%5 == 1){
            cat('------iteration ',iloop,'-------\n')
        }
        
        dist.mat<-matrix(0, nrow=nrow(dataset), ncol=nCluster)
        
        for(i in 1:nCluster){
            this.data<-dataset[,path[[i]]]
            this.d<-this.data[,2:ncol(this.data)]-this.data[,1:(ncol(this.data)-1)]
            this.d<-abs(this.d)
            dist.mat[,i]<-apply(this.d,1,sum)
        }
        
        indCluster<-max.col(-dist.mat, "last")		
		p.rec<-pnorm((dist.mat[cbind(1:nrow(dist.mat), indCluster)]-null.mat[,1])/null.mat[,2], lower.tail=TRUE)
		
		indCluster[p.rec>p_val]<- -1      
        if(p_val>p_min) p_val <- p_val- (p_max-p_min)/min(50, maxIter)
        
        for(i in 1:nCluster){
            if(sum(indCluster==i)>0){
                d= dist(t(dataset[indCluster==i,]))
                tsp= TSP(d)
				set.seed(999999)
                tour= solve_TSP(tsp, method="nn")
                path[[i]]= as.integer(tour)
            }
        }

        if(sum(preindCluster == indCluster)/nrow(dataset) > 0.99 || iloop> maxIter) break
    }

    r<-new("list")
    r$cluster= indCluster
    r$p.list=p.rec
    return (r)
}


# find the best K of k-means clustering by Gap Statistic
findK <- function(dataset){
    source('clusGap.R')
    n.genes<-1000
    n.samples<-c(100)
    n.grps<-5
    aver.grp.size<-200
    n.fun.types<-c(4)
    epsilon=c(0.2)
    n.depend= 2

    dataset=data.gen(n.genes, n.samples, n.grps, aver.grp.size, n.fun.types, epsilon, n.depend)

    if (0 %in% unique(dataset$grps)){
        ngrp= length(unique(dataset$grps))-1
    }else{
        ngrp= length(unique(dataset$grps))
    }
    print(ngrp)

    gskmn <- clusGap(dataset$data, FUN = data.cluster, grps=dataset$grps, K.max = ngrp+3, B = 60)
    gskmn
    plot(gskmn, main = "clusGap(., FUN = kmeans, n.start=20, B= 60)")
}

test <- function(){
	sink("log.txt", split="TRUE")
	f= file('result.txt','a')
	
	for (n.depend in 0:2){
		for (i in 2:7){             #epsilon ranges from 0.2 to 0.7
		n.samples<-c(100)
	    n.grps<-5                   #total number of genes equal to n.grps * aver.grp.size
	    aver.grp.size<-200
	    n.genes<- n.grps * aver.grp.size
	    
	    n.fun.types<-c(4)
	    epsilon=c(0.1*i)
	     
	
	    dataset=data.gen(n.genes, n.samples, n.grps, aver.grp.size, n.fun.types, epsilon, n.depend)
		cat('epsilon: ',epsilon, '\n')
		cat('depend: ',n.depend, '\n')
		writeLines(sprintf("epsilon %f, depend %f\n", epsilon, n.depend),f)
		flush(f)
		
		dr<- c()
		kr<- c()
		hr<- c()
		#change number of iterations for each parameter setting here
		for(j in 1:20){
			write(sprintf("iteration %d", j), stderr())

			if (0 %in% unique(dataset$grps)){
		        ngrp= length(unique(dataset$grps))-1
		    }else{
		        ngrp= length(unique(dataset$grps))
		    }
		    print(ngrp)
		    cat ("truth\n",dataset$grps,'\n')
		
		    results= data.cluster(dataset$data, ngrp, dataset$grps, 10 )$cluster
		    #cat ('DCOL\n', '\n') #results,
		    write("DCOL", stderr())
		    print(table (dataset$grps, results))
		    #print(RRand(dataset$grps +1 , results+1))
		    dr<- append(dr, rand.index(dataset$grps +1 , results+1))
			print(dr[length(dr)])
			
						
            #cat("kmeans\n")
		    write("kmeans", stderr())
		    kresult= kmeans(dataset$data, ngrp)
		    print(table (dataset$grps, kresult$cluster))
		    #print(RRand(dataset$grps +1 , kresult$cluster))
		    kr<- append(kr, rand.index(dataset$grps +1 , kresult$cluster))
		    print(kr[length(kr)])
			
			
		    source('GDHC.R')
		    #cat("GDHC\n")
		    write("GDHC", stderr())
		    hc= gdhc(dataset$data)
		    hgrps= cutreeDynamic(hc, minClusterSize=50)
		    #print(hgrps)
		    print(table (dataset$grps, hgrps))
		    #print(RRand(dataset$grps + 1, hgrps)) 
		    hr<- append(hr, rand.index(dataset$grps + 1, hgrps))
		    print(hr[length(hr)])	
		}
        #last column is the mean of all runs
	   	write.table(cbind( rbind(dr, kr, hr), rowMeans(rbind(dr, kr, hr))),f)
	   	flush(f)
	}
	}
	
	
	close(f)
    sink()

}



