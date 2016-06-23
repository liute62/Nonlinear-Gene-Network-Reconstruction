library(org.Mm.eg.db)
setwd("~/R Files/co/data-1.23")

co_com<-output.co$community
ov_co_1<-read.table("1 ov_co 1 .txt")
ov_co_1_id<-as.character(ov_co_1[[1]])
ov_co_2<-read.table("1 ov_co 2 .txt")
ov_co_2_id<-as.character(ov_co_2[[1]])

co_cluster_1<-vector()
for(i in 1 : length(co_com)){
  if(co_com[i] == 1){
    tmp<-as.character(result.co[[1]][i])
    tmp<-mget(tmp, org.Mm.egSYMBOL2EG,ifnotfound = NA)
    if(is.na(tmp[[1]])){
       
    }else{
      co_cluster_1<-append(co_cluster_1,tmp)
    }
  }
}


co_cluster_2<-vector()
for(i in 1 : length(co_com)){
  if(co_com[i] == 2){
    tmp<-as.character(result.co[[1]][i])
    tmp<-mget(tmp, org.Mm.egSYMBOL2EG,ifnotfound = NA)
    if(is.na(tmp[[1]])){
      
    }else{
      co_cluster_2<-append(co_cluster_2,tmp)
    }
  }
}

for(i in 1 : length(ov_co_1_id)){
  genes<-mget(ov_co_1_id[i], org.Mm.egGO2EG,ifnotfound = NA)[[1]]
  tmp<-intersect(co_cluster_1, genes)
  if(length(tmp) != 0){
    tmp2<-mget(tmp, org.Mm.egSYMBOL,ifnotfound = NA)
    tmp3<-cbind(tmp[1:length(tmp)],tmp2[1:length(tmp)])
    str<-paste("~/R Files/co/data-1.23/output/geneID",i)
    str<-paste(str,".txt")
    write.table(tmp3,str)
  }
}

for(i in 1 : length(ov_co_2_id)){
  genes<-mget(ov_co_2_id[i], org.Mm.egGO2EG,ifnotfound = NA)[[1]]
  tmp<-intersect(co_cluster_2, genes)
  if(length(tmp) != 0){
    tmp2<-mget(tmp, org.Mm.egSYMBOL,ifnotfound = NA)
    tmp3<-cbind(tmp[1:length(tmp)],tmp2[1:length(tmp)])
    str<-paste("~/R Files/co/data-1.23/output-co-2/geneID",i)
    str<-paste(str,".txt")
    write.table(tmp3,str)
  }
}



ute_com<-output.ute$community

##########################################
ov_ute_1<-read.table("1 ov_ute 1 .txt")
ov_ute_1_id<-as.character(ov_ute_1[[1]])
ute_cluster_1<-vector()
for(i in 1 : length(ute_com)){
  if(ute_com[i] == 1){
    tmp<-as.character(result.ute[[1]][i])
    tmp<-mget(tmp, org.Mm.egSYMBOL2EG,ifnotfound = NA)
    if(is.na(tmp[[1]])){
      
    }else{
      ute_cluster_1<-append(ute_cluster_1,tmp)
    }
  }
}

for(i in 1 : length(ov_ute_1_id)){
  genes<-mget(ov_ute_1_id[i], org.Mm.egGO2EG,ifnotfound = NA)[[1]]
  tmp<-intersect(ute_cluster_1, genes)
  if(length(tmp) != 0){
    tmp2<-mget(tmp, org.Mm.egSYMBOL,ifnotfound = NA)
    tmp3<-cbind(tmp[1:length(tmp)],tmp2[1:length(tmp)])
    str<-paste("~/R Files/co/data-1.23/output-ute-1/geneID",i)
    str<-paste(str,".txt")
    write.table(tmp3,str)
  }
}

##########################################
ov_ute_2<-read.table("1 ov_ute 2 .txt")
ov_ute_2_id<-as.character(ov_ute_2[[1]])

ute_cluster_2<-vector()
for(i in 1 : length(ute_com)){
  if(ute_com[i] == 2){
    tmp<-as.character(result.ute[[1]][i])
    tmp<-mget(tmp, org.Mm.egSYMBOL2EG,ifnotfound = NA)
    if(is.na(tmp[[1]])){
      
    }else{
      ute_cluster_2<-append(ute_cluster_2,tmp)
    }
  }
}

for(i in 1 : length(ov_ute_2_id)){
  genes<-mget(ov_ute_2_id[i], org.Mm.egGO2EG,ifnotfound = NA)[[1]]
  tmp<-intersect(ute_cluster_2, genes)
  if(length(tmp) != 0){
    tmp2<-mget(tmp, org.Mm.egSYMBOL,ifnotfound = NA)
    tmp3<-cbind(tmp[1:length(tmp)],tmp2[1:length(tmp)])
    str<-paste("~/R Files/co/data-1.23/output-ute-2/geneID",i)
    str<-paste(str,".txt")
    write.table(tmp3,str)
  }
}


##########################################
ov_ute_3<-read.table("1 ov_ute 3 .txt")
ov_ute_3_id<-as.character(ov_ute_3[[1]])

ute_cluster_3<-vector()
for(i in 1 : length(ute_com)){
  if(ute_com[i] == 3){
    tmp<-as.character(result.ute[[1]][i])
    tmp<-mget(tmp, org.Mm.egSYMBOL2EG,ifnotfound = NA)
    if(is.na(tmp[[1]])){
      
    }else{
      ute_cluster_3<-append(ute_cluster_3,tmp)
    }
  }
}

for(i in 1 : length(ov_ute_3_id)){
  genes<-mget(ov_ute_3_id[i], org.Mm.egGO2EG,ifnotfound = NA)[[1]]
  tmp<-intersect(ute_cluster_3, genes)
  if(length(tmp) != 0){
    tmp2<-mget(tmp, org.Mm.egSYMBOL,ifnotfound = NA)
    tmp3<-cbind(tmp[1:length(tmp)],tmp2[1:length(tmp)])
    str<-paste("~/R Files/co/data-1.23/output-ute-3/geneID",i)
    str<-paste(str,".txt")
    write.table(tmp3,str)
  }
}

##########################################
ov_ute_4<-read.table("1 ov_ute 4 .txt")
ov_ute_4_id<-as.character(ov_ute_4[[1]])

ute_cluster_4<-vector()
for(i in 1 : length(ute_com)){
  if(ute_com[i] == 4){
    tmp<-as.character(result.ute[[1]][i])
    tmp<-mget(tmp, org.Mm.egSYMBOL2EG,ifnotfound = NA)
    if(is.na(tmp[[1]])){
      
    }else{
      ute_cluster_4<-append(ute_cluster_4,tmp)
    }
  }
}

for(i in 1 : length(ov_ute_4_id)){
  genes<-mget(ov_ute_4_id[i], org.Mm.egGO2EG,ifnotfound = NA)[[1]]
  tmp<-intersect(ute_cluster_4, genes)
  if(length(tmp) != 0){
    tmp2<-mget(tmp, org.Mm.egSYMBOL,ifnotfound = NA)
    tmp3<-cbind(tmp[1:length(tmp)],tmp2[1:length(tmp)])
    str<-paste("~/R Files/co/data-1.23/output-ute-4/geneID",i)
    str<-paste(str,".txt")
    write.table(tmp3,str)
  }
}

##########################################
ov_ute_5<-read.table("1 ov_ute 5 .txt")
ov_ute_5_id<-as.character(ov_ute_5[[1]])

ute_cluster_5<-vector()
for(i in 1 : length(ute_com)){
  if(ute_com[i] == 5){
    tmp<-as.character(result.ute[[1]][i])
    tmp<-mget(tmp, org.Mm.egSYMBOL2EG,ifnotfound = NA)
    if(is.na(tmp[[1]])){
      
    }else{
      ute_cluster_5<-append(ute_cluster_5,tmp)
    }
  }
}

for(i in 1 : length(ov_ute_5_id)){
  genes<-mget(ov_ute_5_id[i], org.Mm.egGO2EG,ifnotfound = NA)[[1]]
  tmp<-intersect(ute_cluster_5, genes)
  if(length(tmp) != 0){
    tmp2<-mget(tmp, org.Mm.egSYMBOL,ifnotfound = NA)
    tmp3<-cbind(tmp[1:length(tmp)],tmp2[1:length(tmp)])
    str<-paste("~/R Files/co/data-1.23/output-ute-5/geneID",i)
    str<-paste(str,".txt")
    write.table(tmp3,str)
  }
}

