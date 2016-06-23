gene.m.name<-c("FGF10","SOX9","Tgfb2","Stat3","Notch2","PTEN","Chrna1","Aldoa","Igf1r","Kif3a","Adam10","Rbpj")
gene.m.id<-c()



for(j in 1:length(gene.m.name)){
  for(i in 1:nrow(gene.name)){
    if(as.character(gene.name[[1]][i]) == gene.m.name[j]){
      print(i)
      gene.m.id<-c(gene.m.id,i)
      next
    }
  }
}

gene.m.name.rel<-c()
for(j in 1:length(gene.m.id)){
    gene.m.name.rel<-c(gene.m.name.rel,as.character(gene.name[[1]][gene.m.id[j]]))
}
gene.m.name.rel

row.co<-row.names(result.co)
row.ute<-row.names(result.ute)
#for co
tmp<-as.integer(row.names(result.co))
co.m.index<-c()
for(j in 1 : length(gene.m.id)){
  for(i in 1 : length(tmp)){
    if(tmp[i] == gene.m.id[j]){
      print(i)
      co.m.index<-c(co.m.index,i)
      next
    }
  }
}

#for ute
tmp<-as.integer(row.names(result.ute))
ute.m.index<-c()
for(j in 1 : length(gene.m.id)){
  for(i in 1 : length(tmp)){
    if(tmp[i] == gene.m.id[j]){
      print(i)
      ute.m.index<-c(ute.m.index,i)
      next
    }
  }
}

graph.co<-output.co$graph
graph.ute<-output.ute$graph

#co
gene.aimed.rel.co<-matrix(nrow=length(co.m.index),ncol=300) #1: only co, 2: only ute, 3: both co and ute
gene.aimed.index.co<-matrix(nrow=length(co.m.index),ncol=300)
num<-1
for(i in 1 : length(co.m.index)){ #length(co.m.index) == length(ute.m.index)
  mat1<-graph.co[co.m.index[i]]
  mat2<-graph.ute[ute.m.index[i]]
  
  for(j in 1 : length(mat1)){
    
    if(mat1[j] != 0){
      real.co<-row.co[j]
      index.ute<-which(row.ute[1:length(row.ute)] == real.co)
      gene.aimed.index.co[i,num]<-real.co
      #means have a link in Co
      if(length(index.ute) == 0){
        #means Ute don't have such gene
        gene.aimed.rel.co[i,num]<- -1
      }else{
        if(mat2[index.ute] != 0){
          #means have a link both in ute and co
          gene.aimed.rel.co[i,num]<-3
        }else{
          #means only have a link in co
          gene.aimed.rel.co[i,num]<-1
        } 
      }
      num<-num+1
    }
  }
}

gene.aimed.rel.ute<-matrix(nrow=length(co.m.index),ncol=300) #1: only co, 2: only ute, 3: both co and ute
gene.aimed.index.ute<-matrix(nrow=length(co.m.index),ncol=300)
num<-1
for(i in 1 : length(co.m.index)){ #length(co.m.index) == length(ute.m.index)
  mat1<-graph.co[co.m.index[i]]
  mat2<-graph.ute[ute.m.index[i]]
  
  for(j in 1 : length(mat2)){
    
    if(mat2[j] != 0){
      real.ute<-row.ute[j]
      index.co<-which(row.co[1:length(row.co)] == real.ute)
      gene.aimed.index.ute[i,num]<-real.ute
      #means have a link in Co
      if(length(index.co) == 0){
        gene.aimed.rel.ute[i,num]<- -2
      }else{
        if(mat1[index.co] != 0){
          #means have a link both in ute and co
          gene.aimed.rel.ute[i,num]<-3
        }else{
          #means only have a link in co
          gene.aimed.rel.ute[i,num]<-2
        } 
      }
      num<-num+1
    }
  }
}


##print the result
setwd("~/R Files/co/data-1.23/graph")
#co
for(i in 1 : length(gene.m.name.rel)){
  suffix<-"co"
  file.name<-suffix
  middle<-gene.m.name.rel[i]
  file.name<-paste(file.name,middle)
  file.name<-paste(file.name,".txt")
  tmp.gene<-c()
  tmp.rel<-c()
  for(j in 1 : ncol(gene.aimed.index.co)){
    tmp<-gene.aimed.index.co[i,j]
    if(!is.na(tmp)){
      tmp.gene<-c(tmp.gene, as.character(gene.name[[1]][as.integer(tmp)]))
      tmp.rel<-c(tmp.rel,gene.aimed.rel.co[i,j])
    }
  }
  write.table(data.frame(tmp.gene,tmp.rel),file.name)
}
#ute
for(i in 1 : length(gene.m.name.rel)){
  suffix<-"ute"
  file.name<-suffix
  middle<-gene.m.name.rel[i]
  file.name<-paste(file.name,middle)
  file.name<-paste(file.name,".txt")
  tmp.gene<-c()
  tmp.rel<-c()
  for(j in 1 : ncol(gene.aimed.index.ute)){
    tmp<-gene.aimed.index.ute[i,j]
    if(!is.na(tmp)){
      tmp.gene<-c(tmp.gene, as.character(gene.name[[1]][as.integer(tmp)]))
      tmp.rel<-c(tmp.rel,gene.aimed.rel.ute[i,j])
    }
  }
  write.table(data.frame(tmp.gene,tmp.rel),file.name)
}







