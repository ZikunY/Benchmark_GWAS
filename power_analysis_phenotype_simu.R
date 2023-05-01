library("mvtnorm")
setwd('/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office/vcf')
fam<-read.table('hapnest_20_prune_chunk_two_way_1_phased_quick.fam')
dosage<-read.table('hapnest_20_prune_chunk_two_way_1_phased_quick.anc0.dosage.txt',
                   header=T)
dosage1<-read.table('hapnest_20_prune_chunk_two_way_1_phased_quick.anc1.dosage.txt',
                   header=T)
snp.info<-dosage[,1:5]
dosage<-dosage[,-(1:5)]
dosage1<-dosage1[,-(1:5)]
maf0<-colSums(t(dosage))/(ncol(dosage)*2)
maf1<-colSums(t(dosage1))/(ncol(dosage)*2)

snp.list<-list()
snp.list[[1]]<-which((0.00<maf0&maf0<0.001))
snp.list[[2]]<-which((0.001<maf0&maf0<0.01))
snp.list[[3]]<-which((0.01<maf0&maf0<0.05))
snp.list[[4]]<-which((0.05<maf0))

ultra.rare.snp<-which((0.00<maf0&maf0<0.001))
rare.snp<-which((0.001<maf0&maf0<0.01))
uncommon.snp<-which((0.01<maf0&maf0<0.05))
common.snp<-which((0.05<maf0))

snp.index<-c('ultra.rare.snp','rare.snp','uncommon.snp','common.snp')

#rel<-read.table('chr_1_peru_plink2_clean_mac.rel')
#rel<-as.matrix(rel)
pca<-read.table('hapnest_20_prune_chunk_two_way_1.eigenvec')
used.pca<-pca[,3:5]

snp=1
tau=1
n=nrow(fam)

#beta=as.matrix(rep(log(alpha),ncol(geno)))
beta=c(1.5,2.0,2.5,3.0)

scenario<-c('balanced','imbalanced','extreme_imbalanced')
scenario.limit<-list()
scenario.limit[[1]]<-c(0.40,0.60)

scenario.limit[[2]]<-c(0.05,0.15)
scenario.limit[[3]]<-c(0.005,0.015)
setwd('/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/office/')
s=1
sc=2
d=4
snp=1
for(snp in 4){
for(d in 1:length(beta)){
  s=1
  if(!dir.exists(paste0('power_phenotype/',snp.index[snp],'/',beta[d]))){
    dir.create(paste0('power_phenotype/',snp.index[snp],'/',beta[d]),recursive = T)
  }
  ass=1
  ti=1
while(s<101){
  ass=ass+1
  if(ass>50){
      break
  }
  X=matrix(c(rnorm(n,0,1),rbinom(n,1,0.5)),ncol=2,nrow=n)
  b=as.matrix(rnorm(n))
  #b<-t(b)
  cc.ratio<-1/(40+d*1+snp*10)#1/30#1/10
  a0=log(cc.ratio)
  
  mu=as.matrix(rep(a0,n))+b+X%*%as.matrix(rep(1,2))+t(as.matrix(dosage[snp.list[[snp]][s],]))*beta[d]
  
  binary.outcome<-rbinom(n,1,exp(mu)/(1+exp(mu)))
  tractor.phe<-data.frame(IID=fam$V2,
                          y=binary.outcome,
                          AGE=X[,1],
                          SEX=X[,2],
                          PC1=pca[,3],
                          PC2=pca[,4],
                          PC3=pca[,5])
       print(sum(binary.outcome))
  
  
  if(sum(binary.outcome)>round(n*scenario.limit[[sc]][1])&sum(binary.outcome)<round(n*scenario.limit[[sc]][2])){
    print(s)
    print(sum(binary.outcome))
    
    write.table(tractor.phe,paste0('power_phenotype/',snp.index[snp],'/',beta[d],'/phe_',scenario[sc],'_',s,'.txt'),
                quote = F,row.names = F,col.names = T,sep = '\t')
    s=s+1
    ass=1
  }
  

  
}

}
}

###########Extreme case control
s=1
sc=3
d=4
snp=4
for(snp in 3:4){
for(d in 1:length(beta)){
  s=1
  if(!dir.exists(paste0('power_phenotype/',snp.index[snp],'/',beta[d]))){
    dir.create(paste0('power_phenotype/',snp.index[snp],'/',beta[d]),recursive = T)
  }
  ass=1
  ti=1
while(s<101){
  ass=ass+1
  if(ass>50){
      break
  }
  X=matrix(c(rnorm(n,0,1),rbinom(n,1,0.5)),ncol=2,nrow=n)
  b=as.matrix(rnorm(n))
  #b<-t(b)
  cc.ratio<-1/(100+d*1+snp*1)#1/30#1/10
  a0=log(cc.ratio)
  
  mu=as.matrix(rep(a0,n))+b+X%*%as.matrix(rep(1,2))+t(as.matrix(dosage[snp.list[[snp]][s],]))*beta[d]
  
  #binary.outcome<-rbinom(n,1,exp(mu)/(1+exp(mu)))
  binary.outcome<-rep(0,n)
  binary.outcome[order(exp(mu)/(1+exp(mu)),decreasing=T)[1:round(n*0.01)]]<-1
  
  
  tractor.phe<-data.frame(IID=fam$V2,
                          y=binary.outcome,
                          AGE=X[,1],
                          SEX=X[,2],
                          PC1=pca[,3],
                          PC2=pca[,4],
                          PC3=pca[,5])
       print(sum(binary.outcome))
  
  
  
  if(sum(binary.outcome)>round(n*scenario.limit[[sc]][1])&sum(binary.outcome)<round(n*scenario.limit[[sc]][2])){
    print(s)
    print(sum(binary.outcome))
    
    write.table(tractor.phe,paste0('power_phenotype/',snp.index[snp],'/',beta[d],'/phe_',scenario[sc],'_',s,'.txt'),
                quote = F,row.names = F,col.names = T,sep = '\t')
    s=s+1
    ass=1
  }
  

  
}

}
}
