args=commandArgs(trailingOnly=TRUE)
d<-as.numeric(args[1])
library("mvtnorm",lib="/mnt/mfs/hgrcgrid/homes/zy2412/R_libs/")
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


snp.list<-which((0.1<maf1&maf1<.2&0.1<maf0&maf0<.2))
snp.index<-c('ultra.rare.snp','rare.snp','uncommon.snp','common.snp')


pca<-read.table('hapnest_20_prune_chunk_two_way_1.eigenvec')

snp=1
tau=1
n=nrow(fam)


beta0=seq(-0.5,0.65,0.05)
beta1=c(0.15)

scenario<-c('balanced','imbalanced','extreme_imbalanced')
scenario.limit<-c(0.2,0.3)
setwd('/mnt/mfs/hgrcgrid/shared/GT_ADMIX/Zikun_peru_project/HAPNEST/HetLanc/')
  snp=4
  sc=1
  s=1
  if(!dir.exists(paste0('power_phenotype/',snp.index[snp],'/',d))){
    dir.create(paste0('power_phenotype/',snp.index[snp],'/',d),recursive = T)
  }
  ass=1
  
while(sc<101){
  ass=ass+1
  if(ass>50){
      break
  }
  X=matrix(c(rnorm(n,0,1),rbinom(n,1,0.5)),ncol=2,nrow=n)
  b=as.matrix(rnorm(n))
  #b<-t(b)
 # cc.ratio<-1/(40+d*1+snp*10)#1/30#1/10
 # a0=log(cc.ratio)
  
  a0=log(1/5)-mean(as.numeric(dosage[snp.list[s],]))*beta0[d]-mean(as.numeric(dosage1[snp.list[s],]))*beta1- sum(colMeans(X))

  
  mu=as.matrix(rep(a0,n))+b+X%*%as.matrix(rep(1,2))+t(as.matrix(dosage[snp.list[s],]))*beta0[d]+t(as.matrix(dosage1[snp.list[s],]))*beta1
  
  binary.outcome<-rbinom(n,1,exp(mu)/(1+exp(mu)))
  tractor.phe<-data.frame(IID=fam$V2,
                          y=binary.outcome,
                          AGE=X[,1],
                          SEX=X[,2],
                          PC1=pca[,2],
                          PC2=pca[,3],
                          PC3=pca[,4])
       print(sum(binary.outcome))
  
  
  if(sum(binary.outcome)>round(n*scenario.limit[1])&sum(binary.outcome)<round(n*scenario.limit[2])){
    print(s)
    print(sum(binary.outcome))
    
    write.table(tractor.phe,paste0('power_phenotype/',snp.index[snp],'/',d,'/phe_',scenario[1],'_',sc,'.txt'),
                quote = F,row.names = F,col.names = T,sep = '\t')
    s=s+1
    sc=sc+1
    ass=1
  }
  if(s==51){
      s=1
  }
  
}



