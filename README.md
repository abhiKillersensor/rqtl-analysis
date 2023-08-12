# rqtl-analysis
#Load required package
library("qtl")
boltbdh23_root_Fwt <- read.cross("csv", ".", "Root_Fwt_R23.csv",genotypes=c("A","B"))
class(boltbdh23_root_Fwt)[1] <- "riself"
#Data checking
summary(boltbdh23_root_Fwt)
boltbdh23_root_Fwt <-jittermap(boltbdh23_root_Fwt)
summary(boltbdh23_root_Fwt) #To view data summary:
#for plot
plot(boltbdh23_root_Fwt)
plotMissing(boltbdh23_root_Fwt)
plot.map(boltbdh23_root_Fwt)
plotGeno(boltbdh23_root_Fwt)
#Explore scatterplots of phenotype vs. phenotype to locate data that may be erroneous:
pairs(jitter(as.matrix (boltbdh23_root_Fwt$pheno)), cex=0.6, las=1)
#Missing genotype information
plotInfo(boltbdh23_root_Fwt,col=c("blue")) #lower graphs
hist(nmissing(boltbdh23_root_Fwt, what="mar"), breaks=50)

#For missing values:
boltbdh23_root_Fwt$pheno[boltbdh23_root_Fwt$pheno == 0] <- NA
#Segregation distortion
gt <-geno.table(boltbdh23_root_Fwt)
gt<- as.data.frame(gt)
p_value <- gt$P.value
gt<-gt[!is.na(p_value) & p_value <0.05,]

#Compare individual genotypes
cg <- comparegeno(boltbdh23_root_Fwt)
Comparegeno_file<-cg
hist(cg, breaks=200,ylim = c(0,500), xlab="proportion of identical genotypes")
rug(cg)

#individuals with very common genotype
which(cg >0.9, arr.ind=TRUE)
#CHECK MARKER ORDER (Pairwise recombination fractions)
boltbdh23_root_Fwt <- est.rf(boltbdh23_root_Fwt)
plotRF(boltbdh23_root_Fwt,main="1")
checkAlleles(boltbdh23_root_Fwt)
#estimate map
estimate_map<-est.map(boltbdh23_root_Fwt,error.prob = 0.001,verbose=TRUE)
plot(estimate_map,main="Estimate Map")
plot(boltbdh23_root_Fwt,estimate_map)
#Examining number of crossovers
nxo <- countXO(boltbdh23_root_Fwt)
#The observed number of crossovers in boltbdh23 can be observed:
plot(nxo, ylab="No. crossovers")
nxo[nxo>100]
nxo[nxo>50]
nxo[nxo>25]
mean(nxo[1:50])
mean(nxo[-(1:50)])
countXO(boltbdh23_root_Fwt, bychr=TRUE)[5,]
countXO(boltbdh23_root_Fwt, bychr=TRUE)[31,]
countXO(boltbdh23_root_Fwt, bychr=TRUE)[46,]

#in table form
z<- plotInfo(boltbdh23_root_Fwt, step=0,alternate.chrid = TRUE)

#For chr 1 and 5:
z[ z[,1]=="C1",]
z[ z[,1]=="C5",]
boltbdh23_root_Fwt <- calc.genoprob(boltbdh23_root_Fwt, step=1, error.prob=0.001)

#Looking for qtls

#scanone em method
out1.em1<-scanone(boltbdh23_root_Fwt,method = "em")
plot(out1.em1, main=" EM LOD Fwt",ylim = c(0,1.5),alternate.chrid = TRUE)
summary(out1.em1)#single largest from each chromosome

# scanone np method non parametric
out.np<-scanone(boltbdh23_root_Fwt,mode = "np")
plot(out.np,main="np LOD Fwt",alternate.chrid=TRUE,ylim = c(0,1.5))
summary(out.np)
max(out.np$lod)

# Perform QTL analysis using Composite Interval Mapping (CIM)
out.cim <- cim(boltbdh23_root_Fwt)
plot(out.cim,main="CIM",alternate.chrid=TRUE)
summary(out.cim)

#Permutation test
operm1000 <- scanone(boltbdh23_root_Fwt,n.perm=1000,mode="np",verbose = TRUE) #with 1000 permutation 
plot(operm1000)
summary(operm1000)
summary(operm1000,alpha=c(0.05,0.10,0.20))
