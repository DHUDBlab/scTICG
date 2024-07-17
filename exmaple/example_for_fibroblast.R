
expData<-read.csv("GSE67310_0.1_MAGIC.csv",row.names = 1)
expData<-t(expData)

obj<-createObject(expData)
obj<-computePCA(obj)
obj<-dimRed(obj,use_pca = F)
obj<-findClusters(obj) 
obj<-calClustersDisance(obj) 
obj<-constructMST(obj,startCluster = 1)
obj<-optimizeLineages(obj)
obj<-calculate_pseudotime(obj)