## Inferring single-cell trajectories through critical cell identification and a greedy strategy

### 1. Environment
+ R version 4.0.3
+ magrittr
+ irlba
+ Rtsne
+ umap
+ NbClust
+ philentropy
+ igraph
+ mclust
+ princurve
+ stats

### 2. Example of scTICG
```
expData<-read.csv("GSE67310_0.1_MAGIC.csv",row.names = 1) # expData is a genes*cells matrix
expData<-t(expData) # transpose expData
obj<-createObject(expData) # create an object of trajectory inference
obj<-computePCA(obj)
# Step 1. construct a cluster-based trajectory
obj<-dimRed(obj,use_pca = F)
obj<-findClusters(obj) 
obj<-calClustersDisance(obj) 
obj<-constructMST(obj,startCluster = 1)

# Step 2. optimize cluster-based on critical cells and a greedy strategy, including identifying critical based on graph centrality algorithm,
          utilizing critical cells to refining trajectory based on a greedy strategy
obj<-optimizeLineages(obj)

# Step 3. calculate pseudotime based on principal curve algorithm
obj<-calculate_pseudotime(obj)
```

