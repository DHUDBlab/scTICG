suppressMessages(library(magrittr))

trajectoryInference<-setClass(
  Class = 'trajectoryInference',
  slots = c(
    expData = 'matrix',
    pca="list", 
    cell.embeddings = 'matrix',
    cell.labels = 'matrix',
    centroid.embeddings = 'matrix',
    distance_matrix = 'matrix',
    paths = 'list', # cluster-based trajectory
    pseudotime = 'data.frame',
    lineages_vertex = 'list', # reconstructed lineages_vertex
    real_order = 'data.frame'
  )
)

# create object for trajectory inference
createObject<-function(data){
  if(!is.matrix(data)) stop('data must be form of matrix!')
  object<-NULL
  object<-new(Class = 'trajectoryInference')
  object@expData<-data
  
  object
}


#' Principal component analysis
computePCA<-function(object, num.pcs=20, seed=NULL){
  data<-object@expData
  # data<-scale(data)
  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  if(ncol(data)<num.pcs){
    num.pcs<-ncol(data)-1
  }
  suppressMessages(library(irlba))
  pca <- irlba(A = data, nv = num.pcs)
  cell.embeddings <- pca$u %*% diag(pca$d)
  
  rownames(cell.embeddings)<-rownames(data)
  
  pca<-list(cell.embeddings=cell.embeddings)
  
  object@pca<-pca
  
  object
}

#' Perform dimensionality reduction using t-SNE
computeTSNE<-function(object,use_pca=FALSE, seed=NULL){
  suppressMessages(library(Rtsne))
  data<-object@expData
  n<-nrow(data)
  perplexity = (ceiling(n/100)*10)
  
  
  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  if (use_pca) {
    data<-object@pca$cell.embeddings
  }
  
  library(Rtsne)
  tsne.done <- Rtsne(data, perplexity=perplexity, check_duplicates = FALSE,dims = 2)
  m<- tsne.done$Y
  rownames(m) <- rownames(data)
  object@cell.embeddings <- m
  
  object
}

#' Dimension Reduction using UMAP
computeUMAP<-function(object,use_pca=FALSE,seed=NULL){
  if (use_pca) {
    data<-object@pca$cell.embeddings
  }else{
    data<-object@expData
  }
  suppressMessages(library(umap))
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (use_pca) {
    data<-object@pca$cell.embeddings
  }else{
    data<-object@expData
  }
  cell.embeddings<-umap(data,config=umap.defaults,method=c("naive"))$layout
  object@cell.embeddings<-cell.embeddings
  object
}


dimRed<-function(object,use_pca=F){
  # set.seed(12345)
  n<-nrow(object@expData)
  if(n<2000) {
    object<-computeTSNE(object,use_pca)
  }else{
    object<-computeUMAP(object,use_pca)
  }
  return(object)
}


#'  Identify clusters based hierarchical.clustering
findClusters<-function(object){
  cell.embeddings<-object@cell.embeddings
  n<-nrow(cell.embeddings)
  suppressMessages(library(NbClust))
  if (!is.null(object@pca$cell.embeddings)) {
    best_k<-NbClust(object@pca$cell.embeddings[,1:3], distance = "euclidean",min.nc	= 2,
                    method = "ward.D2", index = "kl")$Best.nc[1]
    
    if(best_k<3) best_k<-3
    
    res<-NbClust(cell.embeddings, distance = "euclidean",max.nc = best_k,min.nc = best_k,
                 method = "ward.D2", index = "kl")
  }else{
    res<-NbClust(cell.embeddings, distance = "euclidean",
                 method = "ward.D2", index = "ch")
  }
  clusterLabels<-res$Best.partition
  
  tmp<-cbind(cell.embeddings,clusterLabels) %>% as.data.frame()
  centroid_vertex<-aggregate(tmp[,-ncol(tmp)],list(tmp[,ncol(tmp)]),mean)[,-1] %>% as.matrix()
  
  object@cell.labels<-clusterLabels %>% as.matrix()
  object@centroid.embeddings<-centroid_vertex
  object
}



#' calculate distance between clusters
calClustersDisance<-function(object){
  suppressMessages(library(philentropy))
  
  centroid_vertex<-object@centroid.embeddings
  k<-nrow(centroid_vertex)
  distance_matrix<-philentropy::distance(centroid_vertex,method = "euclidean")
  
  rownames(distance_matrix)<-1:k
  colnames(distance_matrix)<-1:k
  object@distance_matrix<-distance_matrix
  
  object
}


#' Construct Minimum Spanning Tree based on clusters
constructMST<-function(object,startCluster=1){
  library(igraph)
  # Minimum spanning tree constructed for identifing lineages
  distance_matrix<-object@distance_matrix
  
  mstree <- ape::mst(distance_matrix)
  forest <- mstree
  # identify sub-trees 
  subtrees <- subtrees.update <- forest
  diag(subtrees) <- 1
  while(sum(subtrees.update) > 0){
    subtrees.new <- apply(subtrees,2,function(col){
      rowSums(subtrees[,as.logical(col), drop=FALSE]) > 0
    })
    subtrees.update <- subtrees.new - subtrees
    subtrees <- subtrees.new
  }
  subtrees <- unique(subtrees)
  trees <- lapply(seq_len(nrow(subtrees)),function(ri){
    colnames(forest)[subtrees[ri,]]
  })
  trees <- trees[order(vapply(trees,length,0),decreasing = TRUE)]
  ntree <- length(trees)
  
  lineages_temp<-list()
  start.clus <- as.character(startCluster)
  for(tree in trees){
    if(length(tree) == 1){
      lineages_temp[[length(lineages_temp)+1]] <- tree
      next
    }
    tree.ind <- rownames(forest) %in% tree
    tree.graph <- forest[tree.ind, tree.ind, drop = FALSE]
    degree <- rowSums(tree.graph)
    g <- graph.adjacency(tree.graph, mode="undirected")

    if(sum(start.clus %in% tree) > 0){
      starts <- start.clus[start.clus %in% tree]
      ends <- rownames(tree.graph)[
        degree == 1 & ! rownames(tree.graph) %in% starts]
      for(st in starts){
        paths <- shortest_paths(g, from = st, to = ends, 
                                mode = 'out', 
                                output = 'vpath')$vpath
        for(p in paths){
          lineages_temp[[length(lineages_temp)+1]] <- as.numeric(names(p))
        }
      }
    }
  }
  object@paths<-lineages_temp
  object
}


sample_in_cluster<-function(samples,nsample = 200){
  N<-dim(samples)[1]
  rownames(samples)<-1:N
  if(nsample>N){
    landmarks<-samples
    return(landmarks)
  }
  sample_ids<-rownames(samples)
  sampled_ids <- sample(sample_ids, size = nsample, replace = FALSE)
  
  landmarks<-samples[sampled_ids,]
  
  landmarks
}


cell_sample<-function(object,nsample = 200){
  cell.embeddings<-object@cell.embeddings
  cell.labels<-object@cell.labels
  k<-max(cell.labels)
  
  cluster_samples<-list()
  for (i in seq(k)) {
    samples<-cell.embeddings[cell.labels==i,]
    cluster_samples[i]<-list(sample_in_cluster(samples,
                                               nsample = nsample))
  }
  
  cluster_samples
}

findCriticalTransitionCells <- function(centrality,minLevel =3,maxLevel = 5){
  cellIdx<-rownames(centrality) %>% as.numeric()
  library(mclust)
  ### fit Gaussian Mixture Model for centrality
  mcl.o <- Mclust(centrality, G = seq(minLevel,maxLevel),verbose = F)
  
  levels <- mcl.o$class
  
  maxLevel<-max(levels)
  criticalCellIndexs<-cellIdx[which(levels==maxLevel)]
  criticalCellIndexs
}

#' reconstruct cell trajectory
optimizeLineages<-function(object,useSample=F,nsample = 200){
  library(igraph)
  lineages<-object@paths 
  cell.embeddings<-object@cell.embeddings
  centroid.embeddings<-object@centroid.embeddings
  nembed<-ncol(cell.embeddings)
  
  
  if(useSample){ # sample cells
    cluster_list<-cell_sample(object,nsample = nsample)
    k<-length(cluster_list)
    cell.embeddings<-matrix(0,nrow = 0,ncol = nembed)
    cell.labels<-c()
    for (i in seq(k)) {
      cell.embeddings<-rbind(cell.embeddings,cluster_list[[i]])
      cell.labels<-c(cell.labels,rep(i,nrow(cluster_list[[i]])))
    }
    
  }
  else{ 
    cell.labels<-object@cell.labels
  }
  
  coordinates<-rbind(centroid.embeddings,cell.embeddings) 
  
  k<-max(cell.labels)
  labels<-c(1:k,cell.labels) 
  
  coordinates_labels<-cbind(coordinates, labels)
  rownames(coordinates_labels)<-1:nrow(coordinates_labels)
  finalPaths<-list()
  
  tt<-list()
  cnt_t<-1
  for (i0 in seq(length(lineages))) {
    path<-c()
    lineage<-lineages[[i0]] %>% as.numeric()
    for (j0 in seq(length(lineage)-1)) {
      clu1<-lineage[j0]
      clu2<-lineage[j0+1]   
      od<-paste0(c(clu1,"_",clu2),collapse = "")
      
      if(length(which(names(tt)==od))>0){ 
        tmpPath<-tt[[od]]
      } else{
        # cell belonging to clu1 and clu2
        cell_points<-coordinates_labels[which(labels==clu1 | labels==clu2),-ncol(coordinates_labels)] 
        cell_points<-cell_points[-c(1,2),] 
        
        get_project_points<-function(points,line_start,line_end,flag = 1){ 
          points_idx<-c()
          if(flag == 0){ 
            line_direction2 <- line_start - line_end
            cos2<-t((t(points)-line_end)) %*% line_direction2/(sum(line_direction2^2))
            points_idx<-which(cos2>=0)
          }else if(flag == 2){
            line_direction1 <- line_end - line_start
            cos1<-t((t(points)-line_start)) %*% line_direction1/(sum(line_direction1^2)) 
            points_idx<-which(cos1>=0)
          }else if(flag == 1){
            line_direction1 <- line_end - line_start
            cos1<-t((t(points)-line_start)) %*% line_direction1/(sum(line_direction1^2)) 
            
            line_direction2 <- line_start - line_end
            cos2<-t((t(points)-line_end)) %*% line_direction2/(sum(line_direction2^2))
            points_idx<-which(cos1>=0 & cos2>=0)
          }
          points_idx
        }
        
        flag<-1
        if(j0 == 1){
          flag = 0
        }else if (j0==(length(lineage)-1)) {
          flag = 2
        }
        cell_ids<-get_project_points(cell_points,coordinates[clu1,],coordinates[clu2,],flag = flag)
        cell_points<-cell_points[cell_ids,]
        n<-nrow(cell_points)
        
        a<-(centroid.embeddings[clu2,]-centroid.embeddings[clu1,]) %>% as.vector()
        mat<-matrix(0,nrow = n,ncol = n)
        for (i1 in 1:n) {
          for (j1 in 1:n) {
            b<-cell_points[j1,]-cell_points[i1,]
            if(i1!=j1 && sum(a*b)/(sqrt(sum(a * a))* sqrt(sum(b * b)))>=0)
              mat[i1,j1]<-sqrt(sum(b * b))
          }
        }
        
        graph<-graph_from_adjacency_matrix(mat,weighted = T)
    
        closeness <- closeness(graph, weights=E(graph)$weight,mode="total",normalized = TRUE)
        
        centrality<-cbind(closeness)
        centrality[is.na(centrality)]<-0
        rownames(centrality)<-rownames(cell_points)
        criticalCells<-findCriticalTransitionCells(centrality)
        
        if(j0==1 || j0==length(lineage)-1){ 
          cellIdx<-rownames(cell_points) %>% as.numeric()
          delta<-0.05
          if (j0==1) { # initial cluster
            in_degree<-igraph::degree(graph,normalized=T,mode="in") %>% as.matrix() # lower in-degree
            maxIn<-min(in_degree)+delta
            final_sources<-cellIdx[which(in_degree<=maxIn)] 
            
            sources<-c(clu1)
            aa<-(coordinates[clu2,]-coordinates[clu1,]) %>% as.vector()
            for (ss in final_sources) {
              bb<-(coordinates[ss,]-coordinates[clu1,]) %>% as.vector()
              if(sum(aa*bb)/(sqrt(sum(aa * aa))* sqrt(sum(bb * bb)))<0) sources<-c(sources,ss)
            }
            if(length(sources)>1) sources<-sources[-1]
            
            targets<-clu2
            
          }else{ # terminal cluster
            out_degree<- igraph::degree(graph,normalized = T,mode="out") %>% as.matrix() # lower out-degree
            maxOut<-min(out_degree)+delta
            final_targets<-cellIdx[which(out_degree<=maxOut)]
            sources<-clu1
            
            targets<-c(clu2)
            aa<-(coordinates[clu1,]-coordinates[clu2,]) %>% as.vector()
            for (tt in final_targets) {
              bb<-(coordinates[tt,]-coordinates[clu2,]) %>% as.vector()
              if(sum(aa*bb)/(sqrt(sum(aa * aa))* sqrt(sum(bb * bb)))<0) targets<-c(targets,tt)
            }
            if(length(targets)>1) targets<-targets[-1]
          }
          
        }else{ 
          sources<-clu1
          targets<-clu2
        }
        cellIdxs<-c(sources,clu1,criticalCells,targets)
        cellIdxs<-cellIdxs[!duplicated(cellIdxs)] 
        
        if (length(cellIdxs)==2) {
          tmpPath<-cellIdxs
        }else{
          cells_coordinates<-coordinates_labels[cellIdxs,-ncol(coordinates_labels)]
          cell_dist<-dist(cells_coordinates) %>% as.matrix()
          nn<-nrow(cells_coordinates)
          rownames(cell_dist)<-1:nn
          colnames(cell_dist)<-1:nn
          
          for (i2 in 1:nn) {
            for (j2 in 1:nn) {
              b<-(cells_coordinates[j2,]-cells_coordinates[i2,]) %>% as.vector()
              if(cell_dist[i2,j2]!=0 && sum(a*b)/(sqrt(sum(a * a))* sqrt(sum(b * b)))<=0) cell_dist[i2,j2]<-0
            }
          }
          cell_dist<-(cell_dist-cell_dist[which.min(cell_dist)])/(cell_dist[which.max(cell_dist)]-cell_dist[which.min(cell_dist)])
          cell_dist[cell_dist==0]<-Inf        
          sourceIdxs<-as.numeric(which(cellIdxs %in% sources))
          targetIdxs<-as.numeric(which(cellIdxs %in% targets))
          
          
          tmpPaths<-list()
          for (ii in sourceIdxs) {
            if(all(is.infinite(cell_dist[ii,])))  break
            
            tmp_dist<-cell_dist 
            tmpPath<-c()
            tmpPath<-c(tmpPath,ii) 
            removeIds<-c()
            while(!any((tmpPath %in% targetIdxs)) && !all(is.infinite(tmp_dist))){  
              if(length(tmpPath)==1 && all(is.infinite(tmp_dist[ii,]))){
                break
              }
              
              currentIdx<-tmpPath[length(tmpPath)] 
              if (all(tmp_dist[currentIdx,]==Inf)) {
                tmpPath<-tmpPath[-length(tmpPath)] 
                removeIds<-c(removeIds,currentIdx)
                
              }else{ 
                nextIdx<-which.min(tmp_dist[currentIdx,]) %>% as.integer() 
                tmp_dist[currentIdx,nextIdx]<-Inf
                if (!(nextIdx %in% removeIds) && !(nextIdx %in% tmpPath)) 
                  tmpPath<-c(tmpPath,nextIdx) 
              }
            }
            if(any((tmpPath %in% targetIdxs))){
              tmpPaths[[ii]]<-cellIdxs[tmpPath]
            }
          }
          tmpPath<-c()
          if(length(tmpPaths)==0){
            tmpPath<-c(clu1,clu2)
          }else{
            for (p in 1:length(tmpPaths))
              if (length(tmpPath)<length(tmpPaths[[p]])) tmpPath<-tmpPaths[[p]]
          }
        }
        
        tt[cnt_t]<-list(tmpPath)
        names(tt)[cnt_t]<-paste0(c(clu1,"_",clu2),collapse = '')
        cnt_t<-cnt_t+1
      }
      path<-c(path,tmpPath)
      
    }
    path<-path[!duplicated(path)] %>% as.numeric()
    finalPaths[i0]<-list(path)
  }
  
  lineages_vertex<-list()
  for (path in 1:length(finalPaths)) {
    lineages_vertex[path]<-list(coordinates_labels[finalPaths[[path]],-ncol(coordinates_labels)])
  }
  object@lineages_vertex<-lineages_vertex
  
  object
  
}


calculate_pseudotime<-function(object,isNormalized=TRUE){
  suppressMessages(library(princurve))
  suppressMessages(library(stats))
  library(slingshot)
  library(graphics)
  
  cell.embeddings<-object@cell.embeddings
  cell.labels<-object@cell.labels
  embeddings_labels<-cbind(cell.embeddings,cell.labels) 
  n_embeddings<-ncol(cell.embeddings)
  n_cells<-nrow(embeddings_labels)
  
  paths<-object@paths
  lineages_vertex<-object@lineages_vertex
  
  trajectories<-list() 
  cnt<-1
  pseudotimes<-matrix(nrow =n_cells,ncol = 0)
  for (i in seq(length(paths))) {
    lineage<-as.numeric(paths[[i]])
    project_vertix<-matrix(nrow = 0,ncol = n_embeddings) 
    cell_order<-c()
    if (length(lineage)<=1) {
      break
    }      
    lineage_vertix<-lineages_vertex[[i]] 
    
    
    for (j in seq(length(lineage))) {
      cluster_idx<-lineage[j]
      project_vertix<-rbind(project_vertix,matrix(embeddings_labels[embeddings_labels[,ncol(embeddings_labels)]==cluster_idx,1:n_embeddings],ncol=n_embeddings))
      cell_order<-c(cell_order,which(embeddings_labels[,ncol(embeddings_labels)]==cluster_idx))
    }
    if(nrow(lineage_vertix)<=4){
      s<-matrix(lineage_vertix,ncol = n_embeddings)
    }else{
      tmp_fit<-principal_curve(matrix(lineage_vertix,ncol =n_embeddings),start = matrix(lineage_vertix,ncol = n_embeddings))
      s<-tmp_fit$s[tmp_fit$ord, ,drop=F]
    }
    fit<-project_to_curve(matrix(project_vertix,ncol = n_embeddings),s)
    trajectories[[cnt]]<-fit
    cnt<-cnt+1
    lambda_save<-as.matrix(fit$lambda) 
    lambda_save<-cbind(lambda_save, as.matrix(cell_order))
    lambda_fit <- matrix(0, nrow = n_cells)
    lambda_fit[lambda_save[, 2], 1] <- lambda_save[, 1]
    if(isNormalized==T){
      max_lambda<-max(lambda_fit)
      min_lambda<-min(lambda_fit)
      lambda_fit<-(lambda_fit-min_lambda)/(max_lambda-min_lambda)
    }
    pseudotimes<-cbind(pseudotimes,lambda_fit)
  }
  
  
  pseudotime <- as.matrix(apply(pseudotimes, 1, max))
  pseudotime<-data.frame(cell=rownames(embeddings_labels),pseudotime)
  
  object@pseudotime<-pseudotime
  
  library(ggplot2)
  themes <- function(base_size=13, base_family="sans serif") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
      + theme(plot.title = element_text(size = rel(1.2)),
              text = element_text(),
              panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.border = element_rect(colour = NA),
              axis.title = element_text(size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text.y = element_text(),
              axis.text.x = element_text(),
              axis.line = element_line(colour="black", size = 1),
              axis.ticks.x = element_line(),
              axis.ticks.y = element_line(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "top",
              legend.direction = "horizontal",
              legend.key.size= unit(0.5, "cm"), 
              legend.margin = margin(0.5,0,0,0,unit = "cm"),
              legend.title = element_text(colour = "black",size = 12,vjust = 1),
              plot.margin=unit(c(3.5,30,3.5,30),"mm"),
              strip.background = element_blank(),

              strip.text = element_text(face="bold", size = 12)
      ))
    
  }
  data<-data.frame(cell.embeddings[,1:2],pseudotime[,2])
  colnames(data)<-c("X","Y","Pseudotime")
  
  tt<-matrix(0,nrow = 0,ncol = 3)
  for(ii in seq(length(lineages_vertex))){
    fit <- trajectories[[ii]]
    df_line <- data.frame(fit$s[fit$ord,1:2],ii)
    colnames(df_line)<-c("X","Y","group")
    tt<-rbind(tt,df_line)
  }
  
  n_points<-nrow(data)
  p_size = -0.58*log(n_points)+6.77
  
  g<-ggplot()  + geom_point(data = data, aes(X, Y, color = Pseudotime), size = p_size) +
    labs(x="tSNE1",y="tSNE2")+
    guides(fill = guide_legend(title = "Legend Title"))+
    theme_classic()+
    theme(
      legend.position = "right", 
      legend.direction = "vertical",
      legend.key.size= unit(0.5, "cm"),
      legend.key.height =  unit(0.9,"cm"),
      legend.title = element_text(colour = "black",face="bold",
                                  margin = margin(b = 3)),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.line = element_line(linewidth = 0.7),
      axis.title = element_text(face = "bold",size=12),
      axis.ticks = element_blank()
    ) + 
    scale_color_gradient2(limits = c(0,1), breaks = seq(0,1,by=0.2),
                          midpoint=0.5, low="#0c027b", mid="#c03984",high = "#f4ec28") + 
    geom_path(data = tt, aes(x = X, y = Y, group = group),linewidth = .7)
  
  print(g)

  object
}



