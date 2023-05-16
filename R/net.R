#' Function to perform network analysis on a distance/dissimilarity matrix of a microbial community.
#' 
#' 
#' @param df Long dataframe such as that returned from \link{\code{asvDist}}.
#' @param metadata Optional dataframe of metadata to attach to nodes. Metadata should belong to data at the node level, so taxonomy for ASV nodes or genotype per sample nodes, but not genotype per ASV since an ASV is present in many samples.
#' @param edge Column name of df to use for edge weighting. Typically this is the same as \code{method} from \code{asvDist} used to make df.
#' @param thresh Optional output from \link{\code{thresh}}. If provided then P values of changepoints will be added as node metadata.
#' @param metadata_join Column name of metadata that identifies nodes. Defaults to "asv".
#' @param thresh_join Column name of thresh that identifies nodes. Defaults to "asv" and there are no cases currently where we expect it to be a good idea to change this.
#' @keywords network, changepoint
#' @import igraph
#' @import data.table
#' @return A named list with three elements:
#' "Nodes" is a dataframe of nodes and their metadata
#' "Edges" is a dataframe of edges connecting nodes.
#' "graph" is the igraph object used to generate the dataframes.
#' .
#' 
#' @examples 
#' 
#' # a<-qc(); b<-cal(a); c<-thresh(b); d<-asvDist(a) ; e<-net(d, thresh = c)
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/threshOutput.rdata"))
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/asvDist_output.rdata"))
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/taxa_rdp.rdata"))
#' taxa<-as.data.frame(taxa_rdp)
#' taxa$asv<-rownames(taxa_rdp)
#' rownames(taxa)<-NULL
#' net_data<-asvNet(asv, taxa, edge="spearman", thresh=threshMods)
#' 
#' @export

asvNet<-function(df, metadata=NULL, edge=NULL, thresh=NULL, metadata_join="asv", thresh_join="asv"){
  if(is.data.frame(df)){
    g<-igraph::graph_from_data_frame(df, directed=F)
  } else if (is.matrix(df)){
    g<-igraph::graph_from_adjacency_matrix(df, "undirected")
  }
  nd<-as.data.frame(igraph::layout.auto(g))
  eg<-igraph::get.data.frame(g)
  #* link metadata to nodes
  if(!is.null(metadata)){
    nd[[metadata_join]]<-igraph::as_data_frame(g, "vertices")$name
    nd<-merge(nd, metadata, by=metadata_join) 
  }
  #* link thresh to nodes
  if(!is.null(thresh)){
    agThresh<-aggregate(p.value. ~ asv+phenotype+model, thresh[thresh$Source!="(Intercept)",], FUN=identity)
    subThresh<-data.table::dcast(data.table::as.data.table(agThresh), asv ~ phenotype+model, value.var='p.value.')
    colnames(subThresh)[2:ncol(subThresh)]<-paste0(colnames(subThresh)[2:ncol(subThresh)], "_p")
    nd<-merge(nd, subThresh, by=thresh_join) 
  }
  #* Calculate network metrics
  nd$betweenness<-igraph::betweenness(g)
  nd$degree<-igraph::degree(g)
  if(!is.null(edge)){
    igraph::E(g)$weight<-eg[[edge]]+0.1
    nd$strength<-igraph::strength(g)
  }
  nd$harmonic_centrality<-igraph::harmonic_centrality(g)
  nd$eigen_centrality<-igraph::eigen_centrality(g)[[1]]
  nd$authority_score<-igraph::authority_score(g)[[1]]
  nd$page_rank<-igraph::page_rank(g)[[1]]
  #* add coordinates for plotting edges
  eg$from.x <- nd$V1[match(eg$from, nd[[metadata_join]])]  
  eg$from.y <- nd$V2[match(eg$from, nd[[metadata_join]])]
  eg$to.x <- nd$V1[match(eg$to, nd[[metadata_join]])] 
  eg$to.y <- nd$V2[match(eg$to, nd[[metadata_join]])]
  return(list("nodes"=nd, "edges"=eg, "graph"=g))
}




