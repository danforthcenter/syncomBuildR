#' Function to perform network analysis on a distance/dissimilarity matrix of a microbial community.
#' 
#' I should inform this function using emd -> net from pcvr.
#' That pipeline is probably really similar to what I want to have happen here.
#' In that case this should be ready to take either a dist object or a dataframe or a matrix similar to dist object
#' That way I can support all of vegdist.
#' 
#' So I need two functions here, asvDist and net
#' 
#' 
#' 
#' 
#' @param data 
#' @keywords changepoint, threshold, regression, phenotype
#' @import chngpt
#' @return A dataframe summarizing changepoint models for individual ASVs vs phenotypes.
#' 
#' @examples 
#' 
#' # a<-qc(); b<-cal(a); c<-thresh(b); d<-asvDist(a) ; e<-net(d, thresh = c)
#' print(load("/home/jsumner/Desktop/stargate/SINC/sincUtils/syncomBuilder/threshOutput.rdata"))
#' 
#' 
#' 
#' @export



