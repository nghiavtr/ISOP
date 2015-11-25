#' Detect sides of models
#'
#' @param dmix.list A list of detected models
#' @param nq.vec A vector of number of components of detected models
#' @return A vector consisting of sides of models
#' @export
#' @examples
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' library("AnnotationDbi")
#' txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
#' #The txdb also can be imported from a gtf file by using package GenomicFeatures
#' data(isoformDataSample)
#' #preprocessing
#' isoformDataSample=ifelse(isoformDataSample <= 3,0,isoformDataSample)
#' isoformDataSample=isoformDataSample[which(rowSums(isoformDataSample)>0),]
#' #tranform read count dataset to log scale
#' isoformDataSample=logC(isoformDataSample)
#' #now data is ready
#' tbreak=round(sqrt(ncol(isoformDataSample)))
#' model.res=doMixtureModelMatrix(isoformLevel.data=isoformDataSample,txdb=txdb,tbreak=tbreak)
#' res=detectSide(dmix.list=model.res$dmix.list,nq.vec=unlist(model.res$nq.list))

detectSide <- 
function(dmix.list,nq.vec){
  res= lapply(c(1:length(nq.vec)), function(x){
    # find null distribution
    d0=which.min(abs(dmix.list[x][[1]][[nq.vec[x]]]$delta));
    mydelta=dmix.list[x][[1]][[nq.vec[x]]]$delta	
    if (nq.vec[x]==1) return (0); 
    if (nq.vec[x]==2) {
      if (mydelta[d0]<0.75) return(1)
      
      res=mydelta[1]*mydelta[2]
      res=ifelse(sign(res)<0,2,1)
      return(res)
    }
    
    if (nq.vec[x]==3) {
      mydelta=dmix.list[x][[1]][[nq.vec[x]]]$delta
      #remove the null distribution
      mydelta=mydelta[-d0];
      res=mydelta[1]*mydelta[2]
      res=ifelse(sign(res)<0,2,1)
      return(res)
    }
  })
  res=unlist(res)
  return(res)
}