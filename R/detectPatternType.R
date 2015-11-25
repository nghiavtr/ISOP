#' Detect a type of pattern for each pair
#'
#' @param dmix.list A list of detected models
#' @param side.nqs A list of type of sides of models
#' @param nq.list A list of number of components of detected models
#' @param isoformLevel.data A data matrix of isoform expression
#' @param geneLevel.data A data matrix of gene expression
#' @return A list consisting of types of patterns and merged types of patterns
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
#' res=detectPatternType(dmix.list=model.res$dmix.list,
#'     nq.list=model.res$nq.list,isoformLevel.data=isoformDataSample)
#' pie(table(res$pattern.type),col=colors()[c(12,11,132,131,137,136)])
detectPatternType <- 
function(dmix.list,side.nqs=NULL,nq.list,isoformLevel.data,geneLevel.data=NULL){
  nq.vec=unlist(nq.list)
  if (is.null(side.nqs)) side.nqs=detectSide(dmix.list,nq.vec)
  pattern.type=lapply(c(1:length(side.nqs)), function(x){
    # 1: 1 components
    if (nq.vec[x]==1){return ("I")}
    if (nq.vec[x]==2){
      #2, 5 & 10: 2 components
      mydmix=dmix.list[x][[1]][[nq.vec[x]]];		
      if (sum(abs(mydmix$delta)<0.75) == 0) return ("V") else{ 
        if (side.nqs[x]==2) return("X") else return("II");
      }
    }else{
      # 4, 10, 11: 3 components
      if (side.nqs[x]==1) return ("VI");
      s.names=unlist(strsplit(names(nq.vec[x]),"\\^"))
      if (length(s.names)==3){
        g.name=s.names[1]	
        iso1.name=s.names[2]
        iso2.name=s.names[3]
      } else{
        cat("\n Warning: weird name ",x," should recheck ... ")
        g.name=paste(s.names[1],"^",s.names[2],sep="")
        iso1.name=s.names[3]
        iso2.name=s.names[4]
      }
      if (g.name != iso2.name){
        x0= isoformLevel.data[which(rownames(isoformLevel.data)==iso1.name),]
        x1= isoformLevel.data[which(rownames(isoformLevel.data)==iso2.name),]
      }else{
        x0= isoformLevel.data[which(rownames(isoformLevel.data)==iso1.name),]
        x1= geneLevel.data[which(rownames(geneLevel.data)==iso2.name),]
      }
      deltaVal=x0-x1;
      if (sum(x0>0 & x1>0)/length(x0) <=0.01) return("X") else return("XI")			
    }
  })
  pattern.type=unlist(pattern.type)

  merged.pattern.type=pattern.type
  merged.pattern.type[which(pattern.type=="II")]="I"
  merged.pattern.type[which(pattern.type=="VI")]="V"
  merged.pattern.type[which(pattern.type=="XI")]="X"

  return(list(pattern.type=pattern.type,merged.pattern.type=merged.pattern.type))
}