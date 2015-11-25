#' Validate the randomness of isoform pair via permutation
#' @param iso.pair.names A vector of names of isoform pairs in the format: gene.isoform1.isoform2
#' @param isoformLevel.data A data matrix of isoform expression
#' @param geneLevel.data A data matrix of gene expression
#' @param per.num A number of permutation
#' @param tbreak The number of breaks to create bins
#' @param esp.val A small value added to denominator of X2 formula
#' @param useExt A boolean value to allow the left-most and right-most bins extending to infinity
#' @param useParallel An option for using parallel (=TRUE)
#' @return A list containing results of X2 tests and p-values
#' @export
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
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
#' #library(doParallel)
#' #registerDoParallel(cores=4)
#' set.seed(2015)
#' iso.pair.names=names(model.res$nq.list)
#' res=validateIsoformPair(iso.pair.names=iso.pair.names,
#'	   isoformLevel.data=isoformDataSample,per.num=100,tbreak=tbreak,useParallel=FALSE)
#' PVAL.X2.list=unlist(res[names(res)=="pval"])
#' hist(PVAL.X2.list)
#' fdr.val=p.adjust(PVAL.X2.list,method="BH")
#' length(which(fdr.val <= 0.05))

validateIsoformPair <-function(iso.pair.names=NULL,isoformLevel.data,geneLevel.data=NULL,per.num=10000,tbreak=10,esp.val=0.5,useExt=FALSE, useParallel=FALSE){	
	getBin=function(dstat,tbreak=NULL,useExt=FALSE,minVal=NULL,maxVal=NULL){
		n = length(dstat)
	    if (missing(tbreak)) {
	        tbreak = floor(sqrt(n))
	    }

	    if (length(tbreak) == 1) {
	    	if (is.null(maxVal)) rr = range(dstat, na.rm = TRUE) else rr=c(minVal,maxVal)	        
	        tbreak = seq(rr[1], rr[2], length = tbreak)
	    }
	    nbreaks = length(tbreak)
	    tbreak = sort(tbreak)
	    if (useExt) {
	        tbreak[1] = -Inf
	        tbreak[nbreaks] = Inf
	    }
	    return (tbreak)
	}

	doValidation=function(myisopairname,isoformLevel.data,geneLevel.data,useExt=FALSE,tbreak=10,esp.val=0.5,per.num=10000){
		s.names=unlist(strsplit(myisopairname,"\\^"))
		if (length(s.names)==3){
			g.name=s.names[1]	
			iso1.name=s.names[2]
			iso2.name=s.names[3]
		} else{
			cat("\n Warning: weird name ",myisopairname," should recheck ... ")
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

		data.iso1=x0
		data.iso2=x1
		deltaVal=data.iso1-data.iso2

		deltaVal.per.mat=lapply(c(1:per.num), function(x) sample(data.iso1)-sample(data.iso2))
		deltaVal.per.mat=do.call("rbind",deltaVal.per.mat)

		minVal=min(deltaVal.per.mat)
		maxVal=max(deltaVal.per.mat)		
		tbreak.vec = getBin(deltaVal,tbreak,useExt=useExt,minVal=minVal,maxVal=maxVal)
		y = table(cut(deltaVal, tbreak.vec))

		bin.per.mat=apply(deltaVal.per.mat,1, function(deltaVal.per) table(cut(deltaVal.per, tbreak.vec)))
		bin.per.mat=t(bin.per.mat)
		mean.y.per=colSums(bin.per.mat)/per.num

		X2.vec=apply(bin.per.mat,1, function(y.per) sum((y.per-mean.y.per)^2/(mean.y.per),na.rm=TRUE))
		mean.X2.val=sum((y-mean.y.per)^2/(mean.y.per),na.rm=TRUE)
		pval=sum(X2.vec>=mean.X2.val)/per.num

		return(list(X2.vec=X2.vec,mean.X2.val=mean.X2.val,pval=pval))
	}

	if(is.null(iso.pair.names)){
		cat("Warning: no isoform-pair names ... QUIT!")
		return;
	}
	
	perRes=NULL;

	if (useParallel){
		perRes=foreach(i=1:length(iso.pair.names),.combine=c) %dopar% {
			myisopairname=iso.pair.names[i]
			res=doValidation(myisopairname=myisopairname,isoformLevel.data=isoformLevel.data,geneLevel.data=geneLevel.data,useExt=useExt,tbreak=tbreak,esp.val=esp.val,per.num=per.num)
		}

	}else{
	
		for (i in 1:length(iso.pair.names)){
			myisopairname=iso.pair.names[i]
			res=doValidation(myisopairname=myisopairname,isoformLevel.data=isoformLevel.data,geneLevel.data=geneLevel.data,useExt=useExt,tbreak=tbreak,esp.val=esp.val,per.num=per.num)

			perRes=c(perRes,res)
		}
	}
	return(perRes)
}
