#' Compute differential patterns of isoform pairs from X2 test between group labels of cells and assignment of cells to components of the mixture models
#'
#' @param cellCompMat.prob A results from assignCellComp function
#' @param group.label A vector of group label of cells
#' @param usePermutation An option to compute empirical p-values
#' @param per.num The number of permutations when usePermutation=TRUE
#' @param useParallel An option for using parallel (=TRUE)
#' @return A list containing model-based theoretical p-values (t.VAL) and permutation-based empirical p-values (e.PVAL)
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
#' map.res=assignCellComp(iso.pair.names,model.res$dmix.list,
#'	       model.res$nq.list, isoformLevel.data=isoformDataSample)
#' group.label=unlist(lapply(colnames(isoformDataSample),function(x) unlist(strsplit(x,"_"))[1]))
#' res=getDPIP(map.res$cellCompMat.prob,group.label,usePermutation=TRUE,per.num=100,useParallel=FALSE)
#' hist(res$e.PVAL, breaks=10)
getDPIP=function(cellCompMat.prob,group.label,usePermutation=TRUE,per.num=10000,useParallel=TRUE){
	t.PVAL=NULL
	e.PVAL=NULL

	if (useParallel){
        res=foreach(i=1:length(cellCompMat.prob),.combine=c) %dopar% {
			mylabel=apply(cellCompMat.prob[[i]],1, function(x) return(which.max(x)))
			table.res=table(mylabel,group.label)
			X2.test.res=suppressWarnings(chisq.test(table.res))
			t_pval=X2.test.res$p.value
			e_pval=NULL
			if (usePermutation){
				observed.statistic=X2.test.res$statistic
				permuted.statistic=NULL
				for (j in 1:per.num){
					permuted.group.label=sample(group.label)
					table.res=table(mylabel,permuted.group.label)
					permuted.statistic=c(permuted.statistic,suppressWarnings(chisq.test(table.res)$statistic))
				}
				e_pval=sum(observed.statistic <= permuted.statistic)/per.num				
			}
			list(t_pval=t_pval,e_pval=e_pval)
        }
        t.PVAL=unlist(res[which(names(res)=="t_pval")])
        e.PVAL=unlist(res[which(names(res)=="e_pval")])
    }else{
		for (i in 1:length(cellCompMat.prob)){
			mylabel=apply(cellCompMat.prob[[i]],1, function(x) return(which.max(x)))
			table.res=table(mylabel,group.label)
			X2.test.res=suppressWarnings(chisq.test(table.res))
			t_pval=X2.test.res$p.value
			t.PVAL=c(t.PVAL,t_pval)
			e_pval=NULL
			if (usePermutation){
				observed.statistic=X2.test.res$statistic
				permuted.statistic=NULL
				for (j in 1:per.num){
					permuted.group.label=sample(group.label)
					table.res=table(mylabel,permuted.group.label)
					permuted.statistic=c(permuted.statistic,suppressWarnings(chisq.test(table.res)$statistic))
				}
				e_pval=sum(observed.statistic <= permuted.statistic)/per.num
				e.PVAL=c(e.PVAL,e_pval)
			}
		}
	}
	names(t.PVAL)=names(cellCompMat.prob)
	if(!is.null(e.PVAL)) names(e.PVAL)=names(cellCompMat.prob)
	return(list(t.PVAL=t.PVAL,e.PVAL=e.PVAL))
}
