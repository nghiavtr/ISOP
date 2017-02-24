#' Do Gaussian mixture modelling for isoform expression differences
#'
#' @param dstat A vector of isoform difference
#' @param s A vector of initital values of standard deviation of components
#' @param delta A vector of initial values of mean of components
#' @param p A vector of initial values of proportion of components
#' @param tbreak The number of breaks to create bins
#' @param useExt Option of binning tails or not
#' @param delta.thres Used for estimated p0 
#' @return A list of results and parameters of the mixture model
#' @export
#' @examples
#' data(isoformDataSample)
#' #preprocessing
#' isoformDataSample=ifelse(isoformDataSample <= 3,0,isoformDataSample)
#' isoformDataSample=isoformDataSample[which(rowSums(isoformDataSample)>0),]
#' #tranform read count dataset to log scale
#' isoformDataSample=ifelse(isoformDataSample==0,0,log2(isoformDataSample))
#' #now data is ready
#' tbreak=round(sqrt(ncol(isoformDataSample)))
#' x0=isoformDataSample[1,]
#' x1=isoformDataSample[2,]
#' deltaVal=x0-x1
#' s=c(1,1,1);
#' delta=c(-5,0,5) 
#' res=doGaussianMixtureModel(deltaVal,s=s,delta=delta,tbreak=tbreak)

doGaussianMixtureModel <-
function (dstat, s=NULL, delta, p=NULL, tbreak=NULL,useExt=FALSE, delta.thres = 0.75) 
{
    n = length(dstat)
    if (is.null(tbreak)) {
        tbreak = floor(sqrt(n))
    }

    if (length(tbreak) == 1) {
        rr = range(dstat, na.rm = TRUE)
        tbreak = seq(rr[1], rr[2], length = tbreak)
    }
    nbreaks = length(tbreak)
    tbreak = sort(tbreak)
    if (useExt) {
        tbreak[1] = -Inf
        tbreak[nbreaks] = Inf
    }
    y = table(cut(dstat, tbreak))

    ff = function(param, df, nq, ng, y, tbreak) {
        delta = c(param[1:nq])
        beta = param[(nq+1):(2 * nq)]
        sigma = param[(nq*2+1):(3 * nq)]
        sigma=sigma*sigma
        sigma = ifelse(sigma <0.01,0.01,sigma)

        p = exp(beta)/(sum(exp(beta)))        
        FF = matrix(0, nrow = length(tbreak), ncol = nq)
        FF[, 1] = pnorm(tbreak,delta[1],sigma[1])

        if (nq>1)
        for (j in 2:nq) {
            FF[, j] = pnorm(tbreak, delta[j],sigma[j])
        }
        ff.val = apply(FF, 2, diff)        
        fit = ng * c(ff.val %*% p)+1.0e-100
        -sum(log(fit) * y)
    }

    nq=length(delta)
    if (is.null(p)) p=rep(1/nq,nq)    
    param0 = c(delta, log(p),s)
    oo = optim(param0, ff, df = df, nq = nq, ng = n, y = y, tbreak = tbreak)
    p = oo$par[(nq+1):(2 * nq)]
    p = exp(p)/(sum(exp(p)))    
    delta = oo$par[1:nq]
    s=oo$par[(nq*2+1):(3 * nq)]
    AIC = 2 * oo$value + 2 * (length(delta)+length(s))    
    delta0.pos = which.min(abs(delta))
    p0.raw = p[delta0.pos]
    p0.est = sum(p[abs(delta - delta[delta0.pos]) < delta.thres])

    list(p = p, delta = delta, s = s, AIC = AIC,opt = oo,p0.est=p0.est, p0.raw=p0.raw,delta0.pos=delta0.pos,param0=param0,tbreak=tbreak)
}
