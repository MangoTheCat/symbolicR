#' analyse one connected component
#' to check if there is any following patterns
#'
#' @param cp connected.component
#' @return found parameters used by multiple critical parameters
analyse.common.used.parameter.in.one.component = function(cp){
    theta <- eta <- list()
    rules = list(
        list( quote( THETA['?a'(PARA.IND)] ) , function(e,dict) { used.theta <<- union(used.theta,dict$PARA.IND);e } ),
        list( quote( ETA['?a'(PARA.IND)] ) , function(e,dict) { used.eta <<- union(used.eta,dict$PARA.IND);e } ) )
    extractor = symbolic.simplify.gigo(rules)
    L1 <- L2 <- list()
    for(i in 1:length(cp)){
        used.theta <- used.eta <- integer(0)
        x = cp[[i]]
        if (x$type=='directly.related.to.nonmem.parameters') {
            y = x$defn
            extractor(y$mean)
            extractor(y$var)
            theta[[i]]=used.theta
            eta[[i]]  =used.eta
            if (i==1) next()
            for(j in 1:(i-1)){
                if (length(tmp<-intersect(used.theta , theta[[j]]))>0) {
                    L1[[ length(L1)+1 ]] = list(
                        pairing=c(asCharacterSymbolCall(cp[[i]]$critical.varname), 
                                  asCharacterSymbolCall(cp[[j]]$critical.varname)),
                        common.par=tmp)
                }
                if (length(tmp<-intersect(used.eta , eta[[j]]))>0) {
                    L2[[ length(L2)+1 ]] = list(
                        pairing=c(asCharacterSymbolCall(cp[[i]]$critical.varname), 
                                  asCharacterSymbolCall(cp[[j]]$critical.varname)),
                        common.par=tmp)
                }
            }
        }
    }
    list(L1=L1,L2=L2)
}

#' find parameter spaces embedding in the result of \code{analyse.0.0}
#' 
#' A possible solution for this issue, is to add one more layer as the embedding from a lower dimensional
#' space into the full dimensional parameter space
#'
#' @param cps connected components, the \code{analysed.defns} of returned object of \code{analyse.0.0}
#' @return if any duplication(one parameter used by multiple critical variables) is found
analyse.0.1 = function(cps) {
    # 
    duplicated = F
    for(i in cps){
        # check if parameters commons
        commons = analyse.common.used.parameter.in.one.component(i)
        if(length(commons$L1)>0) {
            warning(sprintf('equations has critical variables depend explicitly on the same THETAs\n %s',
            paste( sapply(commons$L1, function(z) sprintf('[ %s ] on : THETA [ %s ]',
                paste(z$pairing, collapse=','),
                paste(z$common.par, collapse=','))), collapse='\n')))
            duplicated = T
        }
        if(length(commons$L2)>0) {
            warning(sprintf('equations has critical variables depend explicitly on the same ETAs\n %s',
            paste( sapply(commons$L2, function(z) sprintf('[ %s ] on : ETA [ %s ]',
                paste(z$pairing, collapse=','),
                paste(z$common.par, collapse=','))), collapse='\n')))
            duplicated = T
        }
    }
    #
    duplicated
}

#' get.parameter.space.embedding.for.cp
#' 
#' The meanning for this function is the following:
#'
#' For a connected component, e.g. 
#'
#' \code{ V = A + ETA[1] }
#' \code{ CL = B + ETA[1] }
#'
#' convert into
#'
#' \code{ V = A + ETA.HAT[1] }
#' \code{ CL = B + ETA.HAT[2] }
#'
#' where an extra embedding from one dimentional space \code{ETA[1]} to the two dimensional space \code{(ETA.HAT[1], ETA.HAT[2])}
#' with the mapping function
#'
#' \code{ ETA[1] --> (ETA[1], ETA[1]) }
#'
#' @param cp one connected component
#' @return a list with modified cp(you may use or not use it) and the theta.hat and eta.hat defined by theta and eta
get.parameter.space.embedding.for.cp = function(cp){
    theta.hat.n = 0
    eta.hat.n = 0
    # theta.hat : for a entry: theta.hat[i] =j 
    # it means theta[j] --> theta.hat[i] , i.e. , is from theta spaces to a higher dimensional theta.hat spaces
    theta.hat = integer(0)
    eta.hat = integer(0)
    rules = list(
        list( quote( THETA['?a'(PARA.IND)] ) , function(e,dict) {
            ind = NULL
            if (dict$PARA.IND %in% current.thetas) {
                ind = which(theta.hat == dict$PARA.IND)
            } else {
                # add one THETA_HAT
                theta.hat.n <<- theta.hat.n + 1
                theta.hat[ theta.hat.n ] <<- dict$PARA.IND
                current.thetas <<- union(current.thetas, dict$PARA.IND)
                ind = theta.hat.n
            }
            CONS( '[', quote(THETA.HAT), as.numeric(ind))
            } ),
        list( quote( ETA['?a'(PARA.IND)] ) , function(e,dict) {
            ind = NULL
            if (dict$PARA.IND %in% current.etas) {
                ind = which(eta.hat == dict$PARA.IND)
            } else {
                # add one THETA_HAT
                eta.hat.n <<- eta.hat.n + 1
                eta.hat[ eta.hat.n ] <<- dict$PARA.IND
                current.etas <<- union(current.etas, dict$PARA.IND)
                ind = eta.hat.n
            }
            CONS( '[', quote(ETA.HAT), as.numeric(ind))
            } ) )
    walker = symbolic.simplify.gigo(rules)
    for(i in 1:length(cp)){
        current.thetas = integer(0)
        current.etas = integer(0)
        o = cp[[i]]
        if(o$type=='directly.related.to.nonmem.parameters'){
               o$defn$mean = walker(o$defn$mean)
               o$defn$var = walker(o$defn$var)
               cp[[i]] = o
        }
    }
    list(cp.hat = cp, theta.hat=theta.hat, eta.hat=eta.hat)
}
