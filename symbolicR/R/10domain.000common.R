#' these types doesn't have pk meaning\cr
#' even tempvar can be some PK parameter, e.g.\cr
#' if we added a line\cr
#'
#' \code{S2 = V}
#'
#' then, \code{V} becomes a tempvar, not necessary a terminal
#'
#' @param G graph of equations
#' @param v variable name
#' @return type of variable in the graph
#' 
guess.variable.type.by.graph=function(G,v='V'){
    if (v %in% c('THETA[]','ETA[]','EPS[]')) return(v)
    if (!(v %in% G$vertex)) {
        return('undef')
    }
    v0 = neighbourhood.of.v(G,v)
    if(length(v0$v.income)==0 && length(v0$v.outcome)==0) {
        return('unused')
    }
    if(length(v0$v.income)==0){
        return('opening')
    }
    if(length(v0$v.outcome)==0){
        return('terminal')
    }
    'tempvar'
}

#' test any entry of the list is a nonmem parameters
#' if any element of \code{le} is a nonmem parameter \code{THETA}, \code{ETA} or \code{EPS},
#' return \code{TRUE}, else \code{FALSE}
#' @param le  list of expressions
#' @return logical value
list.any.is.nonmem.par = function(le){
    for(e in le){
        if (is.call(e) && length(e)>1 && e[[1]]=='[' 
            && (e[[2]]=='THETA' || e[[2]]=='ETA' || e[[2]]=='EPS')) return(TRUE)
    }
    FALSE
}

#' test if the expression has some \code{ETA[i]} or \code{EPS[i]}
#' @param e expression
#' @return logical value 
expressions.has.random.variables = function(e) {
    foo = extract.vertices(e)
    for(j in foo$rhs.symbol){
        if (is.call(j) && length(j)>2 && j[[1]]=='['
            && (j[[2]]=='ETA' || j[[2]]=='EPS')) return(TRUE)
    }
    FALSE
}

#' guess.variable.meaning.by.name.default
#' 
#' for a list or character, output the corresponding possible type of the variable
#'
#' @param v variable name
#' @return type
guess.variable.meaning.by.name.default = function(v){
    if ( (is.list(v) || is.character(v)) && length(v)>1) {
        L = character(length(v))
        for(i in 1:length(v)){
            L[i]=Recall(v[[i]])
        }
        names(L)=v
        return(L)
    }
    if (is.symbol(v)) v = as.character(v)
    if (is.call(v)) {
        if (v[[1]]=='[') {
            if (is.symbol(v[[2]])) {
                return(Recall(sprintf('%s[]', as.character(v[[2]]))))
            } 
        }
        return('unknown')
    }
    if (v %in% c('CL','KA','K','KM'))  return('pk.basic')
    if (length(grep('^K\\d+$',v))>0) return('pk.basic.elimination.rate')
    if (v %in% c('VM','VC','V')) return('pk.by.trans')
    if (length(grep('^V\\d+$',v))>0) return('pk.basic.volumn')
    if (v %in% c('F1','F2','F3')) return('pk.output.fraction')
    if (v %in% c('S1','S2','S3','SC','S0')) return('pk.scaling')
    if (v %in% c('ALAG1','ALAG2')) return('pk.absorb.lags')
    if (v %in% c('THETA[]','ETA[]','EPS[]')) return('nonmem.parameter')
    if (grepl('^MU_\\d+$', v)) return('pk.mu.referencing')
    if (v %in% c('SEX','RACE','WT','AGE')) return('covariant')
    'unknown'
}

#' default filter for set of critical points
#'
#' People may write their own filter on purpose. However, this filtered is only used for picking up those \cr
#' known PK parameter (by \code{\link{guess.variable.meaning.by.name.default}}) which are also left hand side variables\cr
#' 
#' @param gg graph object created by \code{\link{create.abstractGraph.from.equations}}
#' @param variables variables to be filtered
#' @param eqns definition equations
#' @return the indices of socp for \code{variables}
socp.filter.default = function(gg, variables=NULL, eqns=NULL){
    if (is.null(variables)) variables = gg$vertex.table
    var.types = guess.variable.meaning.by.name.default(variables)
    intersect(grep('^pk\\.', var.types), gg$lhs.ind)
}

#' This filter uses the following idea: 
#' if a lhs variable is defined by some appearant pattern like lognormal or normal
#' It is probably being a critial variable
#' @param gg graph
#' @param variables variables list
#' @param eqns definition equations
#' @return the indices of variables
socp.filter.known.pattern = function(gg, variables=NULL, eqns=NULL){
    if (is.null(variables)) variables = gg$vertex.table
    var.types = guess.variable.meaning.by.name.default(variables)
    total.inds = intersect(grep('^pk\\.', var.types), gg$lhs.ind)
    for(i in 1:length(gg$lhs.ind)) {
        if (gg$lhs.ind[[i]] %in% total.inds) next()
        if (!is.null(pattern.distrib.normal(eqns[[i]][[3]])) ||
            !is.null(pattern.distrib.lognormal(eqns[[i]][[3]])) ) {
            total.inds = union(total.inds, gg$lhs.ind[[i]])
        }
    }
    total.inds
}
