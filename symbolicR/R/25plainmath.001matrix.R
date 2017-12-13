# not symbolic math, but funny function for operating plain (atomic) matrix

#' direct sum of linear mapping
#' @param a matrix
#' @param b matrix
#' @return diag block matrix a,b
matrix.direct.sum2 = function(a=as.matrix(a),b=as.matrix(b)){
    re = matrix(0, NROW(a)+NROW(b), NCOL(a)+NCOL(b))
    re[1:NROW(a),1:NCOL(a)] = a
    re[NROW(a)+1:NROW(b), NCOL(a)+1:NCOL(b)] = b
    re
}

#' direct sum of linear mapping
#' @param l list of l
#' @return diag block matrix
matrix.direct.sum.list = function(l){
    l = lapply(l, as.matrix)
    dims = sapply(l, dim)
    re = array(0, rowSums(dims))
    last = c(0,0)
    for(i in 1:NCOL(dims)){
        re[ last[1]+(1:dims[1,i]), last[2]+1:dims[2,i] ] = l[[i]]
        last = last + dims[,i]
    }
    re
}

#' handy interface
#' useful for combine block diagonal matrix to a full matrix
#' @param ... any matrix
#' @return matrix
matrix.direct.sum.n = function(...){
   matrix.direct.sum.list(list(...))
}

#' inverse changes
#' change a matrix into block diagonal format, if possible
#' Not necessary to be symmetric, but need to be square
#'
#' @param m0 square matrix
#' @param eps threshold of zero
#' @param return.split.pt choose if return the list or matrix or just split points
#' @return a list of diagonal blocks or split
matrix.decompose.to.list = function(m0, eps = .Machine$double.eps * 1e4, return.split.pt = FALSE){
    stopifnot(NROW(m0)==NCOL(m0))
    m = abs(m0) >= eps
    m[row(m)==col(m)] = FALSE
    N = NROW(m)
    # divide it up
    split.pt = list()
    walker = function(x){
        k = x
        while(k<N){
            if ( any(m[(k+1):N, x:k ]) ||
                 any(m[ x:k, (k+1):N]) ) {
                # off diagonal non-zero
                k = k + 1
            } else {
                # x:k is a block
                split.pt[[length(split.pt)+1]] <<- c(x,k) 
                Recall(k+1)
                return(NULL)
            }
        }
        split.pt[[length(split.pt)+1]] <<- c(x,k) 
    }
    walker(1)
    if (return.split.pt) {
        return(split.pt)
    }
    lapply(split.pt, function(x) m0[x[1]:x[2],x[1]:x[2]])
}


#' find a permutation for better organization for blockrize
#' Here we assume \code{m} is symmetric
#' @param m matrix
#' @return permutation
matrix.blockrize.permutation = function(m, eps = .Machine$double.eps * 1e4) {
    N = NROW(m)
    if (N==1) return(1)
    G = list(vertex=1:N)
    edges = list()
    for(i in 1:(N-1)) {
        for(j in (i+1):N) {
            if (abs(m[i,j]) >= eps) {
                edges[[ length(edges)+1 ]] = c(i,j)
            }
        }
    }
    G$edges = edges
    L = divide.graph(G)
    unlist(L)
}

#' Give more concise description of matrix list
#' assume each entry in ml are square
#' @param ml matrix list
#' @param prefix prefix for constructing names
#' @return list of diagonal and block
collect.matrix.list.to.diagonal.and.block = function(ml,prefix='ETA'){
    N = length(ml)
    L = list()
    total.etas = 0
    #
    vec = integer(0)
    dim.names = character(0)
    #
    for( i in 1:N ) {
        if ((block.size<-NROW(ml[[i]]))>1) {
        # a REAL block , need to do two things:
        # first,  record previous possiblely merged vec and dim.names
        # second, store the current block
            if (length(vec)>0) {
                L[[ length(L)+1 ]] = list(type='diagonal',val=vec, dim.names=dim.names)
                #
                vec = integer(0)
                dim.names = character(0)
                #
            }
            tmp.names = attr(ml[[i]],'dim.names')
            if (is.null(tmp.names)) {
                tmp.names = sprintf('%s[%s]',prefix, (total.etas+1):(total.etas+block.size))
            }
            total.etas = total.etas + block.size
            L[[ length(L)+1 ]] = list(type='block', val=ml[[i]], dim.names=tmp.names)
            next()
        }
        # block with only one entry, merge them
        total.etas = total.etas + 1
        vec = c(vec, ml[[i]])
        if (is.null(attr(ml[[i]], 'dim.names'))) {
            dim.names = c(dim.names, sprintf('%s[%s]', prefix, total.etas))
        } else {
            dim.names = c(dim.names, attr(ml[[i]], 'dim.names'))
        }
    }
    if (length(vec)>0) {
        L[[ length(L)+1 ]] = list(type='diagonal',val=vec, dim.names = dim.names)
    }
    L
}

#' assign matrix names
#' @param ml matrix list
#' @param dim.names overall names
#' @return matrix list
inject.names.into.list.of.matrix = function(ml,dim.names) {
    N = length(ml)
    if (N<1) return(ml)
    ind = 0
    for(i in 1:N){
        ml[[i]] = as.matrix(ml[[i]])
        sz = NROW(ml[[i]])
        nms = dim.names[ (ind + 1):(ind + sz) ]
        attr(ml[[i]], 'dim.names') = nms
        colnames(ml[[i]])=nms
        rownames(ml[[i]])=nms
        ind = ind + sz
    }
    ml
}

#' get dim names from a matrix list
#' @param ml matrix list
#' @param prefix optional prefix
#' @return names
overall.dim.names.matrix.list = function(ml, prefix='ETA') {
    txt = character(0)
    ind = 0
    for(i in seq_along(ml)){
        sz = NROW(ml[[i]])
        if (is.null(nms<-attr(ml[[i]],'dim.names'))){
            nms = sprintf('%s[%s]', prefix, (ind+1):(ind+sz))
        }
        ind = ind + sz
        txt = c(txt, nms)
    }
    txt
}
