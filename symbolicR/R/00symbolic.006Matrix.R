#' convert TUPLE to list which means a symbolic Vector
#' \code{ TUPLE( x1,x2,x3) => list(x1,x2,x3) }
#' \code{
#' TUPLE( 
#'    TUPLE(x1,x2,x3),
#'    TUPLE(x4,x5,x6)) => list( list(x1,x2,x3), list(x4,x5,x6)) }
#' @param e0 quoted expression started by TUPLE
#' @return list
tuple.to.list = function(e0){
    walker = function(e){
        if (is.call(e) && e[[1]]=='TUPLE') {
            N = length(e) - 1 
            if (N>0) {
                for(i in 1:N){
                    e[[i+1]] = Recall(e[[i+1]])
                }
                return(as.list(e[-1]))
            }
            return(list())
        }
        if (is.call(e)) {
            for(i in 1:length(e)){
                e[[i]] = Recall(e[[i]])
            }
        }
        e
    }
    walker(e0)
}

#' enforce expression as a list
#'
#' 1:10 => list(1,2,3,4,...,10)
#' e => list(e)
#' @param e expression or an atomic vector
#' @return list of the above
vectorize.symbol.and.atomic = function(e) {
    if (is.call(e) || is.symbol(e)) return(list(e))
    if (is.atomic(e)) return(as.list(e))
    stop(sprintf('Check if %s can be though as a variable', paste(deparse(e),collapse='')))
}

#' inner product for symbolic vectors
#'
#' @param le1 list of expression, as a symbolic vector
#' @param le2 list of expression, as a symbolic vector
#' @return one expression which is the inner product
symbolic.inner.product = function(le1=vectorize.symbol.and.atomic(le1),le2=vectorize.symbol.and.atomic(le2)) {
    n1 = length(le1)
    n2 = length(le2)
    N = max(n1,n2)
    l = vector(N, mode='list')
    for(i in 0:(N-1)) {
    # simulate R's replication of vector
        l[[i+1]] = simplify.2(CONS('*',le1[[ i %% n1 + 1]], le2[[i %% n2 + 1]]))
    }
    simplify.2(expression.from.list.of.sum(l))
}

#' use \code{\link{symbolic.inner.product}} we get \code{(A==1)*a1 + (A==2)*a2}
#' however, use this function, we get \code{ifelse(A==1,a1,ifelse(A==2,a2,NA))}
#' @param le1 list of logical expression, as a symbolic vector
#' @param le2 list of expression, as a symbolic vector
#' @return one ifelse expression which is the inner product
symbolic.inner.product.logical = function(le1,le2){
   walker = function(v1,v2) {
        if (length(v1)==1) {
            return( CONS('ifelse', v1[[1]], v2[[1]], NA) )
        }
        CONS('ifelse',v1[[1]], v2[[1]], Recall(v1[-1],v2[-1]))
   }
   walker(le1,le2)
}

#' get column i
#' @param A symbolic matrix (list of list)
#' @param i column index
#' @return list as column
symbolic.matrix.col.i = function(A,i){
    lapply(A, function(x) x[[i]])
}

#' sub matrix defined by i,j
#' @param A matrix
#' @param i row index
#' @param j col index
#' @return sub matrix
symbolic.matrix.sub = function(A,i,j){
    lapply(i, function(x) lapply(j, function(y) A[[x]][[y]]))
}

#' compute \code{Ax + b}
#' @param A symbolic matrix
#' @param x symbolic vector
#' @param b right vector
#' @return result vector
symbolic.affine.mapping = function(A, x, b=NULL) {
    M = length(A)
    re = vector(M, mode='list')
    for(i in 1:M){
            re[[i]] = symbolic.inner.product(A[[i]], x) 
    }
    if (!is.null(b)) {
        re = symbolic.vector.addition(re, b)
    }
    re
}

#' multiply two symbolic RMatrix
#' @param A RMatrix
#' @param B RMatrix
#' @return RMatrix
symbolic.RMatrix.mul.RMatrix = function(A,B){
    dm1 = dim.matrix(A)
    dm2 = dim.matrix(B)
    if (dm1[2]!=dm2[1]) stop('un-compatible dimensions')
    C = symbolic.RMatrix(0, dm1[1], dm2[2])
    BT = symbolic.transpose.RMatrix(B)
    for(i in 1:dm1[1]){
        for(j in 1:dm2[2]){
            C[[i]][[j]] = symbolic.inner.product(A[[i]], BT[[j]])
        }
    }
    C
}

#' multiply two symbolic Matrix
#' @param A Matrix
#' @param B Matrix
#' @return RMatrix
symbolic.Matrix.mul = function(A,B) {
    if (is(A,'CMatrix')) A = asRMatrix.CMatrix(A)
    if (is(B,'CMatrix')) B = asRMatrix.CMatrix(B)
    symbolic.RMatrix.mul.RMatrix(A,B)
}

#' symbolic Vector addition
#' @param le1 list of expression, as a symbolic vector
#' @param le2 list of expression, as a symbolic vector
#' @return one list which is the point-wise addition
symbolic.vector.addition = function(le1=vectorize.symbol.and.atomic(le1),le2=vectorize.symbol.and.atomic(le2)){
    n1 = length(le1)
    n2 = length(le2)
    N = max(n1,n2)
    l = vector(N, mode='list')
    for(i in 0:(N-1)) {
    # simulate R's replication of vector
        l[[i+1]] = simplify.2(CONS('+',le1[[ i %% n1 + 1]], le2[[i %% n2 + 1]]))
    }
    l
}

#' symbolic list of list matrix
#' @param le list of expressions
#' @param nrow row number
#' @param ncol column number
#' @param byrow by row or by column
#' @return  list of list representation
symbolic.RMatrix = function( le , nrow=length(le), ncol=1, byrow=TRUE) {
    if (missing(nrow) && !missing(ncol)) {
        nrow = max(round(length(le) / ncol),1)
    } else if (!missing(nrow) && missing(ncol)) {
        ncol = max(round(length(le) / nrow),1)
    }
    N = length(le)
    A = vector(nrow,mode='list')
    for(i in 1:nrow){
        A[[i]] = vector(ncol, mode='list')
        for(j in 1:ncol){
            if (byrow) {
                ind = (i-1)*ncol + j - 1
            } else {
                ind = i+(j-1)*nrow - 1
            }
            A[[i]][[j]]=le[[ ind %% N + 1 ]]
        }
    }
    class(A)=c('RMatrix','Matrix','listoflist')
    A
}

#' convert between listoflist and symbolic.RMatrix
#' @param A a list of list
#' @return RMatrix
symbolic.listoflist.to.RMatrix = function(A){
    stopifnot(is.list(A) && is.list(A[[1]]))
    class(A)=c('RMatrix','Matrix','listoflist')
    A
}

#' similar to R matrix, columns first
#' @param le list of expressions
#' @param nrow row number
#' @param ncol column number
#' @param byrow by row or by column
#' @return  list of list representation
symbolic.CMatrix = function(le , nrow=length(le), ncol=1, byrow=FALSE) {
    if (missing(nrow) && !missing(ncol)) {
        nrow = max(round(length(le) / ncol),1)
    } else if (!missing(nrow) && missing(ncol)) {
        ncol = max(round(length(le) / nrow),1)
    }
    if (byrow) {
        le = unlist(symbolic.transpose.RMatrix(symbolic.RMatrix(le, nrow, ncol, byrow)))
    }
    class(le)=c('CMatrix','Matrix')
    N = length(le)
    if (N < nrow * ncol) {
        for(i in (length(le)+1):(nrow*ncol)){
            le[[i]] = le[[ (i-1) %% N + 1 ]]
        }
    } else if (N > nrow*ncol) {
        le = le[1:(nrow*ncol)]
    }
    attr(le, '.Dim')=c(nrow,ncol)
    le
}

#' convertion, listoflist to matrix
#' @param A RMatrix
#' @return CMatrix
asCMatrix.RMatrix = function(A) {
    symbolic.CMatrix(unlist(A), length(A), length(A[[1]]), byrow = TRUE)
}

#' convert list to listoflist matrix
#' @param A CMatrix
#' @return RMatrix
asRMatrix.CMatrix = function(A) {
    dm = dim.matrix(A)
    symbolic.RMatrix(unlist(A), dm[1],dm[2], byrow = TRUE)
}

#' transpose of RMatrix
#' @param A RMatrix
#' @return RMatrix
symbolic.transpose.RMatrix = function(A) {
    symbolic.RMatrix( unlist(A), length(A[[1]]), length(A), byrow=FALSE)
}

#' return dim of symbolic Matrix
#' @param A Matrix
#' @return dimension vector
dim.matrix = function(A){
    if (is(A,'RMatrix')) {
        return(c(length(A),length(A[[1]])))
    }
    if (is(A,'CMatrix')) {
        return(attr(A,'.Dim'))
    }
    NULL
}

#' formatter for A
print.CMatrix = function(A,...) {
    print(asRMatrix.CMatrix(A))
}

#' construct the string formaterring A
print.Matrix0 = function(A) {
    if (!(is(A,'RMatrix') || is(A,'CMatrix'))) return(NULL)
    if (is(A,'RMatrix')) byrow=TRUE
    else byrow=FALSE
    dm = dim.matrix(A)
    a = matrix(sapply(unlist(A), function(x) paste(deparse(x),collapse='')) , dm[1], dm[2], byrow=byrow)
    a = as.table(cbind('[',a,']'))
    rownames(a) = rep('',dm[1])
    colnames(a) = rep('',dm[2]+2)
    b = capture.output(print(a))
    paste(b[-grep('^\\s*$',b)], collapse='\n')
}

#' S3 version
print.Matrix = function(A) {
    cat(print.Matrix0(A),'\n')
}

#' transpose CMatrix
#' @param A CMatrix
#' @return CMatrix
symbolic.transpose.CMatrix = function(A) {
    dm = dim.matrix(A) 
    symbolic.CMatrix(A, dm[2] , dm[1] , byrow=TRUE)
}
