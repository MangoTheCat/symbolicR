#' Of type, we deal with one variable only
#' \code{ifelse(var1==a1, .. , ifelse( var1==a2, .. ,..))}
#' @param e expression
#' @param vt factor type variable
#' @param lvs also can give directly the levels
#' @return design matrix of \code{vt} 
design.matrix.of.one.factor = function(e, vn, vt=NULL, lvs = levels(vt)) {
    if (is.integer(lvs)) {
    # we usually only use numeric in our symbolic expression
        lvs = as.numeric(lvs)
    }
    re = list()
    for(v0 in lvs) {
        v1 = simplify(eval.symbolic( list( CONS('=', vn, v0 ), e ) ))
        re[[ as.character(v0) ]] = v1
    }
    re
}

#' match a simpliest pattern
#' \code{ifelse(var1,..., ifelse(var1,...,ifelse))}
pattern.factorizable.one.ifelse = function(e0) {
    if (symbolic.match(quote(ifelse( '?s'(var.name) == '?a'(level.name), '?'(arg1),  '?'(arg2) )), e0) ) {
        var0 = var.name
        l = list()
        l[[ as.character(level.name) ]] = arg1
        e = arg2
        while(!(is.atomic(e)&&is.na(e))){
            if (symbolic.match(quote(ifelse( '?s'(var.name) == '?a'(level.name), '?'(arg1),  '?'(arg2) )), e) ) {
                if (var0 != var.name) {
                    warning('interaction found')
                    return(NULL)
                }
                if (!is.null(l[[ as.character(level.name) ]])) {
                    warning('previous appearing level')
                }else{
                    l[[ as.character(level.name)]] = arg1
                    e = arg2
                }
            } else {
                break()
            }
        }
        if (!(is.atomic(e) && is.na(e))) {
#            warning('The final else part is not NA, however, 
#            dropped in the output factor, the result may be nonsense')
        # we do not record the ifelse(A==1,1,0) type ifelse because they do not correspond to factor pattern
            return(NULL)
        }
        re = list(var.name=var.name, levels=l)
        class(re) = 'factor.covariate'
        return(re)
    }else{
        return(NULL)
    }
}

#' pattern.affine.mapping.factor.effect
#'
#' match a the following pattern 
#' \code{ifelse(var1,..., ifelse(var1,...,ifelse))}
#' and which can be converted into affine mapping
#' @param e0 expression
#' @return a \code{affine.mapping.factor.effect} type object which contain both the levels and covariate name and the morphism(affine mapping) of expression corresponding to levels
pattern.affine.mapping.factor.effect = function(e0){
    pat = pattern.factorizable.one.ifelse(e0)
    if (is.null(pat)) return(NULL)
    affine = pattern.affine.mapping.of.theta(pat$levels)
    if (is.null(affine)) return(NULL)
    lvs = names(pat$levels)
    re = list(var.name=pat$var.name, level.name = lvs, affine.mapping=affine)
    class(re)='affine.mapping.factor.effect'
    re
}

#' we do not want to make things too complex, hence no nested will be detected
#' Only those like 
#' \code{ ifelse(,,ifelse(,,ifelse(,,NA))) * something * THETA[i] } will de decomposed
#' @param e0 expression
#' @return if matched, type \code{Product.of.Morphisms} morphism
pattern.product.expression.containing.ifelse = function(e0){
    le = expand.as.product(e0)
    if (length(le)==1) return(NULL)
    mor = vector(length(le),mode='list')
    has.ifelse=FALSE
    for(i in seq_along(le)){
        tmp = explain.expression.as.morphism(le[[i]])
        mor[[i]]=tmp
        if(!has.ifelse && tmp$defn$type=='affine.mapping.factor.effect') {
            has.ifelse=TRUE
        }
    }
    if (!has.ifelse) return(NULL)
    morphism.concrete.new('Product.of.Morphisms', mor)
}

#' pattern.affine.function.theta.or.ifelse 
#'
#' similar as above \code{\link{pattern.product.expression.containing.ifelse}}
#' should be a morphism level of \code{\link{expand.as.sum.of.product}} .
#' For example, \code{ S1 + S2 + S3 } will be decomposed, 
#' if those expression \code{Si} is matched by 
#' \code{\link{pattern.product.expression.containing.ifelse}}.
#'
#' @param e0 expression
#' @return if matched, a \code{Sum.of.Morphisms} morphism
pattern.affine.function.theta.or.ifelse = function(e0) {
    le = expand.as.sum(e0)
    if (length(le)==1) return(NULL)
    mors = list()
    can.simplify = FALSE
    for(i in seq_along(le)){
        tmp = explain.expression.as.morphism(le[[i]])
        mors[[i]] = tmp
        if (!can.simplify && 
            tmp$defn$type %in% 
                c('affine.mapping.factor.effect',
                  'Product.of.Morphisms',
                  'Sum.of.Morphisms')){
             can.simplify=TRUE
        }
    }
    if (!can.simplify) return(NULL)
    morphism.concrete.new('Sum.of.Morphisms', mors)
}

#' match \code{Sum c_i * THETA[i]}
#' @param e expression
#' @return row matrix and intercept
pattern.affine.function.of.theta = function(e){
    le = create.polynomial.of.theta(e) 
    A = list()
    b = 0
    parameter.space=integer(0)
    for(i in seq_along(le)){
        if( is.atomic(le[[i]][[2]]) && le[[i]][[2]] == 1) {
            b = le[[i]][[1]]
            next()
        }
        if (symbolic.match(quote(THETA['?a'(t.ind)]), le[[i]][[2]])) {
            A[[length(A)+1]] = le[[i]][[1]]
            parameter.space = c(parameter.space, t.ind)
            next()
        }
        #non linear part found
        return(NULL)
    }
    re =list(A=A,b=b,parameter.space=parameter.space)
    class(re)='affine.function'
    return(re)
}

#' multiple expressions
#' @param mes multiple expressions as a vector mapping
#' @return matrix and right vector if the whole mapping is an affine mapping
pattern.affine.mapping.of.theta = function(mes) {
    if (length(mes)==0) return(NULL)
    if (length(mes)==1) return( pattern.affine.function.of.theta(mes[[1]]))
    # get the whole parameter.space
    l0 = lapply(unname(mes), pattern.affine.function.of.theta)
    if (any(sapply(l0, function(x) is.null(x)))) return(NULL)
    parameter.space = sort(unique(unlist(sapply(l0, function(x) x$parameter.space))))
    N = length(parameter.space)
    # Merge row matrix into a large matrix
    L = vector(length(l0),mode='list')
    for(i in 1:length(l0)){
        v = vector(N,mode='list')
        for(j in 1:N){
            ind = which(parameter.space[j] == l0[[i]]$parameter.space)
            if(length(ind)==0) {
                v[[j]] = 0
            } else {
                v[[j]] = l0[[i]]$A[[ind]]
            }
        }
        L[[i]]=v
    }
    re = list(A=symbolic.listoflist.to.RMatrix(L),b=lapply(l0, function(x) x$b), parameter.space=parameter.space)
    class(re)='affine.mapping'
    re
}

#' convert atomic matrix to symbolic affine mapping
#' @param A matrix
#' @param b optional vector
#' @return affine.mapping
convert.matrix.to.affine.mapping = function(A,b=NULL) {
    if (is.matrix(A)) {
        M = NROW(A)
        N = NCOL(A)
        if (is.null(b)) b = matrix(0, nrow=M, ncol=1)
        re = list(A = apply(A, 1, as.list), b = b)
    } else if (is(A,'CMatrix')) {
        if (is.null(b)) b = matrix(0, nrow=dim.matrix(A)[1], ncol=1)
        re = list(A = asRMatrix.CMatrix(A) , b = b)
    } else if (is(A,'RMatrix')) {
        if (is.null(b)) b = matrix(0, nrow=dim.matrix(A)[1], ncol=1)
        re = list(A = A, b = b)
    }
    class(re)='affine.mapping'
    re
}

#' account the THETA's or THETA_HAT's in a list of expressions
number.of.free.parameters = function(le){
    gg = create.abstractGraph.from.equations(le)
    tol = 0
    for(i in gg$rhs.ind){
        if (symbolic.match(quote(THETA['?a'(t.ind)]), gg$vertex.table[[ i]]) ||
            symbolic.match(quote(THETA_HAT['?a'(t.ind)]), gg$vertex.table[[ i]])) {
            tol = tol+1
        }
    }
    tol
}

#' pretty printer for the affine mapping
#' the working horse
print.affine.mapping0 = function(mapping){
    A = matrix(sapply(mapping$A, function(x) sapply(x, function(y) paste(deparse(y),collapse=''))),
        nrow=length(mapping$A),byrow=TRUE)
    b = matrix(sapply(mapping$b, function(x) paste(deparse(x),collapse='')),nrow=length(mapping$b),byrow=T)
    re = as.table(cbind('[',A,']','','','','','','[',b,']'))
    rownames(re)=rep('',NROW(re))
    colnames(re)=c('',sprintf('A%s',seq_along(mapping$A[[1]])),rep('',7),'b','')
    re = paste(capture.output(print(re)),collapse='\n')
    re
}

#' pretty printer for the affine mapping
#' the working horse
print.affine.mapping.factor.effect0 = function(factors.affine){
    mapping = factors.affine$affine.mapping
    A = matrix(sapply(mapping$A, function(x) sapply(x, function(y) paste(deparse(y),collapse=''))),
        nrow=length(mapping$A),byrow=TRUE)
    b = matrix(sapply(mapping$b, function(x) paste(deparse(x),collapse='')),nrow=length(mapping$b),byrow=T)
    re = as.table(cbind('[',A,']','','','','','','[',b,']'))
    rownames(re)=sprintf('%s.%s',paste(deparse(factors.affine$var.name),collapse=''),factors.affine$level.name)
    colnames(re)=c('',sprintf('A%s',seq_along(mapping$A[[1]])),rep('',7),'b','')
    re = paste(capture.output(print(re)),collapse='\n')
    re
}

#' the S3 version
print.affine.mapping = function(mapping) {
    cat(print.affine.mapping0(mapping),'\n')
}

#' the S3 version
print.affine.mapping.factor.effect = function(factors.affine){
    cat(print.affine.mapping.factor.effect0(factors.affine),'\n')
}
