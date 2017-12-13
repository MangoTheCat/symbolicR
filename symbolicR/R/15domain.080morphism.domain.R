#' create one dimensional interval description
#' @param a left end
#' @param b right end
#' @return domain description expression
morphism.domain.one.dimensional = function(a,b) {
    a = unname(a)
    b = unname(b)
    if (a==-Inf && b==Inf) return( quote( RealNumber^1 ) )
    # empty set
    if (a > b) return( quote( EmptySet^0 ) )
    CONS('^', CONS('[', quote(ClosedInterval), a, b ), 1)
}

simplify.cartesian.product.rules = list(
# combining rule
list(quote( '?'(baseDomain)^'?'(dim1) * '?'(baseDomain)^'?'(dim2) ), 
        CONS('^', quote(':'(baseDomain)), 
            CONS('+', quote(':'(dim1)), quote(':'(dim2))))),
# we can not use simplify because it removes ^1
# need another rule for atomic evaluation
list(quote( '?'(baseDomain)^'?'(dim.expression) ), function(e,dict) {
        new.dim = simplify.2(dict$dim.expression)
        CONS('^', dict$baseDomain, new.dim)
}),
# associating rule: (a*b)*b -> a*(b*b)
list(quote( '?'(X1) * '?'(X2) * '?'(X3)), CONS('*', quote(':'(X1)), quote(':'(X2) * ':'(X3))))
)

#' e.g. RealNumber^1 * RealNumber^1 => RealNumber^2
simplify.cartesian.product = symbolic.simplify.gigo( simplify.cartesian.product.rules)

calculate.ambience.domain.rules = list(
list( quote('?'(baseDomain1)^'?'(dim1) * '?'(baseDomain2)^'?'(dim2)),
    function(e,dict) {
        CONS('^',quote(RealNumber), simplify(CONS('+',dict$dim1,dict$dim2)))}))
#' A^1 * B^2 * C^3 => ambienceDomain^6
#' Which means, the domain A^1 * B^2 * C^3 has to be embeded into a six dimensional Euclidean space
#' Of course, we assume baseDomains should be real, not complex numbers
calculate.ambience.domain = symbolic.simplify.gigo(calculate.ambience.domain.rules)

#' setting functions for the domain of a morphism
#' assume \code{m} has 3 columns low, est, upper corresponding to NONMEM statements
#' return a domain description
#' @param m matrix
#' @param parameter.space the row number to use, or indices for the THETAs
#' @return expression for domain
as.morphism.domain.matrix = function(m, parameter.space=1:NROW(m)) {
    if (length(parameter.space)==0) return(quote( EmptySet^0 ) )
    if (!is.null(colnames(m))) {
        if (is.na(lower<-match('Lower',colnames(m)))) lower = 1
        if (is.na(upper<-match('Upper',colnames(m)))) upper = 3
    } else {
        lower = 1
        upper = 3
    }
    for(i in seq_along(parameter.space)) {
        a = m[parameter.space[i], lower]
        b = m[parameter.space[i], upper]   
        if (i==1) {
            e = morphism.domain.one.dimensional(a,b)
            next()
        }
        e = simplify.cartesian.product(CONS('*', e, morphism.domain.one.dimensional(a,b)))
    }
    e
}

#' return dimension (number) of a domain expression
#' @param domain the domain expression
#' @return the dimension of ambience space
dimension.of.domain = function(domain){
    calculate.ambience.domain(domain)[[3]]
}

#' check if a vector belongs to domain
#' @param m morphism
#' @param vec vector
#' @return TRUE or FALSE
does.morphism.domain.contain = function(m, vec){
    lmat = as.low.upper.list.morphism.domain(m)
    if (length(vec) != length(lmat)) return(FALSE)
    for(i in seq_along(lmat)){
        if (vec[i] <= lmat[[i]][1] || vec[i] >= lmat[[i]][2]) {
            return(FALSE)
        }
    }
    TRUE
}

#' convert morphism
#' to a list of upper and lower ranges of the parameters
#' @param m morphism
#' @return list
as.low.upper.list.morphism.domain = function(m){
    if (is(m,'morphism')) {
        m = m$domain
    } 
    le = expand.as.product(m)  
    m0 = list()
    push = function(v) m0[[ length(m0)+1 ]] <<- v
    eps = .Machine$double.eps
    for(e in le){
        if (symbolic.match(quote( RealNumber ), e)) {
            push( c(-Inf,Inf))
            next()
        }
        if (symbolic.match(quote( PositiveRealNumber ), e)) {
            push( c(eps,Inf))
            next()
        }
        if (symbolic.match(quote( RealNumber^'?'(dim1) ), e)) {
            for(i in 1:dim1){
                push( c(-Inf,Inf))
            }
            next()
        }
        if (symbolic.match(quote( PositiveRealNumber^'?'(dim1) ), e)) {
            for(i in 1:dim1){
                push( c(eps,Inf))
            }
            next()
        }
        if (symbolic.match(quote( ClosedInterval[ '?'(lower), '?'(upper) ] ), e)) {
            push(c(lower,upper))
            next()
        }
        if (symbolic.match(quote( ClosedInterval[ '?'(lower), '?'(upper) ]^'?'(dim1) ), e)) {
            for(i in 1:dim1){
                push(c(lower,upper))
            }
            next()
        }
        stop(sprintf('unknown domain specification %s ', paste(deparse(e),collapse='')))
    }
    m0
}

#' referring to original THETA's
#' We store the \code{THETA[i]}'s index in \code{m$defn$parameter.space}
#' or \code{m$defn$defn$parameter.space} if \code{m} is a \code{groupeffect} type
#' when we need a value or range by lookup from original nmModel
#' this pointer is needed
#' @param morphism
#' @return integer vector for backreferring to original index of THETAs
backref.parameter.of.morphism = function(m){
    if (m$defn$type=='groupeffect') {
        re = m$defn$defn$parameter.space
    } else if (m$defn$type=='morphism.composition') {
        # recursively output the last morphism
        re = Recall(m$defn$defn[[ length(m$defn$defn) ]])
    } else if (m$defn$type=='affine.mapping.factor.effect'){
        re = m$defn$defn$affine.mapping$parameter.space
    } else if (m$defn$type=='Product.of.Morphisms') {
        re = integer(0)
        for(mor in m$defn$defn) {
            re = c(re,Recall(mor))
        }
    } else {
        re  = m$defn$parameter.space
    }
    if(is.null(m)){
        re = rep(NA, degree.of.freedom(m))
    }
    re
}
