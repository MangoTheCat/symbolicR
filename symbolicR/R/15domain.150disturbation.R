# Not easy to describe variable or pk parameter as following a simple distribution
# for describe a random effect

# random.disturbation
# it can be used as a attribute of PK parameter

# The basic idea is to separate the Combination of random effect to fixed effect and the Distribution of 
# random vector themselves

#' enforce a input symbolic element to be a list
#' @param e expression
#' @return a list if success
as.symbolic.vector = function(e){
    if (is.list(e)) return(e)
    if (is.call(e) && e[[1]] == 'TUPLE') {
        return(tuple.to.list(e))
    }
    if (is.call(e) || is.symbol(e)) {
        return(list(e))
    }
    if (is.atomic(e)) {
        return(as.list(e))
    }
    stop(sprintf('Can not convert %s to vector', paste(deparse(e),collapse='')))
}

#' construct a disturbation
#' This is a handy interface for defining your own distribution class
#'
#' \code{ FIXED[1] * exp( RANDOM[1] ) }
#' The code will try to unapply it and determine the index of the fixed and random part, also the dimension of both
#' part. And distribution should contain the information of random part
#'
#' Of course, we are losing the original information on the combinated parameter - the whole one, e.g.
#' \code{ a * exp( ETA[1] ) } , we are not recognizing it as a whole is log normal now
#'
#' @param e expression
#' @param fixed fixed variables in expression
#' @param random random variables in expression
#' @param distribution distribution for random vectors
#' @return a disturbation object which can be added to some variable
random.disturbation.new = function(e, fixed, random, distribution) {
    # always assume input is a vector
    e = as.symbolic.vector(e)
    # get symbols
    syms = create.abstractGraph.from.equations(e)$vertex.table
    # guess the definition if not supplied
    if (missing(fixed)) {
        fixed.ind = unlist(lapply(syms, function(x) if(symbolic.match(quote(FIXED['?a'(III)]), x))
                    III else integer(0)))
        fixed.ind = sort(fixed.ind)
        fixed = lapply(fixed.ind, function(x) CONS('[', quote(FIXED), as.numeric(x)) )
    } else {
        fixed = as.symbolic.vector(fixed)
        fixed.ind = seq_along(fixed)
    }
    if (missing(random)) {
        random.ind = unlist(lapply(syms, function(x) if(symbolic.match(quote(RANDOM['?a'(III)]), x))
                    III else integer(0)))
        random.ind = sort(random.ind)
        random = lapply(random.ind, function(x) CONS('[', quote(RANDOM), as.numeric(x)) )
    } else {
        random = as.symbolic.vector(random)
        random.ind = seq_along(random)
    }
    if (missing(distribution)){
        distribution = distribution.common.new(dim=length(random.ind))
    }
    mor = morphism.unapply( e, c( fixed, random ) )
    re = list(  FIXED = seq_along(fixed.ind),
                RANDOM = length(fixed.ind)+seq_along(random.ind),
                defn = mor,
                distribution = distribution
    )
    class(re) = 'random.disturbation'
    re
}

#' generate a disturbation from expression containing EPS
#' non fixed
#' EPS[i] as random
#' extra backref
#' @param e expression
#' @return disturbation
random.disturbation.new.eps = function(e) {
    le = as.symbolic.vector(e)
    gg = create.abstractGraph.from.equations(le)
    syms = gg$vertex.table
    eps.ind = which(sapply(syms, function(x) symbolic.match(quote(EPS['?a'(III)]), x) ) )
    if (length(eps.ind)>0) {
        eps.sym = syms[ eps.ind ]
        # eps are EPS[i] 's i
        eps = sort(unique(sapply(eps.sym, function(x) x[[3]])))
        re = random.disturbation.new(le, fixed=NULL,
            random = lapply(eps, function(x) CONS('[', quote(EPS), as.numeric(x))) )
        class(re) = union('random.observation', class(re))
        attr(re,'backref') = eps
    } else {
        re = e
    }
    return(re)
}

#' assume THETA[i] are FIXed and ETA[i] are RANDOM
#' @param e expression
#' @return random disturbation
random.disturbation.new.from.theta.eta = function(e) {
    le = as.symbolic.vector(e)
    gg = create.abstractGraph.from.equations(le)
    syms = gg$vertex.table
    theta.ind = which(sapply(syms, function(x) symbolic.match(quote(THETA['?a'(III)]), x) ) )
    eta.ind =   which(sapply(syms, function(x) symbolic.match(quote(ETA['?a'(III)]), x) ) )
    #
    backref = integer(0)
    f.para = NULL
    if (length(theta.ind)>0) {
        theta.sym = syms[ theta.ind ]
        theta = sort(unique(sapply(theta.sym, function(x) x[[3]])))
        backref = c(backref, theta)
        f.para = lapply(theta, function(x) CONS('[',quote(THETA),as.numeric(x)))
    } 
    r.para = NULL
    if (length(eta.ind)>0) {
        eta.sym = syms[ eta.ind ]
        eta = sort(unique(sapply(eta.sym, function(x) x[[3]])))
        backref = c(backref, eta)
        r.para = lapply(eta, function(x) CONS('[',quote(ETA),as.numeric(x)))
    } 
    re = random.disturbation.new(le, fixed=f.para, random = r.para )
    if (length(backref)>0) {
        attr(re,'backref') = backref
    }
    return(re)
}

#' apply a disturbation on given fixed and random parameters
#' @param d disturbation
#' @param fixed fixed part
#' @return expression
disturbation.apply = function(d, fixed, random) {
    if (missing(fixed)) {
        fixed = lapply(seq_along(d$FIXED), function(x) CONS('[', quote(FIXED), as.numeric(x)) )
    } else {
        fixed = as.symbolic.vector(fixed)
    }
    if (missing(random)) {
        random = lapply(seq_along(d$RANDOM), function(x) CONS('[', quote(RANDOM), as.numeric(x)) )
    } else {
        random = as.symbolic.vector(random)
    }
    m = instantiate.morphism(d$defn, c(fixed, random))
    m
}

#' S3 version printer for disturbation
#' @param d disturbation
print.random.disturbation = function(d) {
    tmp = as.symbolic.vector(disturbation.apply(d))
    tmp = symbolic.CMatrix(tmp)
    re = print.Matrix0(tmp)
    re1 = print.distribution.common0(d$distribution)
    if (dim.matrix(tmp)[1]<=1 && nchar(re)+nchar(re1)<80) {
        s = sprintf('Random Disturbation : %s , by %s', re, re1)
    } else {
        s = sprintf('Random Disturbation by %s\n%s', re1, re)
    }
    fdomain = random.disturbation.extract.fixed.domain(d)
    if (!is.null(fdomain)) {
        s = sprintf('%s\nwhere FIXED is restricted in %s',s, paste(deparse(fdomain),collapse=''))
    }
    cat(s,'\n')
}

#' extract the domain of the fixed parameters
#' @param d disturbation
#' @return if fixed domain found, return domain expression
random.disturbation.extract.fixed.domain = function(d){
    if (length(d$FIXED)>0) {
        mat = matrix(unlist(as.low.upper.list.morphism.domain(d$defn$domain)[ d$FIXED ]), ncol=2, byrow=T)
        colnames(mat)=c('Lower','Upper')
        return(as.morphism.domain.matrix(mat))
    } 
    NULL
}

#' only for nonmem importing purpose
#' we may save hidden attribute backref referring to origianl index of the THETA, ETA's
#' this function extract those values
#' @param rd random.disturbation
#' @return the found references
lookup.random.part.backref.random.disturbation = function(rd){
    if ((N<-N.random.effects(rd)) <=0 ) return(NULL)
    backref = attr(rd,'backref')
    if (is.null(backref)) return(rep(NA,N))
    backref[ N.fixed.effects(rd)+1:N ]
}

#' update the domain part of the fixed parameters
#' @param d disturbation
#' @param new.domain new domain expression
random.disturbation.update.fixed.domain = function(d, new.domain) {
    if (length(d$FIXED)>0) {
        #
        mat = matrix(unlist(as.low.upper.list.morphism.domain(d$defn$domain)), ncol=2, byrow=T)
        colnames(mat)=c('Lower','Upper')
        #
        rdomain = as.morphism.domain.matrix(mat[d$RANDOM, ,drop=FALSE])
        d$defn$domain = simplify.cartesian.product( CONS('*',new.domain,rdomain) )   
    }
    d
}

#' class for recording the distribution common distributions
#' @param family distribution family
#' @param dim dimension
#' @param mean mean vector
#' @param cov. covariance matrix
#' @param ...  any thing you want to put to the detail
#' @return a distribution
distribution.common.new = function(family='normal', 
                                    dim=1, 
                                    mean=numeric(dim),
                                    cov.=diag(rep(1,dim)),
                                    ...) {
    re = list(
          family = family,
          dim = dim,
          detail = list( mean = mean, cov. = cov., ...))
    class(re)=c('distribution.common','distribution')
    re
}

#' printer working horse
#' @param d distribution
print.distribution.common0 = function(d){
    dim = d$dim
    mean = d$detail$mean
    cov. = d$detail$cov.
    family=d$family
    if (dim==1) {
        if (family=='normal' && isTRUE(mean==0) && isTRUE(all(cov.==1))) {
            re = 'Standard Normal random variable'
        } else {
            re = sprintf('%s random variable of %s and variance %s',family, 
                ifelse(isTRUE(mean==0), 'zero mean',sprintf('mean value equal to %s', mean)), 
                cov.)
        }
    } else {
        re = sprintf('%s dimensional %s random vector', dim, family)
        if (isTRUE(all(mean==0))) re = sprintf('%s, zero mean',re)
        if (all( cov.[lower.tri(cov., diag=F)]==0 )) {
            if (family=='normal') {
                if (all(diag(cov.)==1)) {
                    re = sprintf('Independent %s dimensional standard normal random vector',dim)
                } else {
                    re = sprintf('Independent %s', re)
                }
            } else {
                re = sprintf('Un-correlated %s', re)
            }
        }
    }
    re
}

#' S3 version
print.distribution.common = function(d){
    cat(print.distribution.common0(d),'\n')
}

#' update
#' @param d distribution
#' @param new.cov new covariance matrix
#' @return updated distribution
update.distribution.covariance = function(d, new.cov ){
    d$detail$cov. = unname(new.cov)
    d
}


#' convertion
#' @param type distribution type
#' @param Var  if atomic, assuming just a RANDOM type, or variation expression \code{coeff*Var(ETA[i])}, then \code{sqrt(coeff)}
#' @return a random.disturbation object
convert.distrib.description.to.random.disturbation = function(type, Var=1) {
    if (is.atomic(Var)) {
        eps = quote(RANDOM[1])
    } else if (is.call(Var)) {
        coeff = extract.coefficient.of.variance(Var)
        eps = CONS('*', coeff, quote(RANDOM[1]))
    } else {
        stop('Don not know how to explain ', Var , 'as random effect.')
    }
    if (type=='log normal') {
        distr = distribution.common.new('normal', cov.=Var)
        re = random.disturbation.new( CONS('exp',CONS('+', quote(FIXED[1]), eps) ) , distribution=distr)
    } else if (type=='normal') {
        distr = distribution.common.new('normal', cov.=Var)
        re = random.disturbation.new( CONS('+', quote(FIXED[1]), eps), distribution=distr )
    } else if (type=='fixed') {
        re = random.disturbation.new( quote(FIXED[1]))
    } else {
        stop(sprintf('Do not know how to convert type %s to disturbation.', type ))
    }
    re
}


#' return number of random effects in a disturbation
#' @param rd random disturbation
#' @return number of RANDOM parameter
N.random.effects = function(rd){
    length(rd$RANDOM)
}

#' return number of fixed in a disturbation, can be 0
#' @param rd random disturbation
#' @return number of FIXed parameter
N.fixed.effects = function(rd){
    length(rd$FIXED)
}

#' The covariance of the random part of a disturbation
#' @param rd random distribution
#' @return the covariance matrix
Cov.random.disturbation = function(rd) {
    rd$distribution$detail$cov.
}

#' S3 for random observation
print.random.observation = function(m) {
    indent.n = function(ss,n) {
       if (n==0) return(ss)
       indent = paste(rep(' ',n),collapse='')
       ss = gsub('^', indent, ss)
       gsub('\n', sprintf('\n%s',indent), ss)
    }
    txt = indent.n(paste(capture.output(print.random.disturbation(m)),collapse='\n'),4)
    txt = sprintf('Random Observation ~\n%s',txt)
    cat(txt,'\n')
}
