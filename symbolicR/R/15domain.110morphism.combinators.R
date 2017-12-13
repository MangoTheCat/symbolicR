#' just a check for groupeffect definition, and throw a exception if necessary
#' any morphism, or group definition should be THETA(parameters) free
#' @param gpdefn \code{grpvar.definition}
#' @return TRUE or simplify throw exceptions
check.groupeffect.definition = function(gpdefn=NULL) {
    if (!is.null(gpdefn)) {
        if (!is.list(gpdefn)){
            gpdefn = list(gpdefn)
        }
        for(i in seq_along(gpdefn)){
            if (symbolic.grep(quote(THETA['?a'(ANYI)]), gpdefn[[i]])) {
                stop('Error happened, group definition should not have THETA symbols')
            }
        }
    }
    TRUE
}

#' test if \code{e} has any pattern recognized by \code{analyse.group.effect}
#' @param e expression
#' @return morphism object if \code{e} has a known group effect
morphism.defined.by.groupeffect = function(e){
    pat = analyse.group.effect(e)
    if (!is.null(pat)) {
        if (pat$grptype=='identicalmapping') {
            return(morphism.concrete.new('identicalmapping',dim=pat$dim, parameter.space=pat$parameter.space))
        }
        check.groupeffect.definition(pat$grpvar.definition)
        if (pat$grptype %in% c('additive','multiplicative','affine')) {
            re = list(domain=quote(RealNumber^1),
                      image =quote(RealNumber^1),
                      defn  =list(type='groupeffect',
                                  defn=pat))
        } else if (pat$grptype %in% c('linear','exponential')) {
            re = list(domain=quote(RealNumber^2),
                      image =quote(RealNumber^1),
                      defn  =list(type='groupeffect',
                                  defn=pat))
        } else if (pat$grptype=='binomial.linear') {
            re = list(domain=quote(RealNumber^4),
                      image =quote(RealNumber^1),
                      defn  =list(type='groupeffect',
                                  defn=pat))
        } else {
            stop(sprintf('Unknown group type %s ', pat$grptype))
        }
        class(re) = c('morphism', class(re))
        return(re)
    }
    NULL
}

#' for a morphism defined by groupeffect
#' instantiate it
#' @param m morphism
#' @param parameter.pool
#' @return expression or equations
morphism.groupeffect.instantiator = function(m, parameter.pool){
    if (m$defn$type=='groupeffect') {
        type = m$defn$defn$grptype
        grpvar.definition = m$defn$defn$grpvar.definition
        if (is.null(grpvar.definition)) {
            stop('Morphism is declared as groupeffect, however, grpvar.definition is null')
        }
        re = 
        switch(type,
            additive=CONS('+', parameter.pool[[1]], grpvar.definition ) ,
            multiplicative=CONS('*', parameter.pool[[1]], grpvar.definition ) ,
            affine= CONS('+', grpvar.definition[[1]], 
                        CONS('*', parameter.pool[[1]], grpvar.definition[[2]])),
            linear= CONS('+', parameter.pool[[1]],
                        CONS('*', grpvar.definition, parameter.pool[[2]])),
            exponential=CONS('*', parameter.pool[[1]], 
                    CONS('^', grpvar.definition, parameter.pool[[2]] )),
            binomial.linear={ 
                foo1 = CONS('+', parameter.pool[[1]],
                            CONS('*', grpvar.definition[[2]], parameter.pool[[2]]));
                foo2 = CONS('+', parameter.pool[[3]],
                            CONS('*', grpvar.definition[[2]], parameter.pool[[4]]))
                CONS('+', CONS('*', grpvar.definition[[1]], foo1 ),
                    CONS('*', CONS('-',1, grpvar.definition[[1]]), foo2 )) } )
        return(re)
    }
    NULL
}

#' THETA oriented function
#' This try to guess a meaning of an expression contain THETAs
#' we presume the THETA's as the points in the preimage
#' however, we discard the informations here, i.e. remove THETAs here
#' @param e expression
#' @return a morphism object
explain.expression.as.morphism = function(e) {
    if (is.atomic(e) || is.symbol(e)) return(morphism.concrete.new('Constant',val=e))
    pat = morphism.defined.by.groupeffect(e)
    if (!is.null(pat)) return(pat)
    # we will look at if the exp(e) has any groupeffect pattern
    pat = morphism.defined.by.groupeffect( simplify.collect.log.exp(CONS('exp', e)) )
    if (!is.null(pat) && pat$defn$type!='identicalmapping') {
    # we won't explain trivil identicalmapping
        return(morphism.composition.2(morphism.concrete.new('log.change'), pat))
    }
    # some common patterns used in control stream
    # exp(THETA[i])
    if (symbolic.match(quote(exp(THETA['?a'(ii)])), e) ) {
        return(morphism.concrete.new('exponential.change',parameter.space=ii))
    }
    # log(THETA[i])
    if (symbolic.match(quote(log(THETA['?a'(ii)])), e) ) {
        return(morphism.concrete.new('log.change',parameter.space=ii))
    }
    # logit^-1(THETA[i])
    if (symbolic.match( CONS('/',1, quote(1+exp(-THETA['?a'(ii)]))) , e ) ||
        symbolic.match( CONS('/',quote(exp(THETA['?a'(ii)])), 
                                quote(1+exp(THETA['?a'(ii)]))) , e )) {
        return(morphism.concrete.new('inverse.logit.change',parameter.space=ii))
    }
    # ifelse(,,ifelse(,,ifelse(,,))) pattern
    pat = pattern.affine.mapping.factor.effect(e)
    if (!is.null(pat)) {
        return(morphism.concrete.new('affine.mapping.factor.effect', pat))
    }
    # Only for simplify ifelse(,,ifelse(,,ifelse(,,)))
    pat = pattern.product.expression.containing.ifelse(e)
    if (!is.null(pat)) {
        return(pat)
    }
    # Only for simplify a1 * ifelse(,,ifelse(,,ifelse(,,))) + a2 * ifelse(,,...)
    pat = pattern.affine.function.theta.or.ifelse(e)
    if (!is.null(pat)) {
        return(pat)
    }
    ## General Case ##
    # not previousely defined type, set as general pointwise definition
    # The most important part is to remove THETA's here
    variables = extract.vertices(e)$rhs.symbol
    if (length(variables)==0) {
        thetas=list()
    } else {
        ind = sapply(variables, function(x) symbolic.match(quote(THETA['?a'(ii)]), x ) )
        thetas = variables[ ind ]
    }
    # covariates = variables[!ind]
    if (length(thetas)==0) {
        # too many warnings, better not print it
        # warning(sprintf('Expression %s not a mapping?', paste(deparse(e),collapse='')))
        return(morphism.concrete.new('Constant', val=e))
    }
    # e.g.  THETA[30] - THETA[15] ==> x[1] - x[2] , but parameter.space = c(30,15) for future use
    re = morphism.unapply(e, thetas, parameter.space=sapply(thetas, function(x) x[[3]]))
    class(re) = c('morphism', class(re))
    return(re)
}

#' Currys
#' \code{ f : (x1,x2,x3) |-> f(x1,x2,x3) }
#' \code{f x1} is then \code{ (x2,x3) |-> f(x1,x2,x3) }
#' @param m0 main morphism
#' @param e0 first argument
#' @param ind to replace which args
#' @return morphism
morphism.curry.point = function(m0, e0, ind=1) {
    if (m0$defn$type=='Constant') return(m0)
    # random things into m0
    args = generate.formal.arguments.morphism(m0, fsymbol='RGw2G')
    if (length(ind)>1) {
        e0 = as.symbolic.vector(e0)
        N.e0 = length(e0)
        for(i0 in seq_along(ind)) {
            args[[ ind[i0] ]] = e0[[ (i0-1) %% N.e0 + 1 ]]
        }
    } else {
        # only one argument, allow both list or single item
        if (is.list(e0)) args[[ ind ]] = e0[[1]]
        else args[[ ind ]] = e0
    }
    e = instantiate.morphism(m0, parameter.pool=args)
    morphism.unapply(e, args[-ind])
}

#' if the curryed one is a morphism
#' \code{ f : (x1,x2,x3) |-> f(x1,x2,x3) }
#' \code{ g : (x1,x2) |-> g(x1,x2) }
#' \code{f g} is then \code{ (x1,x2,x3,x4) |-> f(g(x1,x2),x3,x4) }
#' @param m0 morphism
#' @param m1 morphism
#' @return new morphism
morphism.curry.morphism = function(m0, m1) {
    if (m0$defn$type=='Constant') return(m0)
    if (dimension.of.image(m1) > degree.of.freedom(m0)) {
        stop('Can not curry a morphism with low dimensional domain with a morphism having higher dimensional image')
    }
    args0 = generate.formal.arguments.morphism(m0, fsymbol='RGw2G')
    args1 = generate.formal.arguments.morphism(m1, fsymbol='UbsCPF')
    args.new = instantiate.morphism(m1, args1)
    if (!is.list(args.new)) args.new = list(args.new)
    for(i in seq_along(args.new)) {
        args0[[ i ]] = args.new[[ i ]]
    }
    args = c(args1,args0[ -seq_along(args.new) ])
    e = instantiate.morphism(m0, args0)
    morphism.unapply(e, args)
}

#' just an application of \code{explain.expression.as.morphism}
#' instantiate \code{m} using THETAs
#' then explain it
#' might be wrong if THETAs have other meaning in \code{m} than the arguments
#' @param m morphism
#' @return morphisms
explain.morphism.as.morphism = function(m){
    if(dimension.of.image(m)>1) stop('Can only explain one dimensional morphism')
    explain.expression.as.morphism(instantiate.morphism(m))
}

#' atomic matrix to affine.mapping morphism
#' 
#' @param A matrix
#' @param b optional vector
#' @return AffineMapping morphism
convert.matrix.to.morphism = function(A,b=NULL){
    morphism.concrete.new('AffineMapping', convert.matrix.to.affine.mapping(A,b))
}

#' combine a affine mapping and a factor (with same number of levels)
#' generate a affine.mapping.factor.effect
#' @param var.name name
#' @param levels
#' @return morphism
combine.affinemorphism.factor = function(m, var.name, levels) {
    affine = m$defn$defn
    re = list(var.name=as.symbol(var.name), level.name = levels, affine.mapping = affine)
    class(re)='affine.mapping.factor.effect'
    morphism.concrete.new('affine.mapping.factor.effect', re)
}
