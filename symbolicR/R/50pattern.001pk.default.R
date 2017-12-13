#' pattern.pk.default
#'
#' Try to match typical control file sytle if we're using NONMEM PREDPP subroutines \cr
#' Assume control stream is PK parameters orientated, i.e. we know the whold model given we know definitions of all PK parameters
#' 
#' For each connected component return by \code{\link{analyse.0.0}}, check if it can be explained by
#' simple relations as additive, multiplicative, affine, linear or binomial.linear combination of parameters and covariates
#'
#' @param eqns equations
#' @param md model returned by \code{importNmMod} , if \code{NULL}, the covariance are left as is, not as numeric values
#' @param L the analysed result by \code{analyse.0.0}, override the \code{md} and \code{eqns}
#' @return list of definition of PK parameters, including their distributions, how the mean (including group or covariant effect) and variance are defined, e.t.c
#' @author Mango solutions
pattern.pk.default = function(eqns=NULL, md = NULL, L=NULL){
    if (!(!is.null(L) && is.list(L) && !is.null(L$analysed.defns))) {
        if (is.null(eqns)) eqns = md$PK
        eqns = nonmem.eqns.to.r.eqns(eqns)
        L = analyse.0.0(eqns, socp.filter=socp.filter.known.pattern)
    }
    if (length(L$other.rels)>0) {
    # more complicated than default
    # NULL means not matched
        warning('other.rels not yet supported by pattern.pk.default')
        return(NULL)
    }
    L = L$analysed.defns
    pkObjs = list()
    for(i in L){
        for(j in i){
            foo = list()
            foo$pkvar  = asCharacterSymbolCall(j$critical.varname)
            foo$type=j$type
            foo$has.random.effect=FALSE
            if (j$type=='directly.related.to.nonmem.parameters') {
                if(!(j$defn$type=='fixed' || equal.symbol.light(j$defn$var, 0))) {
                    foo$has.random.effect=TRUE
                    foo$distrib = j$defn$type
                    foo$Var = eval.nonmem.parameter(j$defn$var,md=md)
                }
                # if we can't create a morphism
                # or the morphism is just constant w.r.t THETA
                tmp = explain.expression.as.morphism(j$defn$mean)
                if (is.null(tmp) || tmp$defn$type=='Constant') foo$mean = j$defn$mean
                else foo$mean = tmp
            }
            if (j$type=='related.on.previous.variables') {
                foo$defn=j$defn
            }
            if (j$type=='unknown') {
                foo$defn=j$defn
            }
            pkObjs[[ foo$pkvar ]] = foo
        }
    }
    pkObjs
}

#' export.pk.default
#'
#' Exporter for the above pattern found by \code{\link{pattern.pk.default}} 
#' 
#' @param pkObjs pattern recognized by \code{pattern.pk.default}
#' @return list containing the text for \code{$PK} block and the mappings for \code{THETA} and \code{ETA}
#' @author Mango solutions
export.pk.default = function(pkObjs){

    if (missing(pkObjs)) { stop("pkObjs is missing") }
    
    # Assign THETA[] to POP_VAR
    # Assign possible Covariants
    theta.mappings = character(0)
    theta.ranges = list()
    theta.backref=integer(0)
    eta.mappings = character(0)
    last.eta.ind = 0
    last.theta.ind=0
    # Output eqns
    output.eqns = list()
    
    for(i in seq_along(pkObjs)){
    
    # TV%s   : typical values, affected by maybe covariant variables
        if (pkObjs[[i]]$type %in% c('related.on.previous.variables','unknown')){
            # just put to the last
            output.eqns = c(output.eqns, CONS('=', as.symbol(pkObjs[[i]]$pkvar), pkObjs[[i]]$defn))
            next()
        }
        ### type == 'directly.related.to.nonmem.parameters'
        if (!is(pkObjs[[i]]$mean,'morphism')) {
            # not a morphism, should be a expression or symbol, wrap it as morphism
            morphism = morphism.concrete.new('Constant', val= pkObjs[[i]]$mean, domain=quote(RealNumber^1))
        } else {
            # Assume mean is defined by morphism
            morphism = pkObjs[[i]]$mean
        }
        morphism.dim = degree.of.freedom(morphism)
        parameter.pool = list()
        if (morphism.dim>0) {
            parameter.pool=lapply( seq(from=last.theta.ind+1, to = last.theta.ind + morphism.dim ),
                function(x) CONS('[',quote(THETA), as.numeric(x)) )
            last.theta.ind = last.theta.ind + morphism.dim 
            theta.mappings = c(theta.mappings, 
                sapply(1:morphism.dim, function(x) sprintf('%s[%s]', pkObjs[[i]]$pkvar, x)))
            theta.ranges = c(theta.ranges, as.low.upper.list.morphism.domain(morphism))
            # in outer program, code can use this information to look up initial estimations
            # from PARAMETER object
            theta.backref = c(theta.backref, backref.parameter.of.morphism(morphism) )
        }
        morphism.expression = instantiate.morphism(morphism, parameter.pool=parameter.pool)
        ### Now, we assume morphism.expression is one-dimensional, i.e. just a expression
        if (pkObjs[[i]]$has.random.effect) {
            # the following variable is used for pretty output of TVKA * EXP(ETA(1))
            # if using plain structure, will output 
            # TVKA = log(THETA(1))
            # KA = EXP(TVKA + ETA(1))
            exp.mean.var = NULL
            if (is.symbol(morphism.expression)) {
            # if just MU_1 + ETA[i] type, do not need to add a typical variable
                typical.var = morphism.expression
            } else {
            # We added typical variable for better looking here if morphism.expression is complex
                typical.var = as.symbol(sprintf('TV%s', pkObjs[[i]]$pkvar))
            # here, we need to deal with log normal specially to make a better looking
                if (pkObjs[[i]]$distrib=='log normal') {
                    if (symbolic.match(quote(log('?'(XXXX))), morphism.expression)) {
                        morphism.expression = XXXX
                        exp.mean.var = typical.var
                    }
                }
                output.eqns = c(output.eqns, CONS('=', typical.var, morphism.expression) )
            }
            ## currently, random part processing and structural part processing framework are not 
            ## consistent, we use morphism for mean and directly process the random part
            ## Maybe changed in future to include everything in morphism
            ## random part
            foo = construct.random.variable(
                pkObjs[[i]]$distrib, 
                pkObjs[[i]]$Var, 
                pkObjs[[i]]$pkvar, 
                typical.var, last.eta.ind,
                exp.mean.var=exp.mean.var)
            output.eqns = c(output.eqns, foo$val)
            last.eta.ind = foo$last.ind
            eta.mappings = c(eta.mappings, foo$mappings)
        } else {
        # Only need the pk part
            output.eqns = c(output.eqns, CONS('=', as.symbol(pkObjs[[i]]$pkvar), morphism.expression ))
        }
    }
    o1 = lapply(output.eqns, rexpression.to.nonmem.expression)
    o2 = list()
    # should not use lapply here, because the output of compile.ifelse.to.if.statements might be list
    # we expect list appending, not nested
    for( obj in o1 ) {
        o2 = c(o2, compile.ifelse.to.if.statements(obj))
    }
    o3 = sapply(o2, rexpression.to.string.nonmem)
    pkblock = c('\n$PK',o3)
    pkblock = paste(pkblock, collapse='\n')
    list(txt=pkblock, 
         type='pkblock',
         nonmem.parameter.mappings=list(
            from=c(lapply(seq_along(theta.mappings), function(x) CONS('[',quote(THETA), as.numeric(x))),
                   lapply(seq_along(eta.mappings), function(x) CONS('[',quote(ETA), as.numeric(x)))),
            to = c(theta.mappings, eta.mappings),
            theta.ranges=theta.ranges,
            theta.backref=theta.backref
            ))
}

#' simplify output equations to NONMEM string
#' @param eqns equations
#' @return one string
reqns.to.nonmem.strings = function(eqns){
    o0 = lapply(eqns, simplify.pretty.01)
    o1 = lapply(o0, rexpression.to.nonmem.expression)
    o2 = list()
    # should not use lapply here, because the output of compile.ifelse.to.if.statements might be list
    # we expect list appending, not nested
    for( obj in o1 ) {
        o2 = c(o2, compile.ifelse.to.if.statements(obj))
    }
    o3 = sapply(o2, rexpression.to.string.nonmem)
    paste(o3,collapse='\n')
}
