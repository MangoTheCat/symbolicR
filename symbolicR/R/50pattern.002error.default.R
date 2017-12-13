#' pattern.error.default
#' 
#' turn set of equations define Error block into errorObjs
#'
#' @param eqns equations in \code{$ERROR} block
#' @param md \code{importNmMod} model
#' @return list of definition of observation and \code{IPRED}
#' @author Mango solutions
pattern.error.default = function(eqns=NULL, md = NULL){
    if (is.null(eqns)) eqns = md$Error
    eqns = nonmem.eqns.to.r.eqns(eqns)
    error.socp.filter = function(gg=NULL, variables=NULL, eqns=NULL){
        which(!is.na(match.symbolic(variables,
                            list(   quote(Y), 
                                    quote(IPRED), 
                                    quote(IPRE), 
                                    quote(IRES),
                                    quote(IWRES)))))
    }
    error.rules = list(
        list(pat=function(defn){
            if (symbolic.match( quote(F * exp(EPS['?'(eps1ind)]) + EPS['?'(eps2ind)]) , defn$final.e) ) {
                return(list(type='ruvprop.ruvadd'))
            }
            NULL
        },
        action=function(defn,pat){
                list(critical.varname=quote(Y), type='directly.related.to.nonmem.parameters', defn=pat)
        }))
    L = analyse.0.0( eqns, socp.filter=error.socp.filter, extra.rules=error.rules)
    if(length(L$other.rels)>0) {
        warning(sprintf('Some other relations found in Error Block, will be ignored'))
    }
    errorObjs = list()
    L = L$analysed.defns
    for(component in L){
        for(variable in component){
            errorObjs[[ asCharacterSymbolCall(variable$critical.varname) ]] = variable
        }
    }
    errorObjs
}

#' export.error.default
#' 
#' Export pattern found by \code{\link{pattern.error.default}}
#' 
#' @param errorObjs pattern found by \code{pattern.error.default}
#' @return list of text of \code{$ERROR} block and mappings for \code{EPS}
#' @author Mango solutions
export.error.default = function(errorObjs){
    last.ind = 0
    output.eqns = list()
    eps.mappings = character(0)
    for(variable in errorObjs){
        if(variable$type=='related.on.previous.variables'){
            foo =   construct.obs.variable(type='asis', 
                        last.ind = last.ind,
                        obs.varname=asCharacterSymbolCall(variable$critical.varname), 
                        defn=variable$defn)
            last.ind = foo$last.ind
            output.eqns[[ length(output.eqns)+1 ]] = foo$val
            eps.mappings = c(eps.mappings, foo$mappings)
            next()
        }
        if(variable$type=='directly.related.to.nonmem.parameters'){
            foo =   construct.obs.variable(type=variable$defn$type,
                        last.ind = last.ind,
                        obs.varname=asCharacterSymbolCall(variable$critical.varname),
                        defn=variable$defn)
            last.ind = foo$last.ind
            output.eqns[[ length(output.eqns)+1 ]] = foo$val
            eps.mappings = c(eps.mappings, foo$mappings)
            next()
        }
    }
    o1 = lapply(output.eqns, rexpression.to.nonmem.expression)
    o2 = lapply(o1, compile.ifelse.to.if.statements)
    o3 = sapply(o2, rexpression.to.string.nonmem)
    txt = c('$ERROR')
    txt = c(txt, o3)
    list(txt = paste(txt, collapse='\n'),
         type='errorblock',
         nonmem.parameter.mappings=
         list(from=lapply(seq_along(eps.mappings), function(x) CONS('[',quote(EPS), as.numeric(x))),
              to  =eps.mappings))
}

#' multiple DVs
#' try to match following:
#' \code{ ifelse( DVID == 2, expression1, ifelse(DVID <=1 , expression2, NA)) }
#' @param e expression
#' @param varname variable name
#' @param var.levels possible levels
#' @return list of expressions by var.levels
pattern.multiple.dimensional.by.factor = function(e, varname=quote(DVID), var.levels=as.character(0:2)) {
    L = list()
    # ensure the input are character to avoid the possibility that the input is numerical 0,1,2
    # which fool the L[ var.levels ] code 
    var.levels = as.character(var.levels)
    varlv = var.levels
    walker = function(e0) {
        if (symbolic.match( quote( ifelse( '?'(logic.fun)('?s'(varn),'?a'(lv)) , '?'(e1),'?'(e2) ) ), e0)) {
            if (as.character(logic.fun) %in% c('<','>','<=','>=','==') && varn==varname && lv %in% var.levels) {
                ind = which(eval(CONS(logic.fun, varlv, lv)))
                if (length(ind)==0) {
                    warning('Found logical expression without possible output.')
                } else {
                    for(i in ind) {
                        L[[ as.character(varlv[i]) ]] <<- e1
                    }
                    # remove used ones
                    varlv <<- varlv[ -ind ]
                }
                if (is.atomic(e2) && is.na(e2) && length(varlv)==0) {
                # the last if(..., ..., NA) form, and assume all possibility are counted
                    return(TRUE)
                }
                return(Recall(e2))
            }
        }
        FALSE
    }
    if (walker(e)) {
        re = list(selecting.variable=varname, levels=symbolic.CMatrix(L[var.levels]))
        class(re)='multiple.dependent.variable'
        return(re)
    }
    NULL
}

multiple.dependent.variable.to.disturbation = function(m) {
    le = m$levels
    gg = create.abstractGraph.from.equations(le)
    syms = gg$vertex.table
    eps.ind = which(sapply(syms, function(x) symbolic.match(quote(EPS['?a'(III)]), x) ) )
    backref.eps = NULL
    levels = names(le)
    if (length(eps.ind)>0) {
        eps.sym = syms[ eps.ind ]
        # eps are EPS[i] 's i
        eps = sort(unique(sapply(eps.sym, function(x) x[[3]])))
        backref.eps = eps
        re = random.disturbation.new(le, fixed = NULL, 
            random = lapply(eps, function(x) CONS('[', quote(EPS), as.numeric(x))))
    } else {
        re = le
    }
    attr(re,'backref')=backref.eps
    attr(re,'selecting.variable') = m$selecting.variable
    attr(re,'levels') = levels
    if (!is(re,'multiple.dependent.variable')) {
        class(re) = c('multiple.dependent.variable', class(re))
    }
    re
}

#' apply a multiple.dependent.variable pattern
#' @param m multiple.dependent.variable object
#' @return expression
multiple.dependent.variable.apply = function(m, ...) {
    m$selecting.variable = attr(m,'selecting.variable')
    level.names = attr(m,'levels')
    if (is(m,'random.disturbation')) {
        m$levels = disturbation.apply(m,...)
    } else if (is(m,'morphism')) {
        m$levels = instantiate.morphism(m, ...)
    }
    if (is.null(level.names)) {
        level.names = names(m$levels)
    }
    v1 = vector(length(level.names),mode='list')
    # try convert into number if possible
    for(i in seq_along(level.names)) {
        tmp = suppressWarnings(as.numeric(level.names[i]))
        if (is.na(tmp)) tmp = level.names[i]
        v1[[i]]=CONS('==', m$selecting.variable, tmp)
    }
    symbolic.inner.product.logical(v1, m$levels)
}

#' S3 print function
print.multiple.dependent.variable = function(m0){
    m0$selecting.variable = attr(m0,'selecting.variable')
    if (is(m0,'random.disturbation')) {
        m0$levels = disturbation.apply(m0)
    } else if (is(m0,'morphism')) {
        m0$levels = instantiate.morphism(m0, generate.formal.arguments.morphism(m0))
    }
    re = sprintf('Multiple variable selected by %s', paste(deparse(m0$selecting.variable),collapse=''))
    if (is(m0,'random.disturbation')) {
        re = sprintf('%s, also disturbed by %s', re, print.distribution.common0(m0$distribution))
    }
    a = matrix(sapply(m0$levels, function(x) paste(deparse(x),collapse='')))
    a = as.table(cbind('[',a,']'))
    colnames(a) = rep('', NCOL(a))
    rownames(a) = sprintf('%s.%s', paste(deparse(m0$selecting.variable),collapse=''),names(m0$levels))
    b = capture.output(print(a))
    foo = paste(b[-grep('^\\s*$',b)], collapse='\n')
    re = sprintf('%s\n%s',re,foo)
    cat(re,'\n')
}
