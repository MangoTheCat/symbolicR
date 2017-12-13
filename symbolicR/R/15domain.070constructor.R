#' typical values constructor
#' @param grptype group effect
#' @param typical.var typical value name
#' @param last.ind the last used index for \code{THETA}
#' @param defn definition for the target pk variable
#' @param pkvar pk variable name
construct.typical.variable = function(grptype, typical.var, last.ind=0, defn, pkvar){
    i = last.ind + 1
    mappings = sprintf('%s',pkvar)
    foo = NULL
    if(grptype=='multiplicative'){
        foo = list(val= CONS('=', typical.var,
                   CONS('*', CONS('[', quote(THETA), as.numeric(i)),
                             defn)))
    } 
    if(grptype=='additive'){
        foo = list(val= CONS('=', typical.var,
                   CONS('+', CONS('[', quote(THETA), as.numeric(i)),
                             defn)))
    } 
    if(grptype=='affine'){
        foo = list(val= CONS('=', typical.var,
                   CONS('+', defn[[1]], 
                       CONS('*', CONS('[', quote(THETA), as.numeric(i)),
                                 defn[[2]])) ))
    } 
    if(grptype=='exponential'){
        foo = list(val=CONS('=', typical.var,
                   CONS('*', CONS('[', quote(THETA), as.numeric(i) ),
                             CONS('^', defn, CONS('[', quote(THETA), as.numeric(i<-i+1))))))

        mappings = c(mappings, sprintf('PWR%s',pkvar))
    }
    if(grptype=='linear') {
        foo = list(val=CONS('=', typical.var,
                    CONS('+', CONS('[', quote(THETA), as.numeric(i)),
                             CONS('*', defn, CONS('[', quote(THETA), as.numeric(i<-i+1))))))
    }
    if(grptype=='binomial.linear') {
        foo1 = CONS('+', CONS('[', quote(THETA), as.numeric(i)),
                         CONS('*', defn[[2]], CONS('[', quote(THETA), as.numeric(i<-i+1))))
        foo2 = CONS('+', CONS('[', quote(THETA), as.numeric((i<-i+1))),
                         CONS('*', defn[[2]], CONS('[', quote(THETA), as.numeric(i<-i+1))))
        foo = list(val=CONS('=', typical.var,
                CONS('+', CONS('*', defn[[1]], foo1),
                        CONS('*', CONS('-',1,defn[[1]]), foo2))))
    }
    if (is.null(foo)) {
        foo = list(val = CONS('=', typical.var, CONS('[', quote(THETA), as.numeric(i))))
    }
    foo$last.ind = i
    foo$mappings=mappings
    foo
}

#'  extract.coefficient.of.variance
#'
#'  assuming the Var part are of pattern coeff * Var(ETA[i]) or just Var(ETA[i])
#'  for expression \code{coeff * Var(ETA[i])} or \code{coeff * Var(EPS[i])}
#'  for generating: \code{ rnd.part = sqrt(coeff) * ETA[i] }, we need the \code{coeff}
#'
#'  @param Var variance representation
#'  @return found coefficient for standard deviation
extract.coefficient.of.variance = function(Var) {
    if (symbolic.match(quote( Var('?'(rndvar))), Var)) {
        return(1)
    }
    if (symbolic.match(quote( '?'(coeff) * Var('?'(rndvar)) ), Var)) {
        return(simplify.sqrt(CONS('sqrt', coeff)))
    } else if (! (is.atomic(Var) || symbolic.match( quote(Var('?'(rndvar))), Var ))) {
        # atomic or just use the above rnd.part = ETA[i] is acceptable
        # and will further set the ETA[i] initial value as Omega[i,i]
        # Or Else are all unknown variance structures
        stop(sprintf('Unknow variance structure: %s', paste(deparse(Var),collapse='')))
    }
    NULL
}

#' construct.random.variable
#' similar to \code{\link{construct.typical.variable}}
#' @param distrib distribution for the PK variable
#' @param Var variation for the PK variable
#' @param pkvar pk variable name
#' @param typical.var typical variable name
#' @param last.ind last used \code{ETA} index
#' @param random.var random symbol
#' @param exp.mean.var we can use \code{exp.mean.var} to override the construct of it from \code{typical.var}
#' @return construct with expression and index information
construct.random.variable = function(distrib, Var, pkvar, typical.var, last.ind=0, random.var=quote(ETA), 
exp.mean.var=NULL) {
    if (distrib=='fixed' || is.null(Var) || (is.atomic(Var) && Var==0)){
    # actually there is a problem, Var==0 and distrib=log normal, the output is just as is, not exp(mean)
        foo = list(val=CONS('=',pkvar, typical.var),last.ind=last.ind, mappings=character(0))
        return(foo)
    }
    # non vanishing Var
    i = last.ind + 1
    rnd.part = CONS('[', random.var, as.numeric(i))
    mappings = sprintf('%s', pkvar)
    foo = NULL
    if (!is.null(coeff.std<-extract.coefficient.of.variance(Var))) {
        rnd.part = simplify.sqrt(CONS('*', coeff.std, rnd.part))
    }
    if (distrib=='normal'){
        foo = list(val=CONS('=',as.symbol(pkvar), CONS('+',typical.var, rnd.part )))
    } 
    if (distrib=='log normal'){
    # in log normal object, mean is the mean of log(RandomVar), hence we should convert it back
        if (is.null(exp.mean.var)) { exp.mean.var = simplify.collect.log.exp(CONS('exp', typical.var)) }
        foo1 = CONS('*', exp.mean.var, CONS('exp', rnd.part ))
        if (symbolic.match( quote( exp('?'(mult1)) * exp('?'(mult2)) ), foo1)) {
            # better looking: exp(A) * exp(B) -> exp(A + B)
            foo1 = CONS('exp',CONS('+',mult1, mult2))
        }
        foo = list(val=CONS('=',as.symbol(pkvar), foo1 ))
    }
    if (is.null(foo)) {
        stop(sprintf('Dont know how to write distribution: %s ', distrib))
    }
    foo$last.ind = i
    foo$mappings = mappings
    foo
}

#' construct.obs.variable
#' similar to \code{\link{construct.typical.variable}}
#' @param type type of observation variable
#' @param last.ind last used \code{EPS} index
#' @param obs.varname the name for observation variable
#' @param defn definition in the errorObject
construct.obs.variable = function(type, last.ind=0, obs.varname, defn) {
    foo = list()
    if (type=='ruvprop.ruvadd') {
        i = last.ind + 1
        foo$val = CONS('=', as.symbol(obs.varname), 
            CONS('+',
                CONS('*', quote(F) , CONS('exp', CONS('[',quote(EPS), i ) ) ),
                CONS('[', quote(EPS), (i<-i+1))))
        foo$last.ind = i
        foo$mappings = c('RUV_PROP','RUV_ADD')
    } else if (type=='asis') {
        foo$val = CONS('=', as.symbol(obs.varname), defn)
        foo$last.ind = last.ind
        foo$mappings = character(0)
    } else if (type=='normal' || type=='log normal' || type=='fixed') {
        foo = construct.random.variable(type, defn$var, obs.varname, defn$mean, last.ind, random.var=quote(EPS))
    } else {
        stop(sprintf('Do not know how to construct [%s] type errorObject ', type ))
    }
    foo
}
