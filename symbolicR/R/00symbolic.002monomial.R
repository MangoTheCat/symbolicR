#' monomial.of.random.variable.from.list.of.product
#'
#' For a list of expressions, convert it as monomial of \code{ETA[i]}'s\cr
#' \code{[ coefficients, ETA[i] ] ==> coefficients * ( ETA[i]^di .... ) }
#'
#' @param le list of product
#' @return list of length 2, like \code{ [coefficients, monomial]}
#' @author Mango solutions
monomial.of.random.variable.from.list.of.product = function(le){
    N = length(le)
    coeffs = list()
    Ncoeff = 0
    mono = list() 
    for(i in 1:N){
        if (symbolic.match(quote( '?s'(rndname)['?a'(rndind)] ), le[[i]]) &&
            (rndname=='ETA' || rndname=='EPS')) {
            var.name = sprintf('%s[%s]', as.character(rndname), rndind )
            if (is.null(mono[[ var.name ]])) {
                mono[[ var.name ]] = list(
                    sym = CONS('[',rndname,rndind),
                    exp = 0
                )
            }
            mono[[ var.name ]]$exp = mono[[ var.name ]]$exp + 1
        } else if ( symbolic.match(quote( '?s'(rndname)['?a'(rndind)] ^ '?a'(rndexp) ), le[[i]]) &&
            (rndname=='ETA' || rndname=='EPS')) {
            var.name = sprintf('%s[%s]', as.character(rndname), rndind )
            if (is.null(mono[[ var.name ]])) {
                mono[[ var.name ]] = list(
                    sym = CONS('[',rndname,rndind),
                    exp = 0
                )
            }
            mono[[ var.name ]]$exp = mono[[ var.name ]]$exp + rndexp
        } else {
            coeffs[[ Ncoeff<-Ncoeff+1 ]] = le[[i]]
        }
    }
    if (length(mono)==0) {
        mono1 = NULL
    } else {
        mono1 = NULL
        for(i in 1:length(mono)){
            if (mono[[i]]$exp==0) next
            if (mono[[i]]$exp==1) {
                tmp = mono[[i]]$sym
            } else {
            # exponential can be negative
                tmp = CONS('^',mono[[i]]$sym, mono[[i]]$exp)
            } 
            mono1 = c(mono1, tmp)
        }
        if (!is.null(mono1)) mono1 = sort.by.object.size(mono1)
    }
    list(expression.from.list.of.product(coeffs),  expression.from.list.of.product(mono1))
}

#' create.polynomial.of.random.variable
#'
#' Convert an expression \code{e} to a polynomial of \code{ETA[i]}'s or \code{EPS[i]}'s as following:\cr
#' \code{ e --> [ [ coeff, main.monomial] ] }
#'
#' @param e expression
#' @return list of list representing sum of list like \code{ [ [coeff1, EPS[1]^1], [coeff2, EPS[1]^2] ] }
#' @author Mango solutions
create.polynomial.of.random.variable = function(e){
    le = expand.as.sum.of.product(e)
    le = lapply(le, monomial.of.random.variable.from.list.of.product)
    N = length(le)
    if (N==1) return(le)
    key.ind = c(1)
    val.inds = list( c(1) )
    for(i in 2:N){
        found=F
        for(j in 1:length(key.ind)) {
            if(equal.symbol.light(le[[i]][[2]],le[[ key.ind[j] ]][[2]])){
                val.inds[[ j ]] = c(val.inds[[j]], i )
                found = T
                break
            }
        }
        if(!found){
            key.ind = c(key.ind,  i)
            val.inds[[ length(key.ind) ]] = c(i)
        }
    }
    N1 = length(key.ind)
    # now combine them
    re = list()
    for(i in 1:N1){
        e1 = simplify.1(expression.from.list.of.sum(lapply( val.inds[[ i ]] , function(x) le[[ x ]][[1]] )))
        e2 = le[[ key.ind[i] ]][[2]]
        re[[i]] = list(e1, e2)
    }
    ord1 = sapply(re, function(x) object.size(x[[2]]))
    ord2 = sapply(re, function(x) deparse(x[[2]]))
    re[ order(ord1,ord2) ]
}

#' pattern.distrib.normal
#'
#' used by higher level pattern matcher, if not matched, then try next distribution pattern
#'
#' try to match a expression to normal distribution, e.g.\cr
#' \code{ expressions + coeffs*ETA[i] } \cr
#' \code{ expressions + coeffs*EPS[i] }
#'
#' the \code{\link{create.polynomial.of.random.variable}} are used here for combining \code{ETA[i]}'s or \code{EPS[i]}'s , for example \cr
#' \code{ TVCL + ETA[1] + 100 * ETA[1] + TVCL*ETA[1] } should be coverted to \code{ TVCL + ( 101 + TVCL)*ETA[1] } \cr
#' and then mean is \code{TVCL} and covariance is \code{(101+TVCL)^2 * Cov(ETA[1])}  \cr
#'
#' @param e expression
#' @param ignore.non.atomic.coefficient do not match if the coefficient random varible is non-atomic
#' @return if matched, a list indicating the distribution, including mean and covariance; if not matched then \code{NULL}
#' @author Mango solutions
pattern.distrib.normal = function(e, ignore.non.atomic.coefficient=FALSE){
    le = create.polynomial.of.random.variable(e)
    # le : [ [ coeff, ETA[i]^di * ETA[j]^dj ] , [ ] ]
    if (length(le)==1) {
        if (symbolic.match(quote( '?'(rndvar)[ '?'(rndind) ] ), le[[1]][[2]] )){
            if (ignore.non.atomic.coefficient && !is.atomic(le[[1]][[1]])) {
                return(NULL)
            }
            distrib = list(type='normal',
                mean = 0,
                var  = simplify.2(
                    CONS('*',
                            CONS('^', le[[1]][[1]],2),
                            CONS('Var',le[[1]][[2]]))))
            return(distrib)
        }
    }
    if (length(le)==2) {
    # assume  [ [ coeff, 1] , [coeff, ETA[i] ] ]
        if (is.atomic(le[[1]][[2]]) && le[[1]][[2]]==1 &&
            symbolic.match(quote( '?'(rndvar)[ '?'(rndind) ] ), le[[2]][[2]] )){
            if (symbolic.grep(quote( EPS['?'(anyind)] ), le[[1]][[1]]) ||
                symbolic.grep(quote( ETA['?'(anyind)] ), le[[1]][[1]]) ){
                # avoid posibility of following:
                #   exp(ETA[i]) + ETA[i]
                return(NULL)
            }
            if (ignore.non.atomic.coefficient && !is.atomic(le[[2]][[1]])) {
                return(NULL)
            }
            distrib = list(type='normal',
                mean = le[[1]][[1]],
                var  = simplify.2(
                            CONS('*',
                                CONS('^',le[[2]][[1]],2),
                                CONS('Var',le[[2]][[2]]))))
            return(distrib)
        }
    }
    NULL
}

#' monomial.of.exp.of.random.variable.from.list.of.product
#'
#' similar to \code{\link{monomial.of.random.variable.from.list.of.product}} \cr
#' For a list of expressions, convert it as monomial of \code{eps(ETA[i])}'s\cr
#' \code{[ coefficients, eps(ETA[i]) ] ==> coefficients * ( ETA[i]^di .... ) }
#'
#' @param le list of product
#' @return list of length 2, like \code{ [coefficients, monomial]}
#' @author Mango solutions
monomial.of.exp.of.random.variable.from.list.of.product = function(le){
    N = length(le)
    coeffs = list()
    Ncoeff = 0
    mono = list() 
    for(i in 1:N){
        # exp( ETA[i] )
        if ((symbolic.match(quote( exp('?s'(rndname)['?a'(rndind)]) ), le[[i]])) &&
            (rndname=='ETA' || rndname=='EPS')) {
            var.name = sprintf('%s[%s]', as.character(rndname), rndind )
            if (is.null(mono[[ var.name ]])) {
                mono[[ var.name ]] = list(
                    sym = le[[i]],
                    exp = 0
                )
            }
            mono[[ var.name ]]$exp = mono[[ var.name ]]$exp + 1
        # exp( ETA[i] ) ^ j
        } else if ((symbolic.match(quote( exp('?s'(rndname)['?a'(rndind)]) ^ '?a'(rndexp) ), le[[i]])) &&
            (rndname=='ETA' || rndname=='EPS')) {
            var.name = sprintf('%s[%s]', as.character(rndname), rndind )
            if (is.null(mono[[ var.name ]])) {
                mono[[ var.name ]] = list(
                    sym = CONS('exp', CONS('[',rndname,rndind)),
                    exp = 0
                )
            }
            mono[[ var.name ]]$exp = mono[[ var.name ]]$exp + rndexp
        # exp( c * ETA[i] )
        } else if ((symbolic.match(quote( exp('?s'(rndname)['?a'(rndind)] * '?a'(rndexp)) ), le[[i]]) ||
                    symbolic.match(quote( exp('?a'(rndexp) * '?s'(rndname)['?a'(rndind)] )), le[[i]]) )
                     && (rndname=='ETA' || rndname=='EPS')) {
            var.name = sprintf('%s[%s]', as.character(rndname), rndind )
            if (is.null(mono[[ var.name ]])) {
                mono[[ var.name ]] = list(
                    sym = CONS('exp', CONS('[',rndname,rndind)),
                    exp = 0
                )
            }
            mono[[ var.name ]]$exp = mono[[ var.name ]]$exp + rndexp
        } else {
            coeffs[[ Ncoeff<-Ncoeff+1 ]] = le[[i]]
        }
    }
    if (length(mono)==0) {
        mono1 = NULL
    } else {
        mono1 = NULL
        for(i in 1:length(mono)){
            if (mono[[i]]$exp==0) next
            if (mono[[i]]$exp==1) {
                tmp = mono[[i]]$sym
            } else {
            # positive or negative
                tmp = CONS('^',mono[[i]]$sym, mono[[i]]$exp)
            }
            mono1 = c(mono1, tmp)
        }
        if (!is.null(mono1)) mono1 = sort.by.object.size(mono1)
    }
    list(expression.from.list.of.product(coeffs),  expression.from.list.of.product(mono1))
}

#' create.polynomial.of.exp.of.random.variable
#' 
#' similar to \code{\link{create.polynomial.of.random.variable}} \cr
#' while monomials here are of form \code{exp(ETA[i])^di} \cr
#' \code{e --> [ [ coeff, main.monomial] ]}
#'
#' @param e expression
#' @return list of list \code{[[coeff, monomial]]}
#' @author Mango solutions
create.polynomial.of.exp.of.random.variable = function(e){
    le = expand.as.sum.of.product(e)
    le = lapply(le, monomial.of.exp.of.random.variable.from.list.of.product)
    N = length(le)
    if (N==1) return(le)
    key.ind = c(1)
    val.inds = list( c(1) )
    for(i in 2:N){
        found=F
        for(j in 1:length(key.ind)) {
            if(equal.symbol.light(le[[i]][[2]],le[[ key.ind[j] ]][[2]])){
                val.inds[[ j ]] = c(val.inds[[j]], i )
                found = T
                break
            }
        }
        if(!found){
            key.ind = c(key.ind,  i)
            val.inds[[ length(key.ind) ]] = c(i)
        }
    }
    N1 = length(key.ind)
    # now combine them
    re = list()
    for(i in 1:N1){
        e1 = simplify.1(expression.from.list.of.sum(lapply( val.inds[[ i ]] , function(x) le[[ x ]][[1]] )))
        e2 = le[[ key.ind[i] ]][[2]]
        re[[i]] = list(e1, e2)
    }
    ord1 = sapply(re, function(x) object.size(x[[2]]))
    ord2 = sapply(re, function(x) deparse(x[[2]]))
    re[ order(ord1,ord2) ]
}

#' pattern.distrib.lognormal
#'
#' similar to \code{\link{pattern.distrib.normal}} \cr
#' For \code{ A = B * exp(ETA[i])} , \cr
#' extract means of \code{log(A)} is \code{log(B)} \cr
#' and  covariance of \code{log(A)} is \code{Cov(ETA[i])}
#'
#' @param e expression
#' @param ... options passed to \code{\link{pattern.distrib.normal}}
#' @return if matched, return the distribution; else \code{NULL}
#' @author Mango solutions
pattern.distrib.lognormal = function(e, ...) {
    log.e = simplify.2( as.call( c(as.symbol('log'), e) ) )
    pt = pattern.distrib.normal(log.e)
    if (!is.null(pt)) {
        pt$type='log normal'
        return(pt)
    }
    e2 = create.polynomial.of.exp.of.random.variable(e)
    if (length(e2)==1){
        e = as.call(c(as.symbol('*'), e2[[1]][[1]], e2[[1]][[2]]))
        log.e = simplify.2( as.call( c(as.symbol('log'), e) ) )
        pt = pattern.distrib.normal(log.e, ...)
        if (!is.null(pt)) {
            pt$type='log normal'
            return(pt)
        }
    }
    NULL
}

#' pattern.distrib.fixed
#'
#' match a expression which not having ETA's
#' @param e expression
#' @return mean and vars if matched; else \code{NULL}
pattern.distrib.fixed = function(e) {
    if (!expressions.has.random.variables(e)) {
        return(list(type='fixed',mean=e ,var=0))
    }
    NULL
}

#' monomial.of.theta.from.list.of.product
#'
#' similar to \code{\link{monomial.of.random.variable.from.list.of.product}} \cr
#' For a list of expressions, convert it as monomial of \code{eps(ETA[i])}'s\cr
#' \code{[ coefficients, THETA[i] ] ==> coefficients * ( THETA[i]^di .... ) }
#'
#' @param le list of product
#' @return list of length 2, like \code{ [coefficients, monomial]}
#' @author Mango solutions
monomial.of.theta.from.list.of.product = function(le){
    N = length(le)
    coeffs = list()
    Ncoeff = 0
    mono = list() 
    for(i in 1:N){
        if ((symbolic.match(quote( '?'(rndname)['?'(rndind)] ), le[[i]])) &&
            rndname=='THETA' ) {
            var.name = sprintf('%s[%s]', as.character(rndname), rndind )
            if (is.null(mono[[ var.name ]])) {
                mono[[ var.name ]] = list(
                    sym = CONS('[',rndname,rndind),
                    exp = 0
                )
            }
            mono[[ var.name ]]$exp = mono[[ var.name ]]$exp + 1
        } else if ((symbolic.match( quote( '?'(rndname)['?'(is.numeric,rndind)] ^ '?'(is.numeric,rndexp) ), 
            le[[i]])) && rndname=='THETA' ) {
            var.name = sprintf('%s[%s]', as.character(rndname), rndind )
            if (is.null(mono[[ var.name ]])) {
                mono[[ var.name ]] = list(
                    sym = CONS('[',rndname,rndind),
                    exp = 0
                )
            }
            mono[[ var.name ]]$exp = mono[[ var.name ]]$exp + rndexp
        } else {
            coeffs[[ Ncoeff<-Ncoeff+1 ]] = le[[i]]
        }
    }
    if (length(mono)==0) {
        mono1 = NULL
    } else {
        mono1 = NULL
        for(i in 1:length(mono)){
            if (mono[[i]]$exp==0) next
            if (mono[[i]]$exp==1) {
                tmp = mono[[i]]$sym
            } else {
                # position or negative
                tmp = as.call(c(as.symbol('^'),mono[[i]]$sym, mono[[i]]$exp))
            }
            mono1 = c(mono1, tmp)
        }
        if (!is.null(mono1)) mono1 = sort.by.object.size(mono1)
    }
    list(expression.from.list.of.product(coeffs),  expression.from.list.of.product(mono1))
}

#' create.polynomial.of.theta
#'
#' similar as \code{\link{create.polynomial.of.random.variable}} \cr
#' \code{ e --> [ [ coeff, main.monomial] ] } \cr
#' Only difference is monomial here are of type \code{ THETA[i]^di }
#' 
#' @param e expression
#' @return list of list as \code{ [ [coeff, main.monomial] ,...] }
#' @author Mango solutions
create.polynomial.of.theta = function(e){
    le = expand.as.sum.of.product(e)
    le = lapply(le, monomial.of.theta.from.list.of.product)
    N = length(le)
    if (N==1) return(le)
    key.ind = c(1)
    val.inds = list( c(1) )
    for(i in 2:N){
        found=F
        for(j in 1:length(key.ind)) {
            if(equal.symbol.light(le[[i]][[2]],le[[ key.ind[j] ]][[2]])){
                val.inds[[ j ]] = c(val.inds[[j]], i )
                found = T
                break
            }
        }
        if(!found){
            key.ind = c(key.ind,  i)
            val.inds[[ length(key.ind) ]] = c(i)
        }
    }
    N1 = length(key.ind)
    # now combine them
    re = list()
    for(i in 1:N1){
        e1 = simplify.1(expression.from.list.of.sum(lapply( val.inds[[ i ]] , function(x) le[[ x ]][[1]] )))
        e2 = le[[ key.ind[i] ]][[2]]
        re[[i]] = list(e1, e2)
    }
    ord1 = sapply(re, function(x) object.size(x[[2]]))
    ord2 = sapply(re, function(x) deparse(x[[2]]))
    re[ order(ord1,ord2) ]
}

