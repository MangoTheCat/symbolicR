#' analyse.group.effect
#'
#' for a expression, determine if there is a recognizable structure of the mean expression
#'
#' the slot \code{parameter.space} in return object is for back referring to the \code{THETA[i]}
#' what can be used for look up the definition domain of the parameters
#'
#' @param e expression
#' @return if the expression has a known pattern as multiplicative, additive, affine and exponential
analyse.group.effect = function(e) {
    le = create.polynomial.of.theta(e)
    if (length(le)==1 && symbolic.match(quote( THETA['?'(XXXX)]), le[[1]][[2]])
        && (!symbolic.grep(quote(THETA['?'(ANYI)]),le[[1]][[1]]))) {
        if (le[[1]][[1]]==1) { 
        # p[1] -> p[1]
            return(list(grptype='identicalmapping',dim=1, parameter.space=XXXX))
        }
        # p[1] -> grpvar * p[1]
        return(list(grptype='multiplicative', grpvar.definition=le[[1]][[1]], parameter.space=XXXX))
    }
    if (length(le)==2 && le[[1]][[2]]==1 && !symbolic.grep(quote(THETA['?a'(ANYI)]), le[[1]][[1]]) &&
        symbolic.match(quote(THETA['?'(XXXX)]), le[[2]][[2]])) {
        if (is.atomic(le[[2]][[1]]) && le[[2]][[1]]==1) {
            # if theta has a coefficient , then it's affine
            # p[1] -> grpvar + p[1]
            return(list(grptype='additive', grpvar.definition=le[[1]][[1]], parameter.space=XXXX))
        }
        # two covariate but one theta
        # p[1] -> grpvar[1] + grpvar[2] * p[1]
        return(list(grptype='affine', grpvar.definition=c(le[[1]][[1]],le[[2]][[1]]), parameter.space=XXXX))
    }
    if (length(le)==2 && le[[1]][[1]]==1) {
        if (symbolic.match(quote(THETA['?'(XXXX)]), le[[1]][[2]]) &&
            symbolic.match(quote(THETA['?'(YYYY)]), le[[2]][[2]])) {
            # one covariate but two thetas(one for intercept, other for coefficient)
            # p[1],p[2] -> p[1] + grpvar * p[2]
            return(list(grptype='linear', grpvar.definition=le[[2]][[1]], parameter.space=c(XXXX,YYYY)))
        }
    }
    if (symbolic.match( quote( THETA['?a'(X)] * '?'(Y) ^ THETA['?a'(Z)] ), e)){
        # p[1], p[2] ->  p[1] * grpvar^ p[2]
        return(list(grptype='exponential', grpvar.definition=simplify(Y), parameter.space=c(X,Z)))
    }
    if (!is.null(pat<-pattern.linear.by.binomial.categorical(e))) return(pat)
    NULL
}

#' match a convex combination of two similar linear pattern
#' \code{ p[1],p[2],p[3],p[4] -> grpvar[1] * (p[1] + grpvar[2]*p[2]) + (1-grpvar[1]) * (p[3] + grpvar[2]*p[4]) }
#' @param e expression
#' @return pattern description
pattern.linear.by.binomial.categorical = function(e){
    pat = pattern.binomial.combination(e)
    if (is.null(pat)) return(NULL)
    # if (!is.binomial.categorical(pat$p)) return(NULL) 
    x1 = analyse.group.effect(pat$val[[1]])
    x2 = analyse.group.effect(pat$val[[2]])
    if (length(intersect(x1$parameter.space,x2$parameter.space))==0 &&
        x1$grptype=='linear' && x2$grptype=='linear' && x1$grpvar.definition==x2$grpvar.definition) {
        re = list(grptype='binomial.linear', grpvar.definition=c(pat$p,x1$grpvar.definition), parameter.space=c(x1$parameter.space,x2$parameter.space))
        return(re)
    }
    NULL
}

#' reverse direction of \code{simplify.2} for \code{log} and \code{exp} expressions
simplify.collect.log.exp.rules = character.rules.to.list.rules(c(
" '?s'(X) * log('?'(Y)) "  ,        " log( ':'(Y)^':'(X) ) ",
" '?a'(X) * log('?'(Y)) "  ,        " log( ':'(Y)^':'(X) ) ",
" '?s'(X)['?'(XX)] * log('?'(Y)) "  ,        " log( ':'(Y)^':'(X)[':'(XX)] ) ",
"log('?'(X)) + log('?'(Y))",            "log(':'(X) * ':'(Y))",
"exp(log('?'(X)))" ,                 "':'(X)",
"exp('?'(X) + '?'(Y))",             " exp(':'(X)) * exp(':'(Y)) "))
#' simplify.collect.log.exp
#' convert \code{ exp(log(A) + C*log(B)) } to \code{ A * B^C }
#' @param e expression
simplify.collect.log.exp = symbolic.simplify.gigo(simplify.collect.log.exp.rules)
