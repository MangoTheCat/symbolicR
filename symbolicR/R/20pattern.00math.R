# This file contains some mathematical patterns for future use.

#' binomial
#' assume \code{X ~ binom(Q)}
#' \code{P(X=X1) = Q1}
#' \code{P(X=X2) = 1-Q1}
#' then \code{E(X) = Q1 * X1 + (1-Q1) * X2}
#' i.e. \code{E(X)} (binomial combination) is the expection of a bernoulli distribution
#' @param e expression
#' @return if matched, the combination is viewed as a random variable, and return the discription of the random variable
pattern.binomial.combination = function(e) {
    check.add.to.one.pair = function(le, ind){
        a = ind[1]
        b = ind[2]
        left.ind = setdiff(1:4,ind)
        if(simplify.2(CONS('+',le[[a]] , le[[b]]))==1) {
            p = le[[a]]
            if (object.size(le[[a]]) > object.size(le[[b]])) {
                p = le[[b]]
                left.ind = rev(left.ind)
            }
            return(list(p=p, val=le[left.ind]))
        }
        else return(NULL)
    }

    if (symbolic.match(quote('?'(AA) * '?'(BB) + '?'(CC) * '?'(DD) ) , e) ) {
        le = c(AA, BB, CC, DD)
        for(i in c(1,2)){
            for(j in c(3,4)){
                if (!is.null(pat<-check.add.to.one.pair(le, c(i,j)))) {
                    return(pat)
                }
            }
        }
    }
    pat = pattern.convex.combination(e)
    if (!is.null(pat) && length(pat$val)==2) {
        pat$p = pat$distrib[[1]]
        pat$distrib=NULL
        return(pat)
    }
    NULL
}

#' similar to bernoulli 
#' see \code{\link{pattern.binomial.combination}}
#' however, \code{Xi} is now thought like a multinomial distribution
#' @param e expression
#' @return \code{NULL} if not matched, random variable's distribution if matched
pattern.convex.combination = function(e){
    le = expand.as.sum.of.product(e)
    s0 = 0
    for(x in le){
        s0 = simplify.2(CONS('+', x[[1]], s0))
    }
    if (s0 == 1) {
        distrib = lapply(le, function(x) x[[1]])
        xi = lapply(le, function(x) expression.from.list.of.product(x[-1]))
        return(list(distrib=distrib, val=xi))
    }
    NULL
}
