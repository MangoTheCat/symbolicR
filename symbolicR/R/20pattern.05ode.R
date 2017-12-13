simplify.pretty.01.rule = list(
    list(CONS('-',1), as.numeric(-1)) , # symbolic to atomic
    list(CONS('*',quote('?'(X)), CONS('^', quote('?'(Y)), -1)), quote( ':'(X)/':'(Y) )),
    list(CONS('*',-1,quote('?'(X))), CONS('-', quote(':'(X)))),
    list(quote(-'?'(X)*'?'(Y)),CONS('-',quote(':'(X)*':'(Y)))),
    list(quote(-'?'(X)/'?'(Y)),CONS('-',quote(':'(X)*':'(Y)))),
    list(quote('?'(X) + -'?'(Y)), quote(':'(X) - ':'(Y) ) )
)

#' simplify.pretty.01
#'
#' translate \code{CONS('-',1)} to atomic -1
#' \code{ X^Y^-1 } to \code{ X/Y }
#' \code{ -X * Y } to \code{ -(X*Y) } which can be print as \code{-X*Y} using \code{\link{rexpression.to.string.nonmem}}
#'
#' @param e expression
#' @return simplied expression
simplify.pretty.01 = symbolic.simplify.gigo( simplify.pretty.01.rule)

#' create.firstorder.linear.ode
#'
#' ODE only consider very simple ones as following
#' \code{ dx/dt = A \%*\% x }
#'
#' @param state.vector a list of symbols or \code{TUPLE(x1,x2,x3)} like expression which will be converted into list automatically
#' @param A matrix should be list of list type
#' @return ODE object
create.firstorder.linear.ode = function(state.vector, A){
    if (is.call(state.vector)) {
        if (state.vector[[1]]=='TUPLE') state.vector = tuple.to.list(state.vector)
        else state.vector = list(state.vector)
    }
    if (is.call(A) && A[[1]]=='TUPLE') A = tuple.to.list(A)
    stopifnot( length(A) == length(state.vector) && length(A[[1]]) == length(state.vector) ) 
    foo = guess.ode.time.variable(state.vector)
    re = list(state.vector = foo$names,
         time.variable = foo$time.variable,
         A = symbolic.listoflist.to.RMatrix(A))
    class(re)=c('firstoder.linear.ode','ODE')
    re
}

#' generate deSolve usable derivatives function
generate.func.firstorder.linear.ode = function(ode){
    L = list()
    st = lapply(1:phase.space.dim(ode), function(i) CONS('[',quote(y),as.numeric(i)))
    for(i in 1:length(ode$A)) {
        L[[i]] = simplify.pretty.01(symbolic.inner.product(ode$A[[i]], st))
    }
    re = function(t, y, param, ... ) {
    }
    body(re) = CONS('list',as.call(c(quote(c), L)))
    re
}

#' generate state vectors
#' the storage of ode , the names and time variable are separated
#' need to combine them together
#' @param ode ode system
#' @return the vector function as state function
state.vectors = function(ode) {
    lapply(ode$state.vector, function(x) CONS(x, ode$time.variable))
}

#' guess time variable 
#' @param vc \code{state.vector} functions list \code{ [ x1(t), x2(t) ] }
#' @return the names of functions and the time variable
guess.ode.time.variable  = function(vc){
    names = list()
    tv = NULL
    for(i in seq_along(vc)) {
        if (!symbolic.match(quote( '?s'(name1)('?'(time.var)) ), vc[[i]])) {
            stop(sprintf('The %s -th entry [%s] seems not a unknown function?', i, paste(deparse(vc[[i]]),collapse='')))
        }
        names[[i]]=name1
        if (is.null(tv)) {
            tv = time.var
            next()
        }
        if (! tv == time.var) {
            stop(sprintf('The %s -th time variable [%s] is equal to previous time variable %s . ',
                              i,paste(deparse(time.var),collapse=''), paste(deparse(tv),collapse='')))
        }
        if (i>1) {
            if (!is.na(ind<-match.symbolic(name1, names[-i]))) {
                stop(sprintf('The %s -th state vector is duplicated with previous %s -th vector. ', i, ind))
            }
        }
    }
    list(names = names, time.variable=tv)
}

#' the length of state.vector
#' the dimensional of state space (or say, phase space)
#' @param ode the ode equations
#' @return number
phase.space.dim = function(ode){
    return(length(ode$state.vector))
}

#' will put more and more know ode equations here
ode.librarys = list(
list(keywords=c('advan2','one compartment'),
        constructor=function() {
             create.firstorder.linear.ode( 
                state.vector = quote( TUPLE( DEPOT(t), CENTRAL(t) ) ),
                A = quote( TUPLE(
                                TUPLE(  -Ka, 0 ),
                                TUPLE(   Ka, -K)))) },
        compile.to.des.eqns = function(ode, Ka=Ka, K=K) {
                N = phase.space.dim(ode)
                A = ode$A
                for(i in 1:N){
                    for(j in 1:N){
                        A[[i]][[j]] = symbolic.gsub(quote(Ka),Ka,A[[i]][[j]])
                        A[[i]][[j]] = symbolic.gsub(quote(K),K,A[[i]][[j]])
                    }
                }
                st = lapply(1:N, function(x) CONS('[',quote(A),as.numeric(x)))
                eqns= vector(N, mode='list')
                for(i in 1:N){
                    eqns[[i]] = CONS('=', CONS('[',quote(DADT),as.numeric(i)), 
                       simplify.pretty.01( simplify(symbolic.inner.product( A[[i]], st)) ))
                }
                eqns })
)

#' change \code{A[[i]][[j]] = v}
#' @param ode ode system
#' @param i x index
#' @param j y index
#' @param v new value
#' @return new system
#' @examples
#'  # tmp = matrix(1:4,2)
#'  # set ode to zero
#'  # ode.update.ijv(ode, row(tmp), col(rmp), 0)
ode.update.ijv = function(ode,i,j,v) {
    if (!is.list(v)) v = list(v)
    NN = c( length(i), length(j), length(v) )
    N = max(NN)
    for(k in 1:N) {
        inds = (k-1) %% NN + 1 
        ode$A[[ i[inds[1]] ]][[ j[inds[2]] ]] = v[[ inds[3] ]]
    }
    ode
}

#' rendering an ode or state.vector of ode to nonmem model description
#' @param stvs state.vectors or ode
#' @return a string
state.vector.nonmem.models.block = function(stvs) {
    if (is(stvs,'ODE')) {
        stvs = stvs$state.vector
    }
    re = paste( sprintf('COMP=(%s)', sapply(stvs, as.character)), collapse='\n')
    indent.n = function(ss,n) {
       if (n==0) return(ss)
       indent = paste(rep(' ',n),collapse='')
       ss = gsub('^', indent, ss)
       gsub('\n', sprintf('\n%s',indent), ss)
    }
    re = indent.n(re,4)
    re = sprintf('$MODEL\n%s\n',re)
    re
}

#' by some criterion(TO DO)
#' determine which ODE solver be used in NONMEM
#' @param ode ode system
#' @return statement
ode.choose.nonmem.subroutine = function(ode) {
    stopifnot(is(ode,'ODE'))
    '$SUBR ADVAN6 TOL=4'
}

#' Compute the matching score for keywords in a list of keywords
#' @param keywords
#' @param list of keywords
score.keywords.match = function(keywords, keywords.list) {
    kws = keywords.list
    kw = unique(unlist(kws,rec=T))
    tb.ind = lapply(kws, function(x) match(x, kw))
    tb = matrix(0, nrow = length(kws), ncol=length(kw))
    for(i in 1:NROW(tb)) tb[i, tb.ind[[i]]] = 1
    attr.vec.ind = pmatch(keywords,kw)
    attr.vec.ind = attr.vec.ind[!is.na(attr.vec.ind)]
    attr.vec = matrix(0, nrow=NCOL(tb), ncol=1)
    attr.vec[ attr.vec.ind, ] = 1
    # dimension reduction if too many distinct keyword in keywords.list
    if (NCOL(tb)>5){
        foo = svd(tb)
        uu = foo$u[,1:5]
        vv = t(foo$v[,1:5])
        score.attr.vec = vv %*% attr.vec
        scores = uu %*% (foo$d[1:5] * score.attr.vec)
    } else {
        scores = tb %*% attr.vec
    }
    scores
}

#' simply find a library definition
#' @param keywords keywords for looking up
#' @return the entry
ode.apropos = function(keywords) {
    kws = lapply(ode.librarys, function(x) x$keywords)
    scores = score.keywords.match(keywords, kws)
    ode.librarys[ which.max(scores) ]
}

#' ode FORMATTERs
print.firstorder.linear.ode0 = function(ode){
    N = phase.space.dim(ode)
    st = state.vectors(ode)
    txt <- tlhs <- trhs <- character(0)
    for(i in 1:N){
        lhs1 = sprintf('d %s', paste(deparse(st[[i]]), collapse=''))
        lhs2 = paste(rep('-',nchar(lhs1)+2),collapse='')
        lhs3 = sprintf('d %s', paste(deparse(ode$time.variable),collapse=''))
        padding = round((nchar(lhs1) - nchar(lhs3))/2) + 1
        if (padding>0) lhs3 = sprintf('%s%s',paste(rep(' ',padding ),collapse=''), lhs3)
        lhs1 = sprintf(' %s', lhs1)
        rhs = paste(deparse(simplify.pretty.01(simplify(symbolic.inner.product(ode$A[[i]], st)))),collapse='')
        tlhs = c(tlhs, 
            lhs1,
            lhs2,
            lhs3
        )
        trhs = c(trhs, '', rhs,'')
    }
    txt = as.table(cbind(tlhs,c('','=',''),trhs))
    colnames(txt) = rep('',NCOL(txt))
    rownames(txt) = rep('',NROW(txt))
    txt = capture.output(print(txt))
    if (regexpr('^\\s*$',txt[1])>0) txt = txt[-1]
    paste(txt,collapse='\n')
}

#' S3 version
print.ODE = function(ode){
    if (is(ode,'firstoder.linear.ode')) cat(print.firstorder.linear.ode0(ode),'\n')
}
