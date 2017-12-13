########################################################################
# general (not nonmem specific) functions for symbolic operations      #
########################################################################

#' match.symbolic
#'
#' match for symbolic expressions, similar to R itself's \code{match} for characters
#' here \code{x} is just expression, not having pattern variables, the equality test
#' just comparing them using R's \code{==} equality test operator
#'
#' This is totally different from \code{symbolic.match}, should not be missly used
#'
#' @param x expression
#' @param table lookup table
#' @param na na symbol
#' @return index found
match.symbolic = function(x, table, na=NA){
    if (length(x)==0) return(integer(0))
    if (is.list(x)) {
        re = integer(length(x))
        for(i in 1:length(re)){
            re[i]=Recall(x[[i]], table, na)
        }
        return(re)
    }
    if (length(table)==0) return(na)
    ind = which(length(x) == sapply(table, length))
    if (length(ind)==0) return(na)
    ind = ind[ sapply(table[ind], function(y) y==x ) ]
    if (length(ind)==0) return(na)
    ind
}

#' union.symbolic
#'
#' union of symbolic expressions, similar to R's \code{union} for characters
#' 
#' @param le list of expressions
#' @param x  expression
#' @return list of unique expressons
union.symbolic = function(le=list(), x){
    if (is.list(x)) {
        for(i in x) {
            le = Recall(le,i)
        }
        return(le)
    }
    if (length(le)==0) {
        le[[1]] = x
        return(le)
    }
    found = F
    for(i in 1:length(le)){
        if (length(le[[i]])==length(x) && le[[i]]==x) {
            found=T
            break
        }
    }
    if (found) return(le)
    c(le,x)
}

#' extract.vertices
#'
#' From one equation \code{LHS = RHS}
#' extract several symbols in the right hand side \cr
#' and symbol in left hand side \cr
#' the array[i] array[j] will be extract just as array[]
#' 
#' Those symbols lhs and rhs are actually vertices in a relation graph.
#' If the expression is a simple equal relation, the graph are the complete bipartite graph, i.e.
#'    all pairs \code{ <rhs[i],lhs[j]> }, for any \code{i,j}.
#'
#' If the expression is a "if" statement, and the "yes" "no" part contain multiple statements,
#'    the edges are not the all pairs as above. Should use \code{create.abstractGraph.from.equations}
#'    to get the actual edges.
#'
#' @param e0 parsed R expression
#' @return list of lhs and rhs, which are just expressions
#' @author jjxie
extract.vertices = function(e0){
    lhs.symbol = list()
    append.lhs = function(x){
        if (length(lhs.symbol)==0) {
            lhs.symbol[[1]]<<-x
            return()
        }
        found=F
        for(i in 1:length(lhs.symbol)){
            if (length(lhs.symbol[[i]])==length(x) && lhs.symbol[[i]]==x) {
                found=T
                break
            }
        }
        if (!found) lhs.symbol[[length(lhs.symbol)+1]] <<- x
    }
    rhs.symbol = list()
    append.rhs = function(x){
        if (length(rhs.symbol)==0) {
            rhs.symbol[[1]]<<-x
            return()
        }
        found=F
        for(i in 1:length(rhs.symbol)){
            if (length(rhs.symbol[[i]])==length(x) && rhs.symbol[[i]]==x) {
                found=T
                break
            }
        }
        if (!found) rhs.symbol[[length(rhs.symbol)+1]] <<- x
    }
    walker = function(e){
        if (is.atomic(e)) return(e)
        if (is.symbol(e)) {
            append.rhs(e)
            return(e)
        }
        if (e[[1]]=='[') {
            # deal with array name specially
            stopifnot(is.symbol(e[[2]]))
            # of course, something very complex as
            #             THETA[ A[i,B[j]] ]
            # cannot work as expected!
            # in that case, extract.edges.quotient.array might be easier because it's dropping some information
            append.rhs(e)
            for(ii in 3:length(e)){
                e[[ii]] = Recall(e[[ii]])
            }
            return(e)
        }
        if (e[[1]]=='=' || e[[1]]=='<-') {
            if (is.atomic(e[[2]])){
                stop('invalid lhs of assignment')
            }
            if (is.symbol(e[[2]])){
                append.lhs(e[[2]])
            } else {
                # should be array, only the array's name is put into lhs
                # The indices if any is a variable, should be put into rhs
                e2 = e[[2]]
                stopifnot(is.language(e2))
                if (e2[[1]] != '[') {
                    stop('invalid left hand side! should only be simple variable or array')
                }
                # Array name must be SYMBOL
                stopifnot(is.symbol(e2[[2]]))
                append.lhs(e2)
                for(ii in 3:length(e2)){
                    e2[[ii]] = Recall(e2[[ii]])
                }
                e[[2]] = e2
            }
            e[[3]] = Recall(e[[3]])
            return(e)
        }
        if (e[[1]]=='if'){
            arity = length(e) - 1
            if(arity==2){
            # if cond {}
                e[[2]] = Recall(e[[2]])
                e[[3]] = Recall(e[[3]])
                return(e)
            } else {
            # if cond {} else {}
                e[[2]] = Recall(e[[2]])
                e[[3]] = Recall(e[[3]])
                e[[4]] = Recall(e[[4]])
                return(e)
            }
            return(e)
        }
        if (e[[1]]=='==' ||
            e[[1]]=='!=' ||
            e[[1]]=='>=' ||
            e[[1]]=='<=' ||
            e[[1]]=='>' ||
            e[[1]]=='<' ||
            e[[1]]=='&&' ||
            e[[1]]=='||'
            ) {
            # #
            e[[2]] = Recall(e[[2]])
            e[[3]] = Recall(e[[3]])
            return(e)
        }
        if (e[[1]]=='+' ||
            e[[1]]=='*' ||
            e[[1]]=='/' ||
            e[[1]]=='^' ){
            # #
            e[[2]] = Recall(e[[2]])
            e[[3]] = Recall(e[[3]])
            return(e)
        }
        if (e[[1]]=='-') {
            e[[2]] = Recall(e[[2]])
            if (length(e)==3) e[[3]] = Recall(e[[3]])
            return(e)
        }
        # a call may have not argument e.g. a()
        if (length(e)<=1) return(e)
        for(i in 2:length(e)){
        # we should not include function names as symbols
            e[[i]] = Recall(e[[i]])
        }
        e
    }
    walker(e0)
    list(lhs.symbol=lhs.symbol,rhs.symbol=rhs.symbol)
}

#' create.abstractGraph.from.equations
#'
#' Note: Each equation should be in a single assignment form as: \code{LHS = f(RHS1,RHS2,...)}
#' This is not correct for a general statement. e.g. if(..) { a1=f1(...);a2=f2(...) } else { a1=f1'(...); a2=f2'(...) } 
#' It's not correct to construct the full binary graph (i.e. connect each lhs and rhs vertex), a new code is added to handle this.
#'
#' And , duplicated edges will be generated, if you put two same lines there:
#' \code{   A = B
#'          A = B }
#' will give two same edges in graph, that from B pointing to A.
#'
#' 
#' similar to \code{\link{create.graph.from.equations}}
#' but in edge representation, the vertex are not characters, just numbers,
#' moreover, we record the vertices table of expressions
#' @param eqns equations expressions
create.abstractGraph.from.equations = function(eqns){
    vertex = list()
    edges = list()
    N.edges = 0
    backref = integer(0)
    eqns.vertex = vector(length(eqns),mode='list')
    eqns.lhs.rhs = vector(length(eqns),mode='list')
    lhs.ind <- rhs.ind <- integer(0)
    for(i in 1:length(eqns)){
        e0 = extract.vertices(eqns[[i]])
        # if rhs is empty, we assume it's some constant
        if (length(e0$rhs.symbol)==0) {
            if (length(e0$lhs.symbol)==0) {
                # ignore if not giving a edges
                next()
            }
        }
        vertex = union.symbolic(vertex, e0$rhs.symbol)
        vertex = union.symbolic(vertex, e0$lhs.symbol)
        lhs = match.symbolic(e0$lhs.symbol, vertex)
        rhs = match.symbolic(e0$rhs.symbol, vertex)
        eqns.vertex[[i]]=c(lhs,rhs)
        eqns.lhs.rhs[[i]]=list(lhs=lhs,rhs=rhs)
        lhs.ind = c(lhs.ind, lhs)
        rhs.ind = c(rhs.ind, rhs)
        # deal with if () {...} else {...} carefully
        if (is.call(eqns[[i]]) && eqns[[i]][[1]]=='if') {
            e1 = eqns[[i]]
            ## CONDITION
            foo = extract.vertices(e1[[2]])
            lhs = match.symbolic(foo$lhs.symbol, vertex)
            rhs = match.symbolic(foo$rhs.symbol, vertex)
            for(j in lhs){
                for(k in rhs){
                    edges[[ N.edges<-N.edges+1 ]] = c(k,j) 
                    backref = c(backref, i)
                }
            }
            ## TRUE and FALSE PART
            for(ary in 3:length(e1)){
                if (is.call(e1[[ary]]) && e1[[ary]][[1]]=='{') {
                    gg = Recall(e1[[ary]][-1])
                } else {
                    gg = Recall(e1[ ary ])
                }
                # merge it with current
                gg.ind = match.symbolic(gg$vertex.table, vertex)
                # add the sub-graph's edge
                for(sub.edge in gg$graph$edges){
                    edges[[ N.edges<-N.edges+1 ]] = gg.ind[ sub.edge ]
                    backref = c(backref, i)
                }
                # add extra edges from conditions rhs to sub.g's lhs
                for(sub.lhs in gg.ind[unique(gg$lhs.ind)]){
                    for(cond.rhs in rhs) {
                        edges[[ N.edges<-N.edges + 1 ]] = c(cond.rhs, sub.lhs)
                        backref = c(backref, i)
                    }
                }
            }
        } else {
        # simple A = f(b1,b2,b3) 's edges are full binary graph
            for(j in lhs){
                for(k in rhs){
                    edges[[ N.edges<-N.edges+1 ]] = c(k,j) 
                    backref = c(backref, i)
                }
            }
        }
    }
    list(vertex.table=vertex, 
        eqns.vertex=eqns.vertex,
        eqns.lhs.rhs=eqns.lhs.rhs,
        lhs.ind=lhs.ind,
        rhs.ind=rhs.ind,
        backref=backref,
        graph=list(vertex=seq_along(vertex),edges=edges))
}

#' Simplify simplify an expression
#'
#' Very elementary mathematical simplification of an R (symbolic) expression
#' \itemize{
#'  \item Arithmetic    \code{0*expression} \samp{->} \code{0}, \code{1*expression} \samp{->} \code{expression}, e.t.c
#'  \item  Logic        \code{expression && F} \samp{->} \code{F},  \code{expression || T} \samp{->} \code{T},   e.t.c        
#'  \item Other         \code{\{one expression\}} \samp{->} \code{one expression}, \code{(one expression)} \samp{->} \code{one expression}
#'  \item Atomic if both operands are atomic, then do the evaluation and replace them by the result
#' }
#' @param e R expression
#' @return simplified expression
#' @author jjxie
simplify = function(e){
    if (is.atomic(e)) return(e)
    if (is.symbol(e)) return(e)
    if (is.call(e) && length(e)==2 && e[[1]]=='-' && is.numeric(e[[2]])) {
    # quote( - 1 ) is CONS('-',1)
    # we should ensure it's atomic
        return( -e[[2]] )
    }
    if (e[[1]]=='=') {
        e[[2]] = Recall(e[[2]])
        e[[3]] = Recall(e[[3]])
        return(e)
    }
    if (e[[1]]=='if'){
        arity = length(e) - 1
        if(arity==2){
        # if cond {}
            e[[2]] = Recall(e[[2]])
            e[[3]] = Recall(e[[3]])
            if(is.atomic(e[[2]]) && as.logical(e[[2]])) {
                return(e[[3]])
            }
            return(e)
        } else {
        # if cond {} else {}
            e[[2]] = Recall(e[[2]])
            e[[3]] = Recall(e[[3]])
            if(is.atomic(e[[2]]) && as.logical(e[[2]])) {
                return(e[[3]])
            }
            e[[4]] = Recall(e[[4]])
            if(is.atomic(e[[2]]) && !as.logical(e[[2]])) {
                return(e[[4]])
            }
            return(e)
        }
        return(e)
    }
    if (e[[1]]=='ifelse' && length(e)==4) {
        e[[2]]=Recall(e[[2]])
        if (is.logical(e[[2]])) {
            if (is.na(e[[2]])) return(NA)
            if (e[[2]]) return( Recall(e[[3]]) )
            else return( Recall(e[[4]]))
        }
        e[[3]]=Recall(e[[3]])
        e[[4]]=Recall(e[[4]])
        return(e)
    }
    if (e[[1]]=='=='){
        e[[2]] = Recall(e[[2]])
        e[[3]] = Recall(e[[3]])
        if (is.atomic(e[[2]]) && is.atomic(e[[3]])){
            return( e[[2]] == e[[3]] )
        }
        return(e)
    }
    if (e[[1]]=='!='){
        e[[2]] = Recall(e[[2]])
        e[[3]] = Recall(e[[3]])
        if (is.atomic(e[[2]]) && is.atomic(e[[3]])){
            return( e[[2]] != e[[3]] )
        }
        return(e)
    }
    if (e[[1]]=='>='){
        e[[2]] = Recall(e[[2]])
        e[[3]] = Recall(e[[3]])
        if (is.atomic(e[[2]]) && is.atomic(e[[3]])){
            return( e[[2]] >= e[[3]] )
        }
        return(e)
    }
    if (e[[1]]=='<='){
        e[[2]] = Recall(e[[2]])
        e[[3]] = Recall(e[[3]])
        if (is.atomic(e[[2]]) && is.atomic(e[[3]])){
            return( e[[2]] <= e[[3]] )
        }
        return(e)
    }
    if (e[[1]]=='>'){
        e[[2]] = Recall(e[[2]])
        e[[3]] = Recall(e[[3]])
        if (is.atomic(e[[2]]) && is.atomic(e[[3]])){
            return( e[[2]] > e[[3]] )
        }
        return(e)
    }
    if (e[[1]]=='<'){
        e[[2]] = Recall(e[[2]])
        e[[3]] = Recall(e[[3]])
        if (is.atomic(e[[2]]) && is.atomic(e[[3]])){
            return( e[[2]] < e[[3]] )
        }
        return(e)
    }
    if (e[[1]]=='&&'){
        e[[2]] = Recall(e[[2]])
        if (is.atomic(e[[2]]) && !as.logical(e[[2]])) return(F)
        e[[3]] = Recall(e[[3]])
        if (is.atomic(e[[3]]) && !as.logical(e[[3]])) return(F)
        if (is.atomic(e[[2]]) && is.atomic(e[[3]])){
            return( e[[2]] && e[[3]] )
        }
        return(e)
    }
    if (e[[1]]=='||'){
        e[[2]] = Recall(e[[2]])
        if (is.atomic(e[[2]]) && as.logical(e[[2]])) return(T)
        e[[3]] = Recall(e[[3]])
        if (is.atomic(e[[3]]) && as.logical(e[[3]])) return(T)
        if (is.atomic(e[[2]]) && is.atomic(e[[3]])){
            return( e[[2]] || e[[3]] )
        }
        return(e)
    }
    if (e[[1]]=='+') {
            e[[2]] = Recall(e[[2]])
            e[[3]] = Recall(e[[3]])

            if (is.atomic(e[[2]]) && e[[2]]==0) {
                if (is.atomic(e[[3]]) && e[[3]]==0) {
                    return(0)
                } else {
                    return(e[[3]])
                }
            } else {
                if (is.atomic(e[[3]]) && e[[3]]==0) {
                    return(e[[2]])
                } 
            }
            if (  is.atomic(e[[2]]) && is.atomic(e[[3]]) ) {
                e = e[[2]] + e[[3]]
                return(e)
            }
            return(e)
    }
    if (e[[1]]=='-') {
            if (length(e)==2) {
                # unity '-'
                e[[2]] = Recall(e[[2]])
                if (is.atomic(e[[2]]) && e[[2]]==0) return(0)
                return(e)
            }
            e[[2]] = Recall(e[[2]])
            e[[3]] = Recall(e[[3]])
            if (is.atomic(e[[3]]) && e[[3]]==0) {
                return(e[[2]])
            } 
            if (  is.atomic(e[[2]]) && is.atomic(e[[3]]) ) {
                e = e[[2]] - e[[3]]
            }
            return(e)
    }
    if (e[[1]]=='*') {
        if (is.atomic(e[[2]]) && e[[2]] == 0) return(0)
        if (is.atomic(e[[3]]) && e[[3]] == 0) return(0)
        e[[2]] = Recall(e[[2]])
        e[[3]] = Recall(e[[3]])

        if (is.atomic(e[[2]]) && e[[2]]==1) {
            if (is.atomic(e[[3]]) && e[[3]]==1) {
                return(1)
            } else {
                return(e[[3]])
            }
        } else {
            if (is.atomic(e[[3]]) && e[[3]]==1) {
                return(e[[2]])
            } 
        }
        if (  is.atomic(e[[2]]) && is.atomic(e[[3]]) ) {
            e = e[[2]] * e[[3]]
        }
        return(e)
    }
    if (e[[1]]=='/') {
        e[[3]] = Recall(e[[3]])
        if (is.atomic(e[[3]]) && e[[3]] == 0) stop('Divided by zero!')
        e[[2]] = Recall(e[[2]])
        if (is.atomic(e[[2]]) && e[[2]] == 0) return(0)
        if (  is.atomic(e[[2]]) && is.atomic(e[[3]]) ) {
            e = e[[2]] / e[[3]]
        }
        return(e)
    }
    if (e[[1]]=='^') {
        e[[3]] = Recall(e[[3]])
        if (is.atomic(e[[3]]) && e[[3]] == 0) return(1)
        e[[2]] = Recall(e[[2]])
        if (is.atomic(e[[2]]) && e[[2]] == 0) return(0)
        if (is.atomic(e[[2]]) && e[[2]] == 1) return(1)
        if (is.atomic(e[[2]]) && is.atomic(e[[3]]) ) {
            e = e[[2]] ^ e[[3]]
            return(e)
        }
        if (is.atomic(e[[3]]) && is.numeric(e[[3]]) && e[[3]]==1) {
        # a^1 => a
            return(e[[2]])
        }
        return(e)
    }
    if (e[[1]]=='{') {
        N.e = length(e)
        if (N.e==2) {
            return(Recall(e[[2]]))
        }
        if (N.e>2) {
            for(ii in 2:N.e){
                e[[ii]] = Recall(e[[ii]])
            }
            return(e)
        }
        return(e)
    }
    if (e[[1]]=='('){
        e[[2]] = Recall(e[[2]])
        if (length(e)==2) return(e[[2]])
        return(e)
    }
    if (e[[1]]=='sqrt' && length(e)==2) {
        e2 = Recall(e[[2]])
        if (is.atomic(e2) && is.numeric(e2)) {
            tmp = round(sqrt(e2),4)
            if (abs(tmp * tmp - e[[2]])<1e-15) {
                return(tmp)
            }
        }
        e[[2]] = e2
        return(e)
    }
    if (e[[1]]=='exp' && length(e)==2) {
        e2 = Recall(e[[2]])
        if (is.atomic(e2) && e2==0) {
            return(as.numeric(1))
        }
        e[[2]] = e2
        return(e)
    }
    # should record logical expression as 
    # return as is , ignore unknown ones
    for(i in 1:length(e)){
        e[[i]] = Recall(e[[i]])
    }
    e
}

########################################################################
#  Symbolic matcher related                                            #
#  Rules operation related                                             #
########################################################################

#' QCONS : quoted CONS
#'
#'  Construct an R call, the first argument is evaluated, the others are quoted
#'
#'  \code{QCONS('*', 1, 2+3)} \samp{->} \code{1 * (2+3)}
#'
#'  this is useful , because 
#'
#'   e.g.   \code{QCONS(is.integer ,a )} \samp{->} \code{is.integer(a)}
#'
#'     where \code{a} is a variable
#'
#'    then we can use \code{eval( QCONS(is.integer, a) )} to the result
#'
#' @param x1 : the function, can be a variable, or quoted symbol
#' @param ... : other arguments
#' @return the call
#' @author jjxie
QCONS = function(x1,...){
    m = match.call()[-1]
    m[[1]]=if (is.character(x1)) as.symbol(x1) else x1
    m
}

#' CONS : evaluate argument other than first one
#'
#' similar as QCONS, however, will evaluated the argument first, as following
#'
#'  \code{CONS('*', 1, 2+3)} \samp{->} \code{1 * 5}
#'
#' @param x1 : the function, can be a variable, or quoted symbol
#' @param ... : other arguments
#' @return the call
#' @author jjxie
#' @seealso \code{\link{QCONS}}
CONS = function(x1,...){
    m = match.call()[-1]
    m[[1]]=if (is.character(x1)) as.symbol(x1) else x1
    l = list(...)
    if (length(l)>0) {
        for(i in 1:length(l)) {
            m[[ i+1 ]] = l[[i]]
        }
    }
    m
}

#' avoid the warnings if length is different
#' @param e1 symbol or call
#' @param e2 symbol or call
#' @return logical
symbolic.check.raw.equal = function(e1,e2){
     length(e1)==length(e2) && e1 == e2
}

#' symbolic.match0
#'
#' working horse for symbolic matcher, match a pattern expression against a given expression
#'
#' 
#' A pattern is a expression might have one or more pattern variables,\cr
#' pattern variable are as \code{'?'(x)}, \code{'?'(y)} , e.t.c\cr
#' note \code{'??'(x)} is equivlent to \code{'?'(x)}\cr
#' also \code{'?a'(x)} can test an atomic term and \code{'?s'(x)} can be used to test a symbol term 
#'
#' Question mark with type-tester, e.g. \cr
#'        \code{'?'(is.atomic,A)} to match \code{A} if \code{A} is atomic, \cr
#'        \code{'?'(is.symbol,A)} to match \code{A} if \code{A} is symbol, \cr
#'        \code{'?'(is.numeric,A)} to match \code{A} if \code{A} is numeric, \cr
#'        \code{'?'(is.integer,A)} to match \code{A} if \code{A} is integer, \cr
#'        \code{'?'(is.character,A)} to match \code{A} if \code{A} is character, \cr
#'        \code{'?'(is.call,A)} to match \code{A} if \code{A} is call 
#' 
#' @param pat : pattern expression
#' @param e0 : expression to be matched against
#' @return an environment if matched, containing the dictionary of matched variables, or \code{NULL} if not matched
#' @examples
#'   symbolic.match0( quote( '?'(A) * exp(ETA[ '?'(B) ])), quote(  TVKA  * exp(ETA[ 2 ] ) ))
#'   # will return : env (dictionary, { A: TVKA, B:2 } )
#' @author jjxie
#' @seealso \code{\link{symbolic.match}}
symbolic.match0 = function(pat, e0) {
    pat.vars = new.env(parent=emptyenv())
    walker = function(pat,e){
        # special case : single entry
        if (is.atomic(pat)) {
        # Need to carefully dealwith both are NA, R will return NA==NA as NA
            if (is.atomic(e)) {
                na.pat = is.na(pat)
                na.e = is.na(e)
                return( (na.pat && na.e) || (!na.pat && !na.e && pat==e) )
            }
            return(F)
        }
        if (is.symbol(pat)) {
            if (is.symbol(e) && pat==e) return(T)
            return(F)
        }
        # process with na, '?na'() to match a na
        if (pat[[1]]=='?na' && length(pat)==1) {
            return(is.atomic(e0) && is.na(e0))
        }
        # pattern case 1: pattern variable without type-tester
        if ((pat[[1]]=='?'||pat[[1]]=='??') && length(pat)==2) {
            # pattern.variable without type
            stopifnot(is.symbol(pat[[2]]))
            varname = as.character(pat[[2]])
            if (exists(varname, env=pat.vars)) {
                return( symbolic.check.raw.equal(get(varname,env=pat.vars) , e) )
            } else {
                assign(varname, e, env=pat.vars)
                return(T)
            }
        }
        # pattern case 2: pattern variable with type-tester
        if (pat[[1]]=='?' && length(pat)==3 && 
            eval(QCONS(pat[[2]], e)) ) {
            # pattern.variable with type
            stopifnot(is.symbol(pat[[3]]))
            varname = as.character(pat[[3]])
            if (exists(varname, env=pat.vars)) {
                return( get(varname,env=pat.vars) == e )
            } else {
                assign(varname, e, env=pat.vars)
                return(T)
            }
        }
        # Just for handy use
        # '?a'(X) == '?'(is.atomic,X)
        if (pat[[1]]=='?a' && length(pat)==2 &&
            is.atomic(e) ) {
            # pattern.variable without type
            stopifnot(is.symbol(pat[[2]]))
            varname = as.character(pat[[2]])
            if (exists(varname, env=pat.vars)) {
                return( get(varname,env=pat.vars) == e )
            } else {
                assign(varname, e, env=pat.vars)
                return(T)
            }
        }
        # Just for handy use
        # '?s'(X) == '?'(is.symbol,X)
        # '?v' is same as '?s'
        # however, in another function symbolic.simplify.gigo, we treat it differently
        if ((pat[[1]]=='?s'||pat[[1]]=='?v') && length(pat)==2 &&
            is.symbol(e) ) {
            # pattern.variable without type
            stopifnot(is.symbol(pat[[2]]))
            varname = as.character(pat[[2]])
            if (exists(varname, env=pat.vars)) {
                return( get(varname,env=pat.vars) == e )
            } else {
                assign(varname, e, env=pat.vars)
                return(T)
            }
        }
        # Just for handy use
        # '?c'(X) == '?'(is.call,X)
        if (pat[[1]]=='?c' && length(pat)==2 &&
            is.call(e) ) {
            # pattern.variable without type
            stopifnot(is.symbol(pat[[2]]))
            varname = as.character(pat[[2]])
            if (exists(varname, env=pat.vars)) {
                return( get(varname,env=pat.vars) == e )
            } else {
                assign(varname, e, env=pat.vars)
                return(T)
            }
        }
        # '?()' : fun of call
        # match when it's a function part of a call
        if (pat[[1]]=='?()' && length(pat)==2 &&
            is.call(e) ) {
            # pattern.variable without type
            stopifnot(is.symbol(pat[[2]]))
            varname = as.character(pat[[2]])
            if (exists(varname, env=pat.vars)) {
                return( get(varname,env=pat.vars) == e[[1]] )
            } else {
                assign(varname, e[[1]], env=pat.vars)
                return(T)
            }
        }
    	# '?[]' : name of array
        # match when it's a name of an array
        if (pat[[1]]=='?[]' && length(pat)==2 &&
            is.call(e) && e[[1]]=='[' ) {
            # pattern.variable without type
            stopifnot(is.symbol(pat[[2]]))
            varname = as.character(pat[[2]])
            if (exists(varname, env=pat.vars)) {
                return( get(varname,env=pat.vars) == e[[2]] )
            } else {
                assign(varname, e[[2]], env=pat.vars)
                return(T)
            }
        }
        # if pat[[1]] is not pattern variable, but the respect e is atomic or symbol, just return FALSE
        if (is.atomic(e) || is.symbol(e)) return(F)
        # general case : do a for loop
        if (length(pat)!=length(e)) return(F)
        for(i in 1:length(pat)){
            if (!Recall(pat[[i]], e[[i]])) return(F)
        }
        T
    }
    re = walker(pat,e0)
    if (re) return(pat.vars)
    NULL
}

#' parse1
#' parse oneline text to a expression (call)
#' @param s oneline string
#' @return an R call
#' @author jjxie
parse1 = function(s){
    parse(text=s)[[1]]
}

#' symbolic.match
#'
#' A more handy interface for \code{\link{symbolic.match0}}. \cr
#' This funciton's return value is only \code{T/F}, but it has side effect on parent environment
#'
#' e.g. for the example in symbolic.match0, it assigns back to the caller's environment  \cr
#'        \code{A = TVKA} \cr
#'        \code{B = 2} 
#'
#' @param pat pattern expression
#' @param e0 expression to be matched against
#' @return \code{TRUE} or \code{FALSE}
#' @author jjxie
symbolic.match = function(pat,e0){
    pat.vars = symbolic.match0(pat,e0)
    if (is.null(pat.vars)) return(FALSE)
    nms = ls(env=pat.vars)
    for(i in nms) {
        assign(i, get(i,env=pat.vars), pos=parent.frame(1))
    }
    TRUE
}

#' symbolic.instantiate 
#'
#' Given a dictionary, instantiate a new expression which might have instance variables
#' The \code{':'(x)} is a example of instance variable.
#' If \code{x} equals to \code{2}, then the following expression
#' \code{sin(':'(x))} will instantiate to \code{sin(2)}
#'
#' @param skeleton expression skeleton to be instantiated
#' @param dict dictionary with definition of instance variables
#' @return the instantiated expression
#' @author jjxie
#' @seealso \code{\link{symbolic.simplify.gigo}}
symbolic.instantiate = function(skeleton, dict){
    stopifnot(is.environment(dict))
    walker = function(e){
        if( is.atomic(e) || is.symbol(e) ) return(e)
        # ':'(x) 
        # arity is necessary, becaue R's : has arity 2, and here we only arity 1!
        if( e[[1]]==':' && length(e)==2 && is.symbol(nm<-e[[2]]) &&
            exists(as.character(nm), dict, inherits=FALSE)) {
            val = get(as.character(nm),dict,inherits=FALSE)
            return(val)
        }
        # e is a call
        for(i in 1:length(e)){
            e[[i]] = Recall(e[[i]])
        }
        e
    }
    walker(skeleton)
}

#' simplify.2
#' math level simplifier \cr
#'
#'      \code{log(A*B)} \samp{->} \code{log(A)+log(B)} \cr
#'      \code{log(A^B)} \samp{->} \code{ B * log(A) } \cr
#'      \code{exp(A)^B} \samp{->} \code{ exp(A*B) } \cr
#'      \code{exp(log(A))} \samp{->} \code{A} \cr
#'      \code{log(exp(A))} \samp{->} \code{A} \cr
#'      evaluate \code{log(atomic)} to a number
#' @param e0 R expression
#' @return simplified expression
#' @author jjxie
#' @seealso \code{\link{simplify}}
simplify.2 = function(e0){
    e0 = simplify.1(e0)
    walker = function(e) {
        if (is.atomic(e)) return(e)
        if (is.symbol(e)) return(e)
        if (length(e)==1) return(e)
        if (symbolic.match(quote(log(exp('?'(anything)))), e)) {
            return( Recall(anything) )
        }
        if (symbolic.match(quote(exp(log('?'(anything)))), e)) {
            return( Recall(anything) )
        }
#        if (e[[1]]=='log' && is.atomic(e[[2]])) {
#            return(log(e[[2]]))
#        }
#        if (e[[1]]=='exp' && is.atomic(e[[2]])) {
#            return(exp(e[[2]]))
#        }
        if (symbolic.match(quote(log( '?'(A) * '?'(B) )), e)) {
            tmp1 = CONS('log',A)
            tmp2 = CONS('log',B)
            tmp3 = CONS('+', tmp1, tmp2)
            return( Recall(tmp3) )
        }
        if (symbolic.match(quote(log( '?'(A) ^ '?'(B) )), e)) {
            tmp1 = CONS('log',A)
            tmp2 = CONS('*', B, tmp1)
            return( Recall(tmp2) )
        }
        if (symbolic.match(quote( exp('?'(A)) ^ '?'(B)), e)) {
            tmp1 = CONS('*', A, B)
            tmp2 = CONS('exp',tmp1)
            return( Recall(tmp2) )
        }
        if (symbolic.match(CONS('^',quote( '?'(A) ^ '?'(B)), quote('?'(C))), e)) {
            tmp1 = CONS('*', B, C)
            tmp2 = CONS('^',A,tmp1)
            return( Recall(tmp2) )
        }
        for(i in 1:length(e)){
            e[[i]]=Recall(e[[i]])
        }
        return(simplify.1(e))
    }
    simplify.1(walker(e0))
}

#################### Other Simplifier ########################################################

#' simplify.sqrt
#'
#' simplify expression with \code{sqrt}
#'
#' @param e expression
#' @return simplified expression
#' @author jjxie
simplify.sqrt = function(e){
    rules = list(
    list(quote( sqrt( ( '?'(X) ) )          ), quote( sqrt(':'(X)) ) ),
    list(quote( sqrt( '?'(X) ^ 2 )           ), quote( ':'(X) ) ),
    list(quote( sqrt( '?'(X) ^ '?a'(Y)     ) ), 
        function(e1,dict){ 
            y = floor(dict$Y / 2)
            if ( y * 2 == dict$Y ) {
                if (y==1) {
                    return(dict$X)
                }
                return( CONS('^', dict$X, y ) )
            }
            return(e1)
        }),
    list(quote( sqrt(exp('?'(X))) )          ,  quote( exp( 1/2 * ':'(X) ) ) ),
    list(quote( log(sqrt('?'(X))))           ,  quote( 1/2 * ':'(X) ) ),
    list(quote( sqrt( '?'(X) ) )             ,  quote( ':'(X) ^ (1/2) ) )
    )
    walker = symbolic.simplify.gigo(rules)
    simplify.2(walker(e))
}

##############################################################################################

#' symbolic.grep
#'
#'   recursively find a pattern in an expression
#'
#' @param pat pattern expression
#' @param e expression from which to find the pattern
#' @return \code{TRUE} or \code{FALSE} indicating if the pattern is found
#' @author jjxie
symbolic.grep = function(pat, e){
    if (symbolic.match(pat,e)) return(T)
    # if not a call, and not matched by above, return false
    if (!is.call(e)) return(F)
    if (length(e)==1) return(F)
    for(i in 2:length(e)){
        if (Recall(pat,e[[i]])) return(T)
    }
    F
}

#' symbolic.sub
#'
#'   recursively replace a pattern \code{pat} in an expression \code{e0} by \code{rep} expression \cr
#'   Only the first occurance is replaced, \cr
#'   If the pattern is not found, \code{e0} is return without changing
#'
#' @param pat pattern expression
#' @param rep expression as a replacement
#' @param e0 expression to be operated
#' @return the replaced expression
#' @author jjxie
symbolic.sub = function(pat, rep, e0){
    subbed = F
    walker = function(e) {
        if (symbolic.match(pat,e)) {
            subbed <<- T
            return(rep)
        }
        # if not a call, and not matched by above, return false
        if (!is.call(e)) return(e)
        # call may also have only 1 argument, e.g. a()
        for(i in 1:length(e)){
            e[[i]] = Recall(e[[i]])
            if (subbed) return(e)
        }
        e
    }
    walker(e0)
}

#' symbolic.gsub
#'
#'   recursively replace a pattern \code{pat} in an expression \code{e0} by an expression \code{rep} \cr
#'   similar to \code{\link{symbolic.sub}} with exception that all the occurance of \code{pat} are replaced
#'
#' @param pat pattern expression
#' @param rep expression as a replacement
#' @param e0 expression to be operated
#' @return the replaced expression
#' @author jjxie
symbolic.gsub = function(pat, rep, e0){
    walker = function(e) {
        if (symbolic.match(pat,e)) return(rep)
        # if not a call, and not matched by above, return false
        if (!is.call(e)) return(e)
        for(i in 1:length(e)){
            e[[i]] = Recall(e[[i]])
        }
        e
    }
    walker(e0)
}

#' symbolic.symbol.extract
#'
#' extract only symbols, excluding function call name and array names
#'
#' @param e0 expression to be looked at
#' @return a character vector containing unique symbols appearing in the expression
#' @author jjxie
symbolic.symbol.extract = function(e0){
    sym = character(0)
    walker = function(e){
        if (is.atomic(e)) return(NULL)
        if (is.symbol(e)) { sym <<- union(sym, as.character(e));return(NULL)}
        if (symbolic.match( quote('?[]'(x)), e)) { return(NULL) } 
        if (length(e)>1) {
            for(i in 2:length(e)) Recall(e[[i]])
        }
        NULL
    }
    walker(e0)
    sym
}

#' character.rules.to.list.rules
#'
#' a syntaxtic surger for translate character written rules into list of list form
#'
#'  Character vector containing rules should look like following: \cr
#'  \code{
#'      c( lhs1, rhs1,
#'         lhs2, rhs2, 
#'         ... , ... ,
#'       ) } \cr
#' And those strings are parsed and paired into list of pair of expressions
#'
#' @param s0 a character vector
#' @return list of pairs of expressions
#' @author jjxie
#' @seealso \code{\link{symbolic.simplify.gigo}}
character.rules.to.list.rules = function(s0){
    apply(matrix(s0,nrow=2),2, function(x) list(parse1(x[1]),parse1(x[2])))
}

#' symbolic.simplify.gigo
#'
#' Garbage-in-Garbage-out symbolic simplifier factory \cr
#' It takes in rules, and according the rules, generated an simplifier, which can operate on a given expression, \cr
#' doing transformations described by rules, and continue on trying the rules on the transformed expression, \cr
#' until no rules applied. \cr
#'
#' rules should be list of pairs as following:
#' \verb{
#'                  [ [ lhs1, rhs1 ],
#'	                  [ lhs2, rhs2 ],
#'                      ...,
#'                    ] }
#' Here, lhs should be a pattern with pattern expression \cr
#' rhs has two possibilities: 
#' \enumerate{
#'         \item  skeleton with \code{':'(x)}
#'         \item  a callback function takes in \code{e} (mathed expression) and dictionary (matched variable)
#'              like \code{function(e,dict) { do something on using e or dict} }
#'  }
#'
#' @param rules list of pairs of which defines transforming rules
#' @return a closure (R function) for doing the transformation (simplification) accordings to the given rules
#' @author jjxie
#' @seealso \code{\link{simplify.03}}, \code{\link{character.rules.to.list.rules}}
symbolic.simplify.gigo = function(rules){

    simplify.exp = function(e) {
        if (is.atomic(e) || is.symbol(e)) return(try.rules(e))
        # when e is a call
        # may need to deal with e[[1]] specially
        # purpose is to treat symbols of function name differently
        for(i in 1:length(e)){
            if (i==1 && is.symbol(e[[1]])) {
                if (length(rules)==0) next
                for(j in 1:length(rules)){
                    if (is.call(rules[[j]][[1]]) && rules[[j]][[1]][[1]]=='?v') {
                        # '?v'(x) wont match function call's name
                        next
                    } else {
                        dict = symbolic.match0(rules[[j]][[1]], e[[1]])
                        if (is.null(dict)) next
                        if (is.function(rules[[j]][[2]])) {
                            tmp = rules[[j]][[2]](e[[1]], dict)
                        } else {
                            tmp = symbolic.instantiate(rules[[j]][[2]], dict)
                        }
                        if (symbolic.match(tmp,e[[1]])) next
                        e[[1]] = Recall(tmp)
                        break
                    }
                }
                next
            }
            e[[i]] = Recall(e[[i]])
        }
        try.rules( e )
    }

    try.rules = function(e) {
        if(length(rules)==0) return(e)
        for(i in 1:length(rules)){
            dict = symbolic.match0(rules[[i]][[1]], e)
            if (is.null(dict)) next
            if (is.function(rules[[i]][[2]])) {
                tmp = rules[[i]][[2]](e, dict)
            } else {
                tmp = symbolic.instantiate(rules[[i]][[2]], dict)
            }
            # avoid infinite loop
            if (symbolic.match(tmp,e)) next
            return(simplify.exp(tmp))
        }
	    e
    }
    simplify.exp
}
###################################################################################################################
# polynomial related
###################################################################################################################

#' sort.by.object.size
#'
#' Internal use \cr
#' Sort a list of expressions by first, the size; then, deparsed result \cr
#'
#' @return a list of terms
#' @author jjxie
sort.by.object.size = function(le){
    sz = sapply(le, object.size)
    nm = sapply(le, function(x) deparse(x)[1])
    ord = order(sz,nm)
    le[ord]
}

#' expand.as.sum
#' 
#' Expand the summation into list, \cr 
#' \code{A + B + C + D} \samp{->} \code{[ A, B, C ,D]}
#'
#' @param e expression
#' @return list of expression
#' @author jjxie
expand.as.sum = function(e) {
    e = simplify(e)
    if (is.atomic(e) || is.symbol(e)) return(list(e))
    if (e[[1]]=='+') {
        return(as.list(c(Recall(e[[2]]),Recall(e[[3]]))))
    }
    if (e[[1]]=='-') {
        if (length(e)==2) {
            if (is.call(e[[2]]) && e[[2]][[1]]=='-' && length(e[[2]]==2)) {
                return(list(Recall(e[[2]])))
            }
            #
            return(list(e))
        }
        if (is.call(e[[3]]) && e[[3]][[1]]=='-' && length(e[[3]]==2)) {
            return(as.list(c(Recall(e[[2]], Recall(e[[3]][[2]])))))
        }
        return(as.list(c(Recall(e[[2]]), as.call(c(as.symbol('-'),Recall(e[[3]]))))))
    }
    # recursive, hence must be a list
    return( list( e ) )
}

#' simplify.03
#'
#' helper simplifier for following purpose \cr
#'  \code{(a^i)^j} \code{=>} \code{a^(i*j)} \cr
#' This simplifier is generated by \code{symbolic.simplify.gigo} factory, hence, also an example for using \cr
#'  \code{\link{symbolic.simplify.gigo}}
#'
#' @param e expression
#' @return simplified expression
#' @author jjxie
simplify.03 = symbolic.simplify.gigo(
    list( list( CONS('^', quote('?'(x) ^ '?'(is.numeric,a)), quote('?'(is.numeric,b)) ) , 
    function(e,d) { CONS('^',d$x, d$a*d$b) } )) )

#' expand.as.product
#' 
#' Expand the production into list, \code{A * B * C * D} \code{->} \code{[ A, B, C ,D]} \cr
#' also some extra translation is applied:   \cr
#'                                          \code{-A -> [-1,A] } \cr
#'                                          \code{A/B -> [A, B^-1] } \cr
#'                                          \code{A^1 -> [A] , (A^i)^j -> [A^i*j]} 
#' @param e expression
#' @return list of expression
#' @author jjxie
expand.as.product = function(e) {
    e = simplify(e)
    if (is.atomic(e) || is.symbol(e)) return(list(e))
    if (e[[1]]=='*') {
        return(as.list(c(Recall(e[[2]]),Recall(e[[3]]))))
    }
    if (e[[1]]=='-' && length(e)==2) {
        return(c( list(-1), Recall(e[[2]])))
    }
    if (e[[1]]=='/') {
        tmp = Recall(e[[3]])
        tmp = lapply(tmp, function(x) CONS('^',x,-1))
        return(as.list(c(Recall(e[[2]]),tmp)))
    }
    if (e[[1]]=='^' && length(e)==3) {
        tmp = Recall(e[[2]])
        if (is.atomic(e[[3]]) && e[[3]]==1) return(tmp)
        # tmp is now a list [x1,x2,...], should simplify each x_i^e[[3]]
        tmp = lapply(tmp, function(x) simplify.03(CONS('^',x,e[[3]])))
        return(tmp)
    }
    # recursive, hence must be a list
    return( list( e ) )
}

#' collect.atomic.product.list
#'
#' compute and combine atomics in a list of expression representing a production\cr
#' \code{[ A, B, 3, 4, -1 ] -> [ -12, A, B ]}
#' 
#' @param le list of expression
#' @return list of expression which atomics are combined
#' @author jjxie
collect.atomic.product.list = function(le){
    if(is.null(le) || length(le)==0) return(list(1))
    N = length(le)
    ind = logical(N)
    for(i in 1:N){
        if (is.atomic(le[[i]])) ind[i]=T
    }
    if (all(ind)) {
        return(prod(unlist(le[ind])))
    }
    if (any(ind)) {
        tmp = prod(unlist(le[ind]))
        if (tmp==0) return(list(0))
        if (tmp==1) return(le[!ind])
        return(c(tmp, le[!ind]))
    }
    le
}

#' simplify.product.list
#'
#' Do following light checking to remove redundant terms in list of production \cr
#' remove 1's from list, if any entry equals to 0 , return 0 directly \cr
#' if \code{A} and \code{A^(-1)} both exists, remove both.
#'
#' @param le  list of expression
#' @return  list of expression
#' @author jjxie
simplify.product.list = function(le){
    if(is.null(le) || length(le)==0) return(list(1))
    le = collect.atomic.product.list(le)
    N = length(le)
    if (N==1) return(le)
    need.remove = rep(F, N)
    # any 0 terms : return 0 
    # any 1       : remove it
    for(i in 1:N){
        if (is.atomic(le[[i]])) {
             if(le[[i]]==0) return( list(0) ) 
             if(le[[i]]==1) need.remove[ i ] = T
        }
    }
    # remove cancelled terms
    for(i in 1:N){
        if (need.remove[i]) next
        if (!(is.call(le[[i]]) && le[[i]][[1]]=='^' && le[[i]][[3]]==-1)) next
        for(j in 1:N){
            if (need.remove[j]) next
            if (i==j) next
            if (equal.symbol.light(le[[j]],le[[i]][[2]])) {
                need.remove[c(i,j)]=T
                break
            }
        }
    }
    if (all(need.remove)) return(list(1))
    le[!need.remove]
}

#' collect.atomic.sum.list
#'
#' compute and combine list of expressions representing a summation, similar to \code{\link{collect.atomic.product.list}} \cr
#' \code{[ A, B, 3, 4, -1 ] -> [ 6, A, B ]}
#' 
#' @param le list of expression
#' @return list of expression which atomics are combined
#' @author jjxie
collect.atomic.sum.list = function(le){
    if(is.null(le) || length(le)==0) return(list(0))
    N = length(le)
    ind = logical(N)
    for(i in 1:N){
        if (is.atomic(le[[i]])) ind[i]=T
    }
    if (all(ind)) {
        return(sum(unlist(le[ind])))
    }
    if (any(ind)) {
        tmp = sum(unlist(le[ind]))
        if (tmp==0) return(le[!ind])
        else return(c(tmp, le[!ind]))
    }
    le
}

#' simplify.sum.list
#' 
#' Similar to \code{\link{simplify.product.list}}, do following light checking to remove redundant terms in list of summation \cr
#' remove 0's from list \cr
#' if \code{A} and \code{-A} both exists, remove both
#'
#' @param le  list of expression
#' @return simplified list of expression 
#' @author jjxie
simplify.sum.list = function(le){
    if(is.null(le) || length(le)==0) return(list(0))
    le = collect.atomic.sum.list(le)
    N = length(le)
    if (N==1) return(le)
    need.remove = rep(F, N)
    # any 0 terms : return 0 
    # any 1       : remove it
    for(i in 1:N){
        if (is.atomic(le[[i]])) {
             if(le[[i]]==0) need.remove[ i ] = T
        }
    }
    # remove cancelled terms
    for(i in 1:N){
        if (need.remove[i]) next
        if (!(is.call(le[[i]]) && le[[i]][[1]]=='-' && length(le[[i]])==2)) next
        for(j in 1:N){ # by writing this, saves writing of swapped comparation
            if (need.remove[j]) next
            if (i == j) next
            if (equal.symbol.light(le[[j]],le[[i]][[2]])) {
                need.remove[c(i,j)]=T
                break
            }
        }
    }
    if (all(need.remove)) return(list(0))
    le[!need.remove]
}

#' collect.algebraical
#' 
#' For a list of expression, match them to coeff*indeterminators, and then combine together the indeterminators\cr
#' \code{ [ A, B*C, 5*A ] -> [6*A, B*C] } \cr
#' Here coefficients are atomic by \code{'?a'}
#'
#' @param le list of expressions
#' @return list of combined expressions
#' @author jjxie
collect.algebraical = function(le){
    if (length(le)==1) return(le)
    coeffs = integer(length(le))
    if (symbolic.match(quote('?a'(coeff) * '?'(terms)), le[[1]])) {
        mono = list( terms )
        coeffs[1] = coeff
    } else {
        mono = list(le[[1]])
        coeffs[1] = 1
    }
    N.mono = 1
    for(i in 2:length(le)){
        if (!symbolic.match(quote('?a'(coeff) * '?'(terms)), le[[i]])) {
            terms = le[[i]]
            coeff = 1
        }
        found = F
        for(j in 1:N.mono){
            if(equal.symbol.light(terms, mono[[j]])) {
                coeffs[j] = coeffs[j] + coeff
                found = T
                break
            }
        }
        if (!found){
            mono[[ (N.mono<-N.mono+1) ]] = terms
            coeffs[ N.mono ] = coeff
        }
    }
    sapply(1:N.mono, function(i) CONS('*',coeffs[i], mono[[i]]))
}

#' expression.from.list.of.sum
#'
#' Convert back list of sum to R expression\cr
#' inverse of \code{\link{expand.as.sum}}, from a list of expression, create a single expression
#'
#' @param le list of expression
#' @return one expression
#' @author jjxie
expression.from.list.of.sum=function(le){
    le = simplify.sum.list(collect.algebraical(le))
    if (length(le)==0) return(0)
    if (length(le)==1) return(le[[1]])
    tmp = le[[1]]
    for(i in 2:length(le)){
        tmp = as.call(c(as.symbol('+'),tmp,le[[i]]))
    }
    simplify(tmp)
}

#' expression.from.list.of.sum
#'
#' Convert back list of product to R expression \cr
#' inverse of \code{\link{expand.as.product}}, from a list of expression, create a single expression
#'
#' @param le list of expression
#' @return one expression
#' @author jjxie
expression.from.list.of.product=function(le){
    le = simplify.product.list(le)
    if (length(le)==0) return(1)
    if (length(le)==1) return(le[[1]])
    tmp = le[[1]]
    for(i in 2:length(le)){
        tmp = as.call(c(as.symbol('*'),tmp,le[[i]]))
    }
    simplify(tmp)
}

#' expression.from.sum.product
#'
#' Convert back list of sum of list of product to R expression \cr
#' e.g. \code{[ [a1,a2,a3] , [b1,b2,b3], [c1] ] -> a1*a2*a3 + b1*b2*b3 + c1}
#'
#' @param le list of expressions
#' @return one expression
#' @author jjxie
#' @seealso \code{\link{expand.as.sum.of.product}}
expression.from.sum.product=function(le){
    simplify(expression.from.list.of.sum(lapply(le, expression.from.list.of.product)))
}

#' simplify.1
#'
#' Simplify level 1: one level higher than simplify, by expand an expression to list of list \cr
#' removing redundant terms, and combine back \cr
#' e.g.    \code{A*B + C - A*B -> C}
#' 
#' @param e expression
#' @return simplified expression
#' @author jjxie
simplify.1 = function(e){
    expression.from.sum.product( expand.as.sum.of.product( e ))
}

#' expand.as.sum.of.product
#'
#' expand an expression into list of list as following \cr
#'   \code{ A*B + C*D + E*F -> [ [ A, B] , [C, D] , [E, F] ] }
#' 
#' non-deep-expand
#'  can do this : \code{(a+b)*c -> [ [a+b, c] ]}
#'  cannot do this : \code{(a+b)*c -> a*c + b*c}
#'
#' @param e expression
#' @return list of expressions
#' @author jjxie
#' @seealso \code{\link{print.listoflist}}
expand.as.sum.of.product = function(e){
    e1 = simplify.sum.list(expand.as.sum(e))
    e2 = lapply(e1, function(x) expression.from.list.of.product(expand.as.product(x)))
    lapply(e2, expand.as.product)
}

#' equal.symbol.light
#'
#' light weighted testing of equality \cr
#' Can ensure if the return value is \code{TRUE}, \code{e1} must equal to \code{e2};  \cr
#' however, it's not able to ensure when \code{e1} equals to \code{e2} mathematically, the return value is \code{TRUE}
#' 
#' @param e1 expression
#' @param e2 expression
#' @return \code{TRUE} or \code{FALSE}
#' @author jjxie
equal.symbol.light = function(e1,e2){
    if (length(e1)==length(e2) && e1==e2) return(T)

    if (is.atomic(e1)) {
        if (!is.atomic(e2)) return(FALSE)
        return( e1 == e2 )
    }
    if (is.atomic(e2)) {
        if (!is.atomic(e1)) return(FALSE)
        return( e1 == e2 )
    }
    if (is.symbol(e1)) {
        if (!is.symbol(e2)) return(FALSE)
        return( e1 == e2 )
    }
    if (is.symbol(e2)) {
        if (!is.symbol(e1)) return(FALSE)
        return( e1 == e2 )
    }
    e11 = sort.by.object.size(expand.as.sum(e1))
    e22 = sort.by.object.size(expand.as.sum(e2))
    if (length(e11)>1 && length(e11)==length(e22)){
        all.eq = T
        for(i in 1:length(e11)){
            if (!Recall(e11[[i]] , e22[[i]])) {
                all.eq=F
                break
            }
        }
        if (all.eq) return(T)
    }
    e11 = sort.by.object.size(expand.as.product(e1))
    e22 = sort.by.object.size(expand.as.product(e2))
    if (length(e11)>1 && length(e11)==length(e22)){
        all.eq = T
        for(i in 1:length(e11)){
            if (!Recall(e11[[i]] , e22[[i]])) {
                all.eq=F
                break
            }
        }
        if (all.eq) return(T)
    }
    F
}

#' expand.interaction
#'
#' more deeper expansion than \code{\link{expand.as.sum.of.product}} \cr
#' can do \code{ (a+b)*c*(d+e) -> acd +ace + bcd + bce }
#'
#' @param e expression
#' @return expanded expression
#' @author jjxie
expand.interaction = function(e){
    l0 = expand.as.product(e)
    le = lapply(l0, expand.as.sum)
    ind = sapply(le, length)
    N = length(ind)
    e0 = 0
    walker = function(top, e1) {
        if (top>N) {
            e0 <<- CONS('+',e0, e1)
            return(NULL)
        }
        for(i in 1:ind[top]){
            Recall(top+1, simplify.2(CONS('*',e1, le[[top]][[i]])))
        }
        return(NULL)
    }
    walker(1,1)
    simplify.2(e0)
}

#' collect.algebraical.1
#'
#' \code{ a*a*a + b*b --> a^3 + b^2 } \cr
#' \code{ (a+b)*(a+b) --> a^2 + 2*a*b + b^2 }
#'
#' @param e expression
#' @return expression
#' @author jjxie
collect.algebraical.1 = function(e){
    m = exponential.matrix.of.monomials.sum.product(expand.as.sum.of.product(expand.interaction(e)))
    e0 = 0
    dg = sapply(m$mono, function(x) 1-as.numeric(is.atomic(x)))
    sz = sapply(m$mono, object.size)
    nm = sapply(m$mono, function(x) deparse(x)[1])
    ord = order(dg,sz,nm)
    for(i in 1:NROW(m$mat)){
        e1 = 1
        for(j in ord) {
            if (m$mat[i,j]==0) next
            if (m$mat[i,j]==1) {
                e1 = simplify(CONS('*',e1, m$mono[[j]]))
                next
            }
            e1 = simplify(CONS('*',e1,CONS('^',m$mono[[j]],m$mat[i,j])))
        }
        e0 = simplify(CONS('+',e0,e1))
    }
    e0
}

#' exponential.matrix.of.monomials.sum.product 
#' 
#' For a list of list as following: \cr
#' \code{ le : [ [ A^2 ,B, A] ,...] } \cr
#' generate the following matrix: \cr
#' \verb{
#'    .   A  B  ...  
#'    1   3  1  0   
#'    2   ......... 
#' }\cr
#' each row for a entry of list \code{le} \cr
#' each column is a monomial appearing in the whole list\cr
#' each entry of matrix is the degree of the monomial
#'
#' @param le list of list of expression
#' @return matrix of exponentials
#' @author jjxie
exponential.matrix.of.monomials.sum.product = function(le){
    mono = list()
    N.mono = 0
    in.mono.ind = function(e){
        if(N.mono<1)return(0)
        for(i in 1:N.mono){
            if(equal.symbol.light(e, mono[[i]])) return(i)
        }
        0
    }
    push = function(e){
        mono[[ (N.mono<<-N.mono+1) ]] <<- e
        N.mono
    }
    # even no common divisor, at most sum of length of each 
    mat = matrix(0, nrow=length(le),ncol=sum(sapply(le,length)))
    for(i in 1:length(le)){
        for(j in 1:length(le[[i]])){
            tmp = le[[i]][[j]]
            if (symbolic.match(quote('?'(XX)^'?'(YY)), tmp)) {
                e0 = XX
                e1 = YY
            } else {
                e0 = tmp
                e1 = 1
            }
            ind = in.mono.ind(e0)
            if (ind==0) {
                ind=push(e0)
            }
            mat[i, ind] = mat[i,ind] + e1
        }
    }
    mat = mat[,1:N.mono,drop=FALSE]
    list(mat=mat,mono=mono)
}

#' factorize.3
#'
#' better than factorize.2,
#' can recognize A A^k has common divisors
#'
#' @param e expression
#' @return an expression as production
#' @author jjxie
factorize.3 = function(e){
    e = simplify(e)
    if (is.atomic(e)) return(e)
    le = expand.as.sum.of.product(e)
    if (length(le)==1 && length(le[[1]])==1) return(le[[1]][[1]])
    m = exponential.matrix.of.monomials.sum.product(le)
    lcd = apply(m$mat,2,min)
    N.mono=length(m$mono)
    mat1 = m$mat - rep(lcd, each=NROW(m$mat))
    lcd.e = simplify.2(expression.from.list.of.product(sapply(1:N.mono, function(i) CONS('^',m$mono[[i]], lcd[i]))))
    l = vector(NROW(mat1),mode='list')
    for(i in 1:NROW(mat1)){
        l[[i]] = expression.from.list.of.product(sapply(1:N.mono, function(j) CONS('^',m$mono[[j]], mat1[i,j])))
    }
    e1 = simplify.2(expression.from.list.of.sum(l))
    if (is.atomic(lcd.e)) {
        if (lcd.e==1) return(e1)
        return( CONS('*',lcd.e, e1) )
    }
    # some non-trivial divisor got, try to factor again
    simplify.2(CONS('*', Recall(e1), lcd.e))
}

#' factorize.2
#'
#' very light weighted factorization \cr
#' just factor out commond terms \cr
#' however, we haven't used any division algorithm's here, \cr
#' hence even can't recognized \code{A|A^k}
#'
#' @param e expression
#' @return an expression as production
#' @author jjxie
factorize.2 = function(e){
    le = expand.as.sum.of.product(e)
    N = length(le)
    if (N==1) {
        return(expression.from.list.of.product(le[[1]]))
    }
    #####
    terms.lengths = sapply(le, length)
    # we test each term in the shortest summand 
    looking.at.from = which.min(terms.lengths)
    # have to deal with following case:
    #   [ [A,A], [A,B],[A,C] ]
    # should carefully count the used terms
    used.term = matrix(F, nrow=N, ncol=max(terms.lengths))
    # least common divisor's indices
    lcd = c()
    for(i0 in 1:length(le[[looking.at.from]])){
        if (is.atomic(term <- le[[looking.at.from]][[ i0 ]])) next
        # of course, do not need to factor out a real number
        found.in.all = T
        # temparary flags: counting which terms are used
        used.term.0 = matrix(F, nrow=N, ncol=max(terms.lengths))
        for(i in 1:N){
            if (i == looking.at.from) next
            found.in.i = F
            for(j in 1:length(le[[i]])){
                if (used.term[i,j]) next
                if (equal.symbol.light(le[[i]][[j]], term)) {
                    used.term.0[i,j]=T
                    found.in.i = T
                    break
                }
            }
            if (!found.in.i){
                found.in.all = F
                break
            }
        }
        if (found.in.all){
            lcd = c(lcd, i0)
            # only if found.in.all other summands, we merge back used.term.0 to used.term
            used.term = used.term | used.term.0
            used.term[looking.at.from, i0] = T
        }
    }
    #######
    if (length(lcd)==0) return(e)
    lcd.exp = expression.from.list.of.product(le[[looking.at.from]][ lcd ])
    le1 = list()
    for(i in 1:N){
        le1[[i]] = expression.from.list.of.product(le[[i]][ ! (used.term[i, 1:terms.lengths[i]]) ])
    }
    e1 = Recall(simplify.2(expression.from.list.of.sum(le1)))
    as.call(c(as.symbol('*'),e1,lcd.exp))
}

##############################################################################################################
# PRINT just for debug
##############################################################################################################

#' print.listoflist
#'
#' Just for debuging purpose, pretty printer for list of list\cr
#' Will display list of list as \code{[ [A11, A12, A13, ... ] , [A21, A22, A23, ... ], ... ]}
#'
#' @param x list of list of expressions
#' @return string
#' @author jjxie
print.listoflist = function(x){
    tmp = 
    sapply(x, FUN=function(y) {
        sprintf('[ %s ]', paste(sapply(y,deparse),collapse=', ')) }
    )
    re = sprintf('[ %s ]\n',paste(tmp , collapse=', '))
    cat(re)
    invisible(re)
}

##############################################################################################################
# eval functions
##############################################################################################################

#' eval.symbolic0
#' evaluate a sequence of equations
#'
#' \code{ TVKA = THETA[1]*WT/70 }\cr
#' \code{ KA = TVKA * exp(ETA[1]) } \cr
#' 
#' can be evaluated to \code{ KA= THETA[1]*WT/70 * exp(ETA[1]) }
#'
#' @param eqns equations (list of expressions)
#' @param eval.env the environment for evaluation
#' @return the last value of evaluation
#' @author jjxie
eval.symbolic0 = function(eqns, eval.env = new.env(parent=emptyenv())){
    walker = function(e) {
        if (is.atomic(e)) return(e)
        if (is.symbol(e)) {
            a = as.character(e)
            if  (!exists(a,envir=eval.env)) return(e)
            return(get(a, env=eval.env))
        }
        if (is.call(e) && is.function(e[[1]])) {
            stop('Found R closure, not yet supported!')
        }
        if (is.call(e) && length(e)==2 && e[[1]]=='-' && is.numeric(e[[2]])) {
        # quote( - 1 ) is CONS('-',1)
        # we should ensure it's atomic
            return( -e[[2]] )
        }
        if (is.call(e) && length(e)==2 && e[[1]]=='(') {
            return( Recall(e[[2]]) )
        }
        if (e[[1]]=='..FROZEN') {
            # add in evaluation engine a frozen token (for further use)
            return(e)
        }
        # There is no '->' token in R, 
        # After parsing, '->' is already changed into '<-'
        if (e[[1]]=='=' || e[[1]]=='<-') {
            # lvalue should not be evaluated!
            lhs = e[[2]]
            rhs = Recall(e[[3]])
            # lhs can be symbol or Array reference
            if (is.symbol(lhs)) {
                assign(as.character(lhs), rhs ,env=eval.env)
                return(rhs)
            }else{
                # R has [ and [[
                # however, to hold the symbol, must be a list
                # and should not allow multple indexing here
                if (is.call(lhs) && (lhs[[1]] == '[' || lhs[[1]]=='[[')) {
                    array.name = as.character(lhs[[2]])
                    dims = length(lhs) - 2
                    stopifnot(dims >=1 )
                    inds = character(dims)
                    for(i in 1:dims) {
                        inds[i] = as.character(Recall(lhs[[2+i]]))
                    }
                    array.name = paste( c('._____',array.name, inds) , collapse=',' )
                    assign(array.name, rhs, env=eval.env)
                    return(rhs)
                }
                stop(sprintf('lhs being %s is not supported', paste(deparse(lhs),collapse='')))
            }
        }
        if (e[[1]]=='[' || e[[1]]=='[['){
            array.name = as.character(e[[2]])
            dims = length(e) - 2
            stopifnot(dims >=1 )
            inds = character(dims)
            for(i in 1:dims) {
                inds[i] = as.character(Recall(e[[2+i]]))
            }
            array.name = paste( c('._____',array.name, inds) , collapse=',' )
            if(!exists(array.name,envir=eval.env)) return(e)
            return(get(array.name,envir=eval.env))
        }
        if (e[[1]]=='simplify'){
            return(simplify.2(Recall(e[[2]])))
        }
        for(i in 1:length(e)){
            # symbol replacement
            e[[i]]=Recall(e[[i]])
        }
        e
    }
    for(e0 in eqns){
        re = walker(e0)
    }
    re
}

#' eval.symbolic
#'
#' More handy interface for \code{\link{eval.symbolic0}} \cr
#' input can now be a expression or several code in \code{{}}'s
#'
#' @param e expression or list of expression
#' @param eval.env evaluation environment
#' @return last value of evaluation
#' @author jjxie
eval.symbolic = function(e, eval.env=new.env(parent=emptyenv())){
    if (is.language(e)) {
        if (is.call(e)) { 
            if ( e[[1]] == '{') {
                return(eval.symbolic0(as.list(e[-1]), eval.env=eval.env))
            }
            return(eval.symbolic0(list(e), eval.env=eval.env))
        } 
        if (is.expression(e)){
            return(eval.symbolic0(as.list(e), eval.env=eval.env))
        }
    }
    if (is.list(e)) {
        return(eval.symbolic0(e, eval.env=eval.env))
    }
    eval.symbolic0(list(e), eval.env=eval.env)
}

#' @title compute.symbolic
#' This is example for compute the symbolic value given several equations\cr
#' implementation is easy, because we have evaluation engine already, we can even \cr
#' write program in that "language" \cr
#' here we just need to evaluate the symbol in the populated environment
#' @param eqns equations, nonmem characters or just R's expressions
#' @return list, with the name
compute.symbolic = function(eqns, v){
    e0 = new.env(parent=emptyenv())
    if (is.character(eqns)) eqns = nonmem.eqns.to.r.eqns(eqns)
    eval.symbolic(eqns,e0)
    L = lapply(v, function(x)eval.symbolic(as.symbol(x),e0))
    names(L)=v
    L
}

#' Change a symbol like CL or a call like A[1] to a character
#' @param x symbol or call
#' @param keep.original if save the original call/symbol in attribute
#' @return the character
asCharacterSymbolCall = function(x, keep.original=FALSE) {
    if (is.character(x)) return(x)
    re = deparse(x)
    if (keep.original) attr(re,'original') = x
    re
}
