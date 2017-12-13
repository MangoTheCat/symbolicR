#' convert if then else statements to function
#' 
#' But we create some loop here as
#' Q = if (cond) A else Q 
#' i.e. Q -> Q
#'
#' hence we need to remove those loops
#'
ifelse.statement.to.function.rules = list(
# if(cond) lhs1=rhs1
# lhs1 = if(cond) rhs1 else lhs1
list(CONS('if',quote('?'(condition)),
        CONS('=',quote('?'(lhs1)), quote('?'(rhs1)))),
    CONS('=',quote(':'(lhs1)), 
             CONS('ifelse', quote(':'(condition)), quote(':'(rhs1)), quote(':'(lhs1))))),
# if(cond) lhs1=rhs1 else lsh1=rhs2
# lhs1 = if(cond) rhs1 else rhs2
list(CONS('if',quote('?'(condition)),
        CONS('=',quote('?'(lhs1)), quote('?'(rhs1))),
        CONS('=',quote('?'(lhs1)), quote('?'(rhs2)))),
    CONS('=',quote(':'(lhs1)), 
             CONS('ifelse', quote(':'(condition)), quote(':'(rhs1)), quote(':'(rhs2))))),
# if(cond) lhs1=rhs1 else lhs2=rhs2
# lhs1 = if(cond) rhs1 else lhs1
# lhs2 = if(cond) lhs2 else rhs2
list(CONS('if',quote('?'(condition)),
        CONS('=',quote('?'(lhs1)), quote('?'(rhs1))),
        CONS('=',quote('?'(lhs2)), quote('?'(rhs2)))),
    CONS('{',
        CONS('=',quote(':'(lhs1)), 
                 CONS('ifelse', quote(':'(condition)), quote(':'(rhs1)), quote(':'(lhs1)))),
        CONS('=',quote(':'(lhs2)), 
                 CONS('ifelse', quote(':'(condition)), quote(':'(lhs2)), quote(':'(rhs2))))))
)

convert.ifelse.statement.to.function.0 = symbolic.simplify.gigo(ifelse.statement.to.function.rules)
convert.ifelse.statement.to.function = function(eqns){
    eqns = convert.ifelse.statement.to.function.0(eqns)
    re = list()
    for(e in eqns){
        re = c(re, expand.block.to.list(e))
    }
    re
}

#' a expression itselt is changed into \code{list(e)}
#' @param e expression
#' @return a list of equation , any block '{ ... }' is expaned into list
expand.block.to.list = function(e) {
    if (is.call(e) && e[[1]]=='{') {
        l = length(e) - 1
        if (l>0){
            re = list()
            for(i in 1:l){
                re = c(re, Recall(e[[i+1]]))
            }
            return(re)
        }
    }
    list(e)
}

#' add temporary variables to remove loops in equations
#' @param eqns equations
#' @return equations which temporary variables are added to break the loops in generated graphs
eliminate.loop.in.eqns = function(eqns){
    has.loop = function(ss){
        gg = create.abstractGraph.from.equations(ss)
        tmp = try(topo.sort(gg$graph),silent=TRUE)
        is(tmp,'try-error') && grepl('Loop found!', tmp)
    }
    temp.variable.count = 0
    generate.next.temp.variable = function(postfix=''){
        temp.variable.count<<-temp.variable.count+1
        as.symbol(sprintf('TEMP[%s][%s]', temp.variable.count, postfix))
    }
    S0 = list()
    ind = 1
    while(ind <= length(eqns)){
        if (!has.loop(eqns[1:ind])) {
            ind <- ind + 1
            next()
        }
        # only assume eqns[[ind]] is lhs = rhs format
        lhs = eqns[[ ind ]][[2]]
        # loop must be caused by lhs
        lhs.char = paste(deparse(lhs),collapse='')
        # First, check if this lhs is used or initialized by previous equations
        lhs.inited = FALSE
        lhs.used = FALSE
        if (ind > 1) {
            for(i in 1:(ind-1)) {
                if( !lhs.inited) {
                    # lhs is initialized only if 
                    # it's appearing in the left hand side of previous equations
                    lhs.inited = symbolic.match(CONS('=',lhs,quote('?'(XXXX))), eqns[[i]])
                }
                if ( !lhs.used) {
                    tmp = extract.vertices( eqns[[i]])
                    lhs.used = !is.na(match.symbolic(lhs, tmp$rhs.symbol))
                }
                if (lhs.used && lhs.inited) break()
            }
        }
        if (!lhs.inited){
            warning(sprintf('variable [%s] causes recursive reference, 
                and it is not initialized before first use.
                To break this looping use (i.e. inducing loop in dependency graph) of the variable,
                we have to assume it is NA and added an equation to initialize it to NA', lhs.char))
            if (!lhs.used) {
            # if previously not used, we only change the use of lhs in current equation: eqn[[ ind]].rhs to ZERO
                lhs.hat = NA
            } else {
            # however, if it's used, then we need to initialize it, and restart the while loop
                eqns = c(CONS('=', lhs, NA) , eqns)
                ind = 2
            # by adding one initilization code, and then restart
                next()
            }
        } else {
                lhs.hat = generate.next.temp.variable( lhs.char )
                for(i in 1:(ind-1)) {
                    eqns[[i]] = symbolic.gsub(lhs, lhs.hat, eqns[[i]] )
                }
        }
        # also need to change the rhs of current equation
        eqns[[ ind ]][[3]] = symbolic.gsub(lhs, lhs.hat, eqns[[ind]][[3]])
        # restart
        ind = 1
    }
    eqns
}

#' remove constant initilization
#' if a initilization for constant is found, 
#' i.e. an equation which has no rhs symbols; 
#' or say, cannot be back referred by any edges
#' we just substitute the contant value of it to all appearence after it 
#' and remove the equation itself from equations
#' @param eqns equations
#' @return equations whose contant initilization are removed
eliminate.constant.initilization.in.eqns = function(eqns){
    while(TRUE){
        if (length(eqns)==1) {
            # should not remove only one constant assignment
            return(eqns)
        }
        gg = create.abstractGraph.from.equations(eqns)
        to.remove = setdiff(seq_along(eqns), gg$backref)
        if (length(to.remove)==0) return(eqns)
        if (length(extract.vertices(eqns[[ to.remove[1] ]])$lhs.symbol)==0) {
            stop(sprintf('Not an assignment for equation %s', paste(deparse(eqns[[to.remove[1]]]), collapse='')))
        }
        eqns = eliminate.one.equation.constant.initilization.in.eqns(eqns, to.remove[1])
    }
    NULL
}

#' used by eliminate.constant.initilization.in.eqns
#' @param eqns equations
#' @param el.ind the index (only one number is allowed) for which equation to be eliminated
#' @return eliminated equations
eliminate.one.equation.constant.initilization.in.eqns = function(eqns, el.ind) {
    lhs = eqns[[ el.ind ]][[2]]
    rhs = eval.symbolic(eqns[1:el.ind])
    N = length(eqns)
    if (el.ind < N) {
        for(j in (el.ind + 1):N) {
            eqns[[j]][[3]] = symbolic.gsub(lhs, rhs, eqns[[j]][[3]])
        }
    }
    eqns[ - el.ind ]
}

#' Compile \code{ifelse} function call to \code{if} statements
#' @param e expression
#' @param simplify.result if call \code{eliminate.variable.alias.in.eqns}
#' @param remove.na.init if call \code{eliminate.na.initilization.in.eqns}
#' @param temp.var.prefix the prefix of temporary varibles used in output code
#' @return equations which are equivlent to e, using non-nested IF/ELSE
compile.ifelse.to.if.statements = function(e0, simplify.result=TRUE, remove.na.init=TRUE, temp.var.prefix='TMPVAR') {
    eqns = list()
    tempvar.ind = 0
    walker = function(e) {
        if (is.atomic(e) || is.symbol(e)) return(e)
        if (!is.call(e)) stop(sprintf('expression of type: %s not supported' , class(e)))
        if (e[[1]]=='=' && is.call(e[[3]]) && e[[3]][[1]]=='ifelse') {
        # A = ifelse(1,2,3)
        # if(1) A=2 else A=3
            return.symbol = e[[2]]
            e3 = e[[3]]
            cond.part = Recall(e3[[2]])
            yes.part = Recall(e3[[3]])
            no.part = Recall(e3[[4]])
            if (is.atomic(no.part) && is.na(no.part)) {
                re = CONS('if', cond.part,
                            CONS('=', return.symbol, yes.part))
            } else {
                re = CONS('if', cond.part,
                                CONS('=', return.symbol, yes.part),
                                CONS('=', return.symbol, no.part))
            }
            return(re)
        }
        if (e[[1]]=='ifelse') {
        # just a function call ifelse(1,2,3)
        # create a temporary variable TEPVAR001, and return it
        # and append to eqns if(1) TMPVAR001=2 else TMPVAR001=3
            cond.part = Recall(e[[2]])
            yes.part = Recall(e[[3]])
            no.part = Recall(e[[4]])
        # should put temporary variable creation after Recall, to make the temporary appearing in sequence
            return.symbol = as.symbol(sprintf('%s%03d', temp.var.prefix, (tempvar.ind<<-tempvar.ind+1)))
            eqns[[ length(eqns) +1 ]] <<- CONS('=', return.symbol, no.part)
            eqns[[ length(eqns) +1 ]] <<- CONS('if', cond.part, CONS('=', return.symbol, yes.part) )
        # we can also use only one if-else statement as following:
        #    eqns[[ length(eqns)+1 ]] <<- CONS('if', cond.part,
        #                    CONS('=', return.symbol, yes.part),
        #                    CONS('=', return.symbol, no.part))
        # however, it is not looking good
            return(return.symbol)
        }
        if (length(e)>1) {
            for(i in 2:length(e)){
                e[[i]] = Recall(e[[i]])
            }
        }
        e
    }
    tmp = walker(e0)
    eqns[[ length(eqns)+1 ]] = tmp
    if (simplify.result) eqns = eliminate.variable.alias.in.eqns(eqns)
    if (remove.na.init) eqns = eliminate.na.initilization.in.eqns(eqns)
    eqns
}

#' eliminate redundant variables
#' if theres is an alias format as following:
#'
#' \code{A = B}
#' \code{C = A}
#' then, eliminate A in the equations
#' we only deal with \code{symbol = symbol} (alias) pattern in equations
#'
#' Of course, reusing some temporary varible will fool this function
#' but it will work well in simplifying the result of \code{compile.ifelse.to.if.statements}
#'
#' @param eqns equations
#' @return eliminated equations
eliminate.variable.alias.in.eqns = function(eqns){
    while(TRUE){
        if (length(eqns)<=1) return(eqns)
        i = 1
        any.can.be.removed = FALSE
        while(i<=length(eqns)) {
            if (symbolic.match(CONS('=',quote('?s'(lhs)),quote('?s'(rhs))), eqns[[i]])) {
                # lhs = rhs
                # eqns[[ i + 1 ]]
                # eqns[[ i + 2 ]]
                # we can eliminate lhs, only if rhs is not used as right hand side in eqns[ (i+1):length(eqns) ]
                if ( i == length(eqns)) {
                    any.can.be.removed = TRUE
                    break()
                }
                gg = create.abstractGraph.from.equations( eqns[ (i+1):length(eqns) ] )
                rhs.ind = match.symbolic(rhs, gg$vertex.table)
                if (! (rhs.ind %in% gg$rhs.ind) ) {
                    any.can.be.removed = TRUE
                    break()
                }
            }
            i = i + 1
        }
        if (!any.can.be.removed) return(eqns)
        eqns = eliminate.one.equation.variable.alias.in.eqns(eqns, i)
    }
    NULL
}

#' recursively replace right hand side
#' if (1) Y=a else Y=b
#' and try to replace b
#' @param pat pattern
#' @param rep replacement
#' @param e0 expression to be replaced
#' @return new expression
symbolic.replace.rhss = function(pat,rep,e0){
    walker = function(e) {
        if (symbolic.match(pat,e)) {
            return(rep)
        }
        if (is.call(e) && e[[1]]=='=' && length(e)==3) {
            e[[3]] = Recall(e[[3]])
            return(e)
        }
        if (is.call(e) && length(e) > 1) {
            for(i in 2:length(e)) {
                e[[i]] = Recall(e[[i]])
            }
            return(e)
        }
        e
    }
    walker(e0)
}

#' used by eliminate.variable.alias.in.eqns
#' both left hand side appearence or right hand side appearence of lhs of the alias will be 
#' substitute by rhs of alias 
#' hence, to ensure it's semanticly correct, we have to ensure the right hand side of \code{eqns[[el.ind]]}
#' is not used as right hand side in the following equations, 
#' or the value of the right hand side is distoried
#' @param eqns equations
#' @param el.ind the index which equation should be removed
#' @return eliminated equations
eliminate.one.equation.variable.alias.in.eqns = function(eqns, el.ind) {
    lhs = eqns[[ el.ind ]][[2]]
    rhs = eqns[[ el.ind ]][[3]]
    N = length(eqns)
    if (el.ind < N) {
        for(j in (el.ind + 1):N) {
        #    if (is.symbol(eqns[[j]])) {
        #        eqns[[j]] = symbolic.gsub(lhs, rhs, eqns[[j]])
        #    } else {
        #        eqns[[j]][[3]] = symbolic.gsub(lhs, rhs, eqns[[j]][[3]])
        #    }
        #   just change everything , can't assume a = b format
            eqns[[j]] = symbolic.gsub(lhs,rhs,eqns[[j]])
        }
    }
    eqns[ - el.ind ]
}

#' eliminate constant initilization to NA
#' \code{ A = NA }
#' will be eliminated
#'
#' This will make code more beautiful, 
#' based on our assumption that nonmem will initialized user variable to NA
#'
#' @param eqns equations
#' @return eliminated equations
eliminate.na.initilization.in.eqns = function(eqns){
    to.remove = sapply(seq_along(eqns), function(ind) symbolic.match(CONS('=', quote('?s'(XXXX)), NA) , eqns[[ind]]) )
    eqns[!to.remove]
}

#' TO DO
#' The above can not deal with '{' blocks
#' e.g.
#'
#' IF (TYPE.EQ.0) THEN
#' ; PK Data
#'     F_FLAG=0
#'     Y=F+F*ERR(1) ; a prediction
#'  ELSE
#' ; Categorical data
#'     F_FLAG=1
#'     A=DEXP(EXPP)
#'     B=1+A
#'     Y=DV*A/B+(1-DV)/B      ; a likelihood
#'  ENDIF
#' 
#' One possible way to solve this is adding the tuple type
#' Vecterize it
