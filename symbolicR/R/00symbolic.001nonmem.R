# nonmem related symbolic operations #
######################################

#' fortran2r
#'
#' translate NONMEM operators to R valid operators\cr
#' this function won't check if the pattern is in the quoted case, if your string has \cr
#' \code{"as .EQ. "} , it will be destoried as \code{"as == "}
#'
#' @param x string
#' @return string
#' @author Mango solutions
fortran2r = function(x){
    x=gsub('.EQ.','==',x,fixed=T)
    x=gsub('.NE.','!=',x,fixed=T)
    x=gsub('.GT.','>',x,fixed=T)
    x=gsub('.GE.','>=',x,fixed=T)
    # avoid the possibility  A<-1, this is parsed into A <- 1, however, we want it to be A< -1
    x=gsub('.LT.','< ',x,fixed=T)
    x=gsub('.LE.','<=',x,fixed=T)
    x=gsub('.AND.','&&',x,fixed=T)
    x=gsub('.OR.','||',x,fixed=T)
    # no need to change this, R will parse ** to ^ automatically
    #x=gsub('**','^',x,fixed=T)
    x
}

#' parse.if
#'
#' parse and convert NONMEM \code{if/then/else} (might be nested) Statements \cr
#' into valid R equivlent one line statement 
#'
#' multiple lines or \code{ELSEIF} will be changed into R's statements using braces \code{{}}
#'
#' @param txts a character vector containing lines of NONMEM abbreviative code
#' @return a list representing parsed structure,  \cr
#'          each entry in the list has a \code{type} entry(recording the type of this entry), \cr
#'          a \code{val} entry (recording a string which is R valid and parseable) \cr
#'          and a \code{pos} (the starting position in the original stream \code{txt}) entry
#'
#' @author Mango solutions
parse.if = function(txts){
    # txts is a character vector
    pos = 0
    NL  = length(txts)
    word = NULL
    # aux functions
    chop.spaces = function(s){
      s = sub('^[[:space:]]+','',s)
      sub('[[:space:]]+$','',s)
    }
    # stream fun
    nextl = function(){
        if (pos + 1 > NL ) {
          return(NULL)
        }
        txts[pos + 1]
    }
    consum= function(){
        if (!is.null(foo<-nextl())) {
          pos <<- pos +1
        }
        foo
    }
    ######
  token.dilim.space = function(){
    x = nextl()
    if (regexpr('^$',x) <= 0) return(FALSE)
    consum()
    TRUE
  }
  token.if.line = function(){
    x = nextl()
    if (regexpr('^\\s*IF\\s*\\(',x) <= 0) return(FALSE)
    y = extractBalanced(x)
    if (regexpr('^\\s*THEN',y[2]) > 0) return(FALSE) 
    word <<- sprintf('if %s %s', fortran2r(y[1]), fortran2r(y[2]))
    consum()
    TRUE
  }
  token.if.then = function(){
    x = nextl()
    if (regexpr('^\\s*IF\\s*\\(',x) <= 0) return(FALSE)
    y = extractBalanced(x)
    if (regexpr('^\\s*THEN',y[2]) < 0) return(FALSE) 
    word <<- sprintf('if %s ', fortran2r(y[1]))
    consum()
    TRUE
  }
  token.else = function(){
    x = nextl()
    if (regexpr('^\\s*ELSE\\s*$',x) <= 0) return(FALSE)
    consum()
    TRUE
  }
  token.elseifthen = function(){
    x = nextl()
    if (regexpr('^\\s*ELSEIF\\s*',x) <= 0) return(FALSE)
    y = extractBalanced(x,leadingChar='ELSEIF\\s*')
    if (regexpr('^\\s*THEN',y[2]) < 0) return(FALSE) 
    word <<-  fortran2r(y[1])
    consum()
    TRUE
  }
  token.endif = function(){
    x = nextl()
    if (regexpr('^\\s*ENDIF\\s*$',x) <= 0) return(FALSE)
    consum()
    TRUE
  }
  # stack fun
  stack <- vector(mode='list',10)
  stack.top <- 0
  push         = function(entry){
    stack.top <<- stack.top + 1
    stack[[ stack.top ]] <<- entry
  }
  pop          = function(){
    if (stack.top < 1) stop('stack overflow')
    re = stack[[ stack.top ]]
    stack.top <<- stack.top - 1
    re
  }
  clear.stack  = function(){
    stack.top <<- 0
  }
#
  #### rules, and rule functions
  rules = list(
    list(lhs='if.cond.then',
         rhs=list( # alternatives:
                # r1
                list(pattern=c('token.if','token.rcode'),
                     action =function(){
                         re = matchedargs[[1]]
                         re$val = sprintf('%s { %s; ',re$val, matchedargs[[2]]$val)
                         re
                     }),
                # r2
                list(pattern=c('if.cond.then','token.rcode'),
                     action =function(){
                         re = matchedargs[[1]]
                         re$val = sprintf('%s %s; ',re$val, matchedargs[[2]]$val)
                         re
                     }),
                # r3
                list(pattern=c('token.if','rcode'),
                     action =function(){
                         re = matchedargs[[1]]
                         re$val = sprintf('%s { %s; ',re$val, matchedargs[[2]]$val)
                         re
                     }),
                # r4
                list(pattern=c('if.cond.then','rcode'),
                     action =function(){
                         re = matchedargs[[1]]
                         re$val = sprintf('%s %s; ',re$val, matchedargs[[2]]$val)
                         re
                     }),
                # r5
                list(pattern=c('if.cond.then','token.elseifthen'),
                     action = function(){
                         re = matchedargs[[1]]
                         re$val = sprintf('%s } else if %s {',re$val, matchedargs[[2]]$val)
                         re
                     })
                # end
         )
    ),
    list(lhs='if.cond.then.else',
         rhs=list( # alternatives:
                # r1
                list(pattern=c('if.cond.then','token.else'),
                     action = function(){
                         re = matchedargs[[1]]
                         re$val = sprintf('%s } else { ', re$val)
                         re
                     }),
                # r2
                list(pattern=c('if.cond.then.else','token.rcode'),
                     action = function(){
                         re = matchedargs[[1]]
                         re$val = sprintf('%s %s;',re$val, matchedargs[[2]]$val)
                         re
                     }),
                # r3
                list(pattern=c('if.cond.then.else','rcode'),
                     action = function(){
                         re = matchedargs[[1]]
                         re$val = sprintf('%s %s;',re$val, matchedargs[[2]]$val)
                         re
                     })
                # r4
        )
    ),
    list(lhs='if.cond.then.endif',
         rhs=list( # alternatives:
                # r1
                list(pattern=c('if.cond.then','token.endif'),
                     action = function(){
                         re = matchedargs[[1]]
                         re$val = sprintf('%s }', re$val)
                         re
                     })
                # ##
        )
    ),
    list(lhs='if.cond.then.else.endif',
         rhs=list( # alternatives:
                # r1
                list(pattern=c('if.cond.then.else','token.endif'),
                     action = function(){
                         re = matchedargs[[1]]
                         re$val = sprintf('%s }', re$val)
                         re
                     })
                # ##
        )
    ),
    list(lhs='rcode',
         rhs=list( # alternatives:
                # r1
                list(pattern=c('if.cond.then.else.endif'),
                     action = function(){
                         re = matchedargs[[1]]
                         re
                     }),
                # r2
                list(pattern=c('if.cond.then.endif'),
                     action = function(){
                         re = matchedargs[[1]]
                         re
                     }),
                # r3
                list(pattern=c('token.rcode'),
                     action = function(){
                         re = matchedargs[[1]]
                         re
                     })
        ))
  )
  ## matcher and simplifier
  matchedargs = list()

  match.rhs = function(pattern){
    np = length(pattern)
    if (np == 0 ) return(TRUE)
    if (stack.top < np) return(FALSE)
    if (!all(sapply(stack[ (stack.top - np + 1) : stack.top ], function(y) y$type) == pattern)) return(FALSE)
    for(i in 1:np) { # for easy access from actions
        matchedargs[[ i ]] <<- stack[[ stack.top - np + i ]]
    }
    TRUE
  }
    
  simplify.to.lhs.by.action = function(lhs,rhs) {
    # pop length(rhs$pattern)
    np = length(rhs$pattern)
    stack.top <<- stack.top - np
    push(rhs$action())
    stack[[ stack.top ]]$type <<- lhs
  }

  match.rule = function() {
    # match rules
    for(rule in rules){
        for(rhs in rule$rhs){
            if (match.rhs(rhs$pattern)) {
                simplify.to.lhs.by.action(rule$lhs, rhs)
                return(TRUE)
            }
        }
    }
    FALSE
  }

  ## reducer
  reduce.stack = function(){
    # try rules until no rules can be used
    # then back to the main loop for next token
    anyreduced = FALSE
    while( stack.top>0 ){
        if (!match.rule()) return(anyreduced)
        anyreduced=TRUE
    }
    anyreduced
  }
  # main loop : turn preToken to Token
  while(!is.null(nextl())){
    # First emit TOKENs
    if (token.dilim.space()) {
      next
    } else if (token.if.line()) {
      push(list(type='token.rcode',val=word,pos=pos))
    } else if (token.if.then()) {
      push(list(type='token.if',val=word,pos=pos))
    } else if (token.elseifthen()) {
      push(list(type='token.elseifthen',val=word,pos=pos))
    } else if (token.else()) {
      push(list(type='token.else',val=NULL,pos=pos))
    } else if (token.endif()) {
      push(list(type='token.endif',val=NULL,pos=pos))
    } else {
      # we just treat it as valid r code
      push(list(type='token.rcode', val=fortran2r(consum()), pos=pos))
    }
    reduce.stack()
  }
  re = stack[1:stack.top]
  re
}

#' remove those constant initialization using ifelse
simplify.nonmem.eqns.constant.initialization = function(e){
    if (is.call(e) && e[[1]]=='=' && symbolic.match(quote(ifelse(NEWIND<=1,'?'(const.exp),NA )),e[[3]])) {
        if( length( extract.vertices(const.exp)$rhs.symbol) == 0) {
            e[[3]] = const.exp
        }
    }
    e
}

#' nonmem.parameters.specialarray.to.r.array.rules
#'
#' rules for changing \code{THETA(i),ETA(i)} to \code{THETA[i],ETA[i]}, \cr
#' both \code{EPS(i),ERR(i)} are changed to \code{EPS[i]} \cr
#' 
#' possible NONMEM(fortran) arrays are changed to R representations, including \cr
#' \code{A, DADT} etc.
#'
#' NONMEM(fortran) 's function are changed to R equivlent, e.g. \code{LOG10(x)} to \code{log10(x)}
#'
#' this is a rule definition
#' @seealso \code{\link{character.rules.to.list.rules}}, \code{\link{symbolic.simplify.gigo}}, \code{\link{nonmem.parameters.specialarray.to.r.array}}
nonmem.parameters.specialarray.to.r.array.rules = 
    character.rules.to.list.rules( c(
        " THETA('?a'(i)) ",                         " THETA[':'(i)] ",
        " ETA('?a'(i)) ",                           " ETA[':'(i)] ",
        " EPS('?a'(i)) ",                           " EPS[':'(i)] ",
        " ERR('?a'(i)) ",                           " EPS[':'(i)] ",
        " A_0('?a'(i)) ",                           " A_0[':'(i)] ",
        " DADT('?a'(i)) ",                          " DADT[':'(i)] ",
        " A('?a'(i)) ",                             " A[':'(i)] ",
        " '?s'(X)('?a'(i))='?'(Y) ",                " ':'(X)[':'(i)] = ':'(Y) ",
        " '?s'(X)('?a'(i),'?a'(j))='?'(Y) ",        " ':'(X)[':'(i),':'(j)] = ':'(Y) ", # two dimensional
        " '?s'(X)('?a'(i),'?a'(j),'?a'(k))='?'(Y) "," ':'(X)[':'(i),':'(j),':'(k)] = ':'(Y) " ,# three dimensional
        " LOG10('?'(X)) " ,                         " log10(':'(X)) ",
        " LOG('?'(X)) " ,                           " log(':'(X)) ",
        " DLOG('?'(X)) " ,                           " log(':'(X)) ",
        " DLOG10('?'(X)) " ,                           " log(':'(X)) ",
        " EXP('?'(X)) " ,                           " exp(':'(X)) ",
        " DEXP('?'(X)) " ,                           " exp(':'(X)) ", # old FORTRAN77 style, should use EXP (generic) now in new FORTRAN's
        " SIN('?'(X)) " ,                           " sin(':'(X)) ",
        " DSIN('?'(X)) " ,                           " sin(':'(X)) ",
        " COS('?'(X)) " ,                           " cos(':'(X)) ",
        " DCOS('?'(X)) " ,                           " cos(':'(X)) ",
        " SQRT('?'(X)) " ,                          " sqrt(':'(X)) ",
        " DSQRT('?'(X)) " ,                          " sqrt(':'(X)) ",
        " ABS('?'(X)) " ,                           " abs(':'(X)) ",
        " TAN('?'(X)) " ,                           " tan(':'(X)) ",
        " DTAN('?'(X)) " ,                           " tan(':'(X)) ",
        " ASIN('?'(X)) " ,                          " asin(':'(X)) ",
        " DASIN('?'(X)) " ,                          " asin(':'(X)) ",
        " ACOS('?'(X)) " ,                          " acos(':'(X)) ",
        " DACOS('?'(X)) " ,                          " acos(':'(X)) ",
        " ATAN('?'(X)) " ,                          " atan(':'(X)) ",
        " DATAN('?'(X)) " ,                          " atan(':'(X)) ",
        " INT('?'(X)) " ,                           " floor(':'(X)) "
))

#' nonmem.parameters.specialarray.to.r.array
#'
#' simplifier generated by \code{\link{nonmem.parameters.specialarray.to.r.array.rules}}
#' @param e expression
#' @return expression
nonmem.parameters.specialarray.to.r.array = symbolic.simplify.gigo(nonmem.parameters.specialarray.to.r.array.rules)

#' nonmem.r.expression.to.nonmem.rule
#'
#' r expression to nonmem, still expression to expression, Array are kept as r Array
#' It's the inverse direct of \code{\link{nonmem.parameters.specialarray.to.r.array}}
nonmem.r.expression.to.nonmem.rule = character.rules.to.list.rules( c(
        " log10('?'(X)) " ,                         " LOG10(':'(X)) ",
        " log('?'(X)) " ,                           " LOG(':'(X)) ",
        " exp('?'(X)) " ,                           " EXP(':'(X)) ",
        " sin('?'(X)) " ,                           " SIN(':'(X)) ",
        " cos('?'(X)) " ,                           " COS(':'(X)) ",
        " sqrt('?'(X)) " ,                          " SQRT(':'(X)) ",
        " abs('?'(X)) " ,                           " ABS(':'(X)) ",
        " tan('?'(X)) " ,                           " TAN(':'(X)) ",
        " asin('?'(X)) " ,                          " ASIN(':'(X)) ",
        " acos('?'(X)) " ,                          " ACOS(':'(X)) ",
        " atan('?'(X)) " ,                          " ATAN(':'(X)) ",
        " floor('?'(X)) " ,                           " INT(':'(X)) "
))

#' batch converting from the NONMEM PK like code to r valid code
#' 
#' some special change is applied that the following analysing code can run correctly 
#'
#' @param eqns strings for nonmem equations
#' @return equation as r expressions
nonmem.eqns.to.r.eqns = function(eqns){
    eqns = sapply(parse.if(eqns), function(x) x$val)
    eqns = lapply(eqns, function(x) {
        simplify(nonmem.parameters.specialarray.to.r.array(parse1(x))) })
    eqns = eliminate.loop.in.eqns( convert.ifelse.statement.to.function( eqns ) )
    eqns = eliminate.constant.initilization.in.eqns(eqns)
    eqns = lapply(eqns, simplify.nonmem.eqns.constant.initialization)
    # eliminate.constant.initilization again : LN2 = log(2), replace log(2) everywhere
    eliminate.constant.initilization.in.eqns(eqns)
}

#' rexpression.to.nonmem.expression
#'
#' generated by \code{\link{nonmem.r.expression.to.nonmem.rule}}
#'
#' @param e expression
#' @return e expression
rexpression.to.nonmem.expression = symbolic.simplify.gigo( nonmem.r.expression.to.nonmem.rule )

#' rexpression.to.string.nonmem
#'
#' R expression to strings which is valid in NONMEM control stream \cr
#' \itemize{
#'  \item logical e.g. \code{ .EQ.  => == }
#'  \item braces \code{ {} } removed
#'  \item if-else re-organized as NONMEM control stream
#'  \item ARRAY \code{[]} is changed into \code{()}
#'  \item exponential \code{^} is changed into \code{**}
#'  \item redundant parentheses are removed, e.g. \code{ a + (b*c) } \samp{->} \code{ a + b*c }
#' }
#' @param e0 R expression
#' @return a single string including <CR>
rexpression.to.string.nonmem = function(e0){
    not.operator.or.array = function(s){
        if (length(s)==1 || is.symbol(s) || is.atomic(s)) return(T)
        # must be a call
        s = s[[1]]
        if (s=='(') return(T)
        if (!is.symbol(s)) return(F)
        s = as.character(s)
        if (nchar(s)<1) return(F)
        s = substr(s,1,1)
        s == '['              ||
        (s >= 'A' && s <='Z') ||
        (s >= 'a' && s <='z')
    }
    wrapper = function(s){
        if (not.operator.or.array(s)) walker(s)
        else sprintf('(%s)',walker(s))
    }
    indent.n = function(ss,n) {
       if (n==0) return(ss)
       indent = paste(rep(' ',n),collapse='')
       ss = gsub('^', indent, ss)
       gsub('\n', sprintf('\n%s',indent), ss)
    }
    walker = function(e,indent.w = 0){
        if (is.atomic(e)) {
            if (is.numeric(e)) {
                if (is.na(e)) e = 'NAN'
                else e = toupper(e)
            }else if(is.logical(e)){
                e = sprintf('.%s.',toupper(e))
            }
            return(e)
        }
        if (is.symbol(e)) return(as.character(e))
        if (e[[1]]=='(') {
            return(sprintf('(%s)',Recall(e[[2]])))
        }
        if (e[[1]]=='=') {
            return(sprintf('%s=%s', Recall(e[[2]]), Recall(e[[3]])))
        }
        if (e[[1]]=='==') {
            return(sprintf('(%s.EQ.%s)', Recall(e[[2]]), Recall(e[[3]])))
        }
        if (e[[1]]=='!=') {
            return(sprintf('(%s.NE.%s)', Recall(e[[2]]), Recall(e[[3]])))
        }
        if (e[[1]]=='>=') {
            return(sprintf('(%s.GE.%s)', Recall(e[[2]]), Recall(e[[3]])))
        }
        if (e[[1]]=='>') {
            return(sprintf('(%s.GT.%s)', Recall(e[[2]]), Recall(e[[3]])))
        }
        if (e[[1]]=='<=') {
            return(sprintf('(%s.LE.%s)', Recall(e[[2]]), Recall(e[[3]])))
        }
        if (e[[1]]=='<' ){
            return(sprintf('(%s.LT.%s)', Recall(e[[2]]), Recall(e[[3]])))
        }
        if (e[[1]]=='&&' ){
            return(sprintf('(%s.AND.%s)', Recall(e[[2]]), Recall(e[[3]])))
        }
        if (e[[1]]=='||' ){
            return(sprintf('(%s.OR.%s)', Recall(e[[2]]), Recall(e[[3]])))
        }
        if (e[[1]]=='+'&& length(e)==3) {
            e3 = if (not.operator.or.array(e[[3]]) || e[[3]][[1]]=='*' || e[[3]][[1]]=='^') Recall(e[[3]])
                 else sprintf('(%s)', Recall(e[[3]]))
            return(sprintf('%s+%s', Recall(e[[2]]), e3))
        }
        if (e[[1]]=='-'&& length(e)==3) {
            e3 = if (not.operator.or.array(e[[3]]) || e[[3]][[1]]=='*'|| e[[3]][[1]]=='^') Recall(e[[3]])
                 else sprintf('(%s)', Recall(e[[3]]))
            return(sprintf('%s-%s', Recall(e[[2]]), e3))
        }
        if (e[[1]]=='*'&& length(e)==3) {
            e2 = if (not.operator.or.array(e[[2]]) || e[[2]][[1]]=='*' || e[[2]][[1]]=='^') Recall(e[[2]])
                 else sprintf('(%s)', Recall(e[[2]]))
            e3 = if (not.operator.or.array(e[[3]]) || e[[3]][[1]]=='*' || e[[3]][[1]]=='^') Recall(e[[3]])
                 else sprintf('(%s)', Recall(e[[3]]))
            return(sprintf('%s*%s', e2,e3))
        }
        if (e[[1]]=='^' && length(e)==3) {
            return(sprintf('%s**%s', wrapper(e[[2]]),wrapper(e[[3]])))
        }
        if (e[[1]]=='/' && length(e)==3) {
            return(sprintf('%s/%s', wrapper(e[[2]]),wrapper(e[[3]])))
        }
        if (e[[1]]=='-' && length(e)==2){
            e2 = if (not.operator.or.array(e[[2]]) || e[[2]][[1]]=='*') Recall(e[[2]])
                 else sprintf('(%s)', Recall(e[[2]]))
            return(sprintf('-%s', e2))
        }
        # Array to NONMEM array
        if (e[[1]]=='[') {
            e2 = Recall(e[[2]])
            tmp = NULL
            for(i in 3:length(e)) {
                tmp = c(tmp, Recall(e[[i]]))
            }
            tmp = paste(tmp, collapse=',')
            return( sprintf('%s(%s)',  e2, tmp))
        }
        # '{'
        if (e[[1]]=='{') {
            if (length(e)==1) return('')
            tmp = NULL
            for(i in 2:length(e)){
                tmp = c(tmp, Recall(e[[i]]))
            }
            tmp = paste(tmp, collapse='\n')
            return(tmp)
        }
        # if
        if (e[[1]]=='if' && length(e)==3) {
            # IF THEN
            # ENDIF
            return(paste(  
                sprintf('IF%sTHEN',Recall(e[[2]])), 
                    indent.n(Recall(e[[3]],indent.w+4),indent.w+4),
                'ENDIF', sep='\n'))
        }
        if (e[[1]]=='if' && length(e)==4) {
            e4 = simplify(e[[4]])
            if (is.call(e4) && e4[[1]]=='if'){
                # IF THEN
                # ELSEIF  THEN
                # ELSE
                # ENDIF
                return(paste(  
                    sprintf('IF%sTHEN',Recall(e[[2]])), 
                        indent.n(Recall(e[[3]],indent.w+4),indent.w+4),
                    sprintf('ELSE%s',
                        Recall(e[[4]],indent.w)), 
                    sep='\n'))
            } else {
                # IF THEN
                # ELSE
                # ENDIF
                return(paste(
                    sprintf('IF%sTHEN',Recall(e[[2]])), 
                        indent.n(Recall(e[[3]],indent.w+4),indent.w+4),
                    'ELSE',
                        indent.n(Recall(e[[4]],indent.w+4),indent.w+4), 
                    'ENDIF', sep='\n'))
            }
        }
        # currently, other things just return as is
        e1 = Recall(e[[1]])
        if (length(e)==1) return(e1)
        tmp = NULL
        for(i in 2:length(e)){
            tmp = c(tmp, Recall(e[[i]]))
        }
        tmp = paste(tmp, collapse=',')
        sprintf('%s(%s)', e1, tmp)
    }
    walker(e0)
}

#' extract.edges.quotient.array
#'
#' From one equation \code{LHS = RHS}
#' extract several symbols in the right hand side \cr
#' and symbol in left hand side \cr
#' the \code{array[i]} \code{array[j]} will be extract just as \code{array[]} \cr
#' i.e. the set of edges are quotient by the equivlent relation of glueing together all \cr
#' array entries into a single array name, hence the vertex.table are only characters 
#'
#' @param e0 parsed R expression
#' @return list of lhs and rhs
#' @author jjxie
extract.edges.quotient.array = function(e0){
    lhs.symbol = character(0)
    rhs.symbol = character(0)
    walker = function(e){
        if (is.atomic(e)) return(e)
        if (is.symbol(e)) {
            x = as.character(e)
            rhs.symbol <<- union(rhs.symbol, x)
            return(e)
        }
        if (e[[1]]=='[') {
            # deal with array name specially
            stopifnot(is.symbol(e[[2]]))
            rhs.symbol <<- union(rhs.symbol, sprintf('%s[]',as.character(e[[2]])))
            for(ii in 3:length(e)){
                e[[ii]] = Recall(e[[ii]])
            }
            return(e)
        }
        if (e[[1]]=='=') {
            if (is.atomic(e[[2]])){
                stop('invalid lhs of assignment')
            }
            if (is.symbol(e[[2]])){
                lhs.symbol <<- union(lhs.symbol,as.character(e[[2]]))
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
                lhs.symbol <<- union(lhs.symbol,sprintf('%s[]',as.character(e2[[2]])))
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
        # should record logical expression as 
        # return as is , ignore unknown ones
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

#' create.graph.from.equations
#'
#' For  system of equations, create dependency equations\cr
#' For  \code{A = f(X1,X2,X3,...XN)} \cr
#'   create edges \code{<X1,A>, <X2,A> , ... , <XN,A>} \cr
#' The whold equations will genearte a Graph.
#'
#' However, there might be more than one \code{<Xi,Yi>} edges in the graph, e.g.  \cr
#' \code{ CL = TVCL * ETA(1) } \cr
#' \code{ CL = TVCL * ETA(2) } \cr
#' will generate two \code{<TVCL,CL>} edges, hence the index is important \cr
#' (of course, we do not need to include it explicitly in the return object) \cr
#' when restrict to subgraph, we need only morphism (see \code{\link{restriction.morphism.of.graph.at.v}})
#' 
#' @param lines text lines from NONMEM control stream, might be PK, PRED, ERROR, DES and e.t.c. blocks
#' @return list where \code{eqns.ori} parsed equations,\cr
#'  \code{eqns.lhs} and \code{eqns.rhs} are left and right hand side of equations, \cr
#'  \code{eqns.rel} simplified equations, and tranformed using \code{\link{nonmem.parameters.specialarray.to.r.array}} \cr
#'  \code{edges.ref.eqn} back reference from indices of edges to indices of equations in \code{eqns.rel}, \cr
#'  \code{graph} the graph
#' @author Mango solutions
create.graph.from.equations = function(lines){
    ##
    vertex <- character(0)
    ##
    edges <- list()
    edges.ref.eqn = integer(0)
    N.edges <- line.no <- 0
    eqns.ori <- eqns.lhs <-eqns.rhs <- eqns.rel <- list()
    ##
    for(s in lines){
        if (is.null(s)) next
        line.no <- line.no + 1
        eqns.ori[[ line.no ]] <- eqn <- parse1(s)
        eqn <- simplify(eqn)
        eqn = nonmem.parameters.specialarray.to.r.array(eqn)
        ## the last simplified result
        eqns.rel[[ line.no ]] = eqn
        ##
        l = extract.edges.quotient.array(eqn)
        # get vertex and edges
        eqns.lhs[[ line.no ]] <- lhs <- l$lhs.symbol
        eqns.rhs[[ line.no ]] <- rhs <- l$rhs.symbol
        vertex = union(vertex, lhs)
        vertex = union(vertex, rhs)
        for(i in lhs) {
            for(j in rhs){
                edges[[ (N.edges<-N.edges+1) ]] =  c(j,i)
                edges.ref.eqn[ N.edges ] = line.no
            }
        }
    }
    list(eqns = list( eqns.ori=eqns.ori, eqns.lhs=eqns.lhs, eqns.rhs=eqns.rhs), 
         eqns.rel=eqns.rel,
         edges.ref.eqn=edges.ref.eqn,
         graph = list(vertex=vertex, edges=edges))
}

#'  extractBalanced
#'
#' extract balanced Tokens : () '' "" ,e.t.c just after the \code{leadingChar} \cr
#' used by \code{\link{parse.if}}
#'
#' When there is any nested \code{()}'s, should it be carefully counted.
#'
#' @param stream a string
#' @param leadingChar character pattern
#' @return a character with length 2, first is the extract string, second is the left stream
#' @author Mango solutions
extractBalanced <- function(stream, leadingChar='IF\\s*'){
    orig = stream
    # Assume comments only happend at the end of stream, and remove all possible comments
    # if we assume comments can happen interior of streams, it will be terrificly hard to deal with

    # to process the following case:
    #   ACCEPT = ( A .NE. 'abc;def')
    # before we remove ;.*$, we have to first clear our strings in ""
    # and we also use this to deal with ")" or "("

    i = 1
    lstream = nchar(stream)
    quote.sym = NULL
    while (i <= lstream){
        if (is.null(quote.sym)){
        # not in quote
            if (substring(stream,i,i)=='"' || substring(stream,i,i)=="'" ) {
                quote.sym = substring(stream,i,i) 
                substring(stream,i,i)=' ' # replace string in quotes by ' '
            }
        } else {
            if (substring(stream,i,i)==quote.sym) { 
                # if in quote state and met a pairing one, exit the quote state
                quote.sym = NULL
            }
            substring(stream,i,i)=' ' # replace string in quotes by ' '
        }
        i = i+1
    }

    # comment position
    pos.comment = regexpr(';.*$','', stream)

    # extract the first appearence of leadingChar=(....) in the stream, 
    # we already protect it from "" or ''
    startpos.tag = regexpr(leadingChar, stream)
    if (startpos.tag <= 0 || (pos.comment>0 && startpos.tag >= pos.comment)) {
        # the leading tag is not found in string, keep the original stream untouched
        return( c('', orig) )
    }

    # Just before '('
    i = startpos.tag + attr(startpos.tag, 'match.length') - 1 
    if (pos.comment > 0) {
        lstream = pos.comment - 1
    } else {
        lstream = nchar(stream)
    }
    nextchar = function(){
        if (i > lstream) return(NULL)
        substring(stream,i+1,i+1)
    }
    consumeone = function(){
        re = nextchar()
        if (!is.null(re)){
            i <<- i+1
            return(re)
        }
        re
    }
    contents = ''
    par.bal = 0
    balanced = TRUE
    # according to NONMEM manual, the contents must be surrounded in () pairs
    while(par.bal==0 && !is.null(nextchar())){
        if (nextchar()=='(') 
            par.bal = par.bal + 1
        consumeone()
    }
    if (!is.null(nextchar())){
        # if we do find the starting '('
        startpos.parenthess = i
        while(par.bal>0 && !is.null(nextchar())){
            if (nextchar()=='(') {
                par.bal = par.bal + 1
            } else {
                if (nextchar()==')') {
                    par.bal = par.bal - 1
                }
            }
            consumeone()
        }
        if (substring(stream,i,i) !=')') {
            balanced = FALSE
            warning(sprintf('Unbalanced parenthess found in %s .', orig))
        }
        contents = substring(orig,startpos.parenthess,i)
        orig = paste( substring(orig,1,startpos.tag-1),
                      substring(orig,i+balanced) , sep='')
    }
    c(contents,orig)
}

#' resolve.to.mapping.chain
#'
#' By Graph of Symbol's, we can section up the equations into several stages\cr
#' convert equations to mappings' chain
#' 
#' This is done by using topological sort on the graph, it can discover very basic,\cr
#' rough sections in your PK/PRED code.
#'
#' @param eqns NONMEM equations
#' @return list representing the discovered sections
resolve.to.mapping.chain = function(eqns){
    gg = create.graph.from.equations(sapply(parse.if(eqns),function(x)x$val))
    G = gg$graph
    ind0 = topo.sort(G)
    N.level = max(ind0$N)
    if (N.level<2) stop('Trivial equations.')
    N.mapper = 0
    re = list()
    for(i in 2:N.level){
        v = names(ind0$N)[ind0$N==i]
        v0 = neighbourhood.of.v(G, v)
        subg.ind = restriction.morphism.of.graph.at.v(G, union(v,v0$v.income))
        eqn.ind = unique(gg$edges.ref.eqn[subg.ind])
        # Now remove all those depending on edges outside this subgraph
        edges.outside.subgraph = setdiff(seq_along(G$edges), subg.ind)
        eqn.ind = setdiff(eqn.ind, gg$edges.ref.eqn[ edges.outside.subgraph ] )
        # #
        N.mapper <- N.mapper + 1
        lhs = gg$eqns.lhs[ eqn.ind ]
        re[[ N.mapper ]] =
        list(
            from = v0$v.income,
            to   = v[order(match(v, lhs))],
            type = 'multivar.function',
            relations= as.call(c(as.symbol('{'),gg$eqns.rel[ eqn.ind ]))
        )
    }
    re
}
