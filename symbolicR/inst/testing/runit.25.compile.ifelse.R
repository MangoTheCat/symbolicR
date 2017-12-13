test.convert.if.else = function(){
    # rule1
    e = quote( if ( cond1 ) A = B )
    re = convert.ifelse.statement.to.function(list(e))
    checkEquals( re[[1]] , CONS('=',quote(A), quote( ifelse(cond1, B, A) )))
    # rule2
    e = quote( if ( cond1 ) A = B else A=C )
    re = convert.ifelse.statement.to.function(list(e))
    checkEquals( re[[1]] , CONS('=',quote(A), quote( ifelse(cond1, B, C) )))
    # rule3
    e = quote( if ( cond1 ) A = B else C=D )
    re = convert.ifelse.statement.to.function(list(e))
    checkEquals( re[[1]] , CONS('=',quote(A), quote( ifelse(cond1, B, A) )))
    checkEquals( re[[2]] , CONS('=',quote(C), quote( ifelse(cond1, C, D) )))
}

test.eliminate.loop  = function(){
    eqns = list(
        parse1(' a = b '),
        parse1(' b = a + 1')
    )
    re = eliminate.loop.in.eqns(eqns)
    checkEquals(re[[1]], CONS('=',as.symbol('TEMP[1][b]'),NA))
    checkEquals(re[[2]], CONS('=',quote(a), as.symbol('TEMP[1][b]')))
    checkEquals(re[[3]], CONS('=',quote(b), quote( a + 1)))
    re1 = eliminate.constant.initilization.in.eqns(re)
    checkEquals(re1[[1]], CONS('=',quote(b), quote( NA + 1)))
}

test.eliminate.const = function(){
    eqns = list(
        parse1(' tmp1 = 100 '),
        parse1(' b = f(tmp1 , EPS[2])')
    )
    re = eliminate.constant.initilization.in.eqns(eqns)
    checkEquals( re[[1]], parse1('b = f(100, EPS[2])'))
}

test.compile.ifelse = function(){
    e = parse1('a = ifelse(cond1,f(x,y),ifelse(cond2,1,2)) ' )
    re = compile.ifelse.to.if.statements(e)
    checkEquals(re[[1]], parse1('TMPVAR001 = 2'))
    checkEquals(re[[2]], parse1('if (cond2) TMPVAR001 = 1'))
    checkEquals(re[[3]], parse1('if (cond1) a = f(x,y) else a = TMPVAR001'))
}

test.eliminate.alias = function(){
    eqns = list(
        parse1(' a = b'),
        parse1(' c = a')       
    )
    re = eliminate.variable.alias.in.eqns(eqns)
    checkEquals(re[[1]], parse1(' c = b' ))
}
