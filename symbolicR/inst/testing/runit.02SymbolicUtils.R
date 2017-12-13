# git revision: $Rev: $
# Date of last change: $LastChangedDate: 28/09/2012 $
# Last changed by: $LastChangedBy: ccampbell $
# 
# Original author: jjxie
# Copyright Mango Solutions, Chippenham, UK
###############################################################################

# test aSymbolic.Utils
test.elementary.simplify = function(){

    
    checkEquals(simplify(quote( compute(1) && FALSE )),          FALSE)
    checkEquals(simplify(quote( compute(1) || TRUE  )),          TRUE)
    checkEquals(simplify(quote( exp(THETA[2]) + 0   )),          quote(exp(THETA[2])))
    checkEquals(simplify(quote( sin(EPS[1])*1       )),          quote(sin(EPS[1])))
    checkEquals(simplify(quote( sin(EPS[1])*0       )),          0)
    checkEquals(simplify(quote( 0/sin(EPS[1])^2     )),          0)
    checkEquals(simplify(quote( tan(EPS[1])^1       )),          quote(tan(EPS[1])))
    checkEquals(simplify(quote( (a+b)^2^3           )),          quote( (a+b)^8 ))
    checkEquals(simplify(quote( {(a+b)^2^3}         )),          quote( (a+b)^8 ))
    checkEquals(simplify(quote( a^(-1)              )),          quote( a^-1) )
    checkEquals(simplify(quote( a^(-1)              )),          quote( a^-1) )
    checkEquals(simplify.2(quote( 1 - A             )),          quote( 1 + -1 * A) )
}

test.lhs.rhs = function(){
    checkEquals(extract.edges.quotient.array( parse1(' A = sin(x) + THETA[2] * exp(WT) ' )) , 
        list(lhs.symbol='A',rhs.symbol=c('x','THETA[]','WT')))
    checkEquals(extract.edges.quotient.array( parse1(' A[2] = sin(x) + THETA[2] * exp(WT) ' )) , 
        list(lhs.symbol='A[]',rhs.symbol=c('x','THETA[]','WT')))
    checkException(extract.edges.quotient.array( parse1(' coef(b) = sin(x) + THETA[2] * exp(WT) ' )) , 
        list(lhs.symbol=character(0),rhs.symbol=c('x','THETA[]','WT')),silent=T)
    checkEquals(extract.edges.quotient.array( parse1(' if (TIME>0) X*Y ' )) , 
        list(lhs.symbol=character(0),rhs.symbol=c('TIME','X','Y')))
    checkEquals(extract.edges.quotient.array( parse1(' if (TIME>0 && SEX==2) C=exp(THETA[2])+TVCL' )) , 
        list(lhs.symbol='C',rhs.symbol=c('TIME','SEX','THETA[]','TVCL')))
}

test.pattern.variable = function(){
    
    
    checkTrue(symbolic.match( quote('?'(is.integer,A) * exp(-Y*t)),  quote( 100L * exp(-Y*t)) ) && A==100L)
    checkTrue(symbolic.match( quote( THETA['?a'(II)] * exp( - TIME / V ) ) , 
        quote( THETA[10] * exp( - TIME/ V)) ) && II == 10)
    checkTrue(symbolic.match( quote( '?a'(II) * exp( - TIME / V ) ) , quote( 2.7 * exp( - TIME/ V)) ) && II == 2.7)
    checkTrue(!symbolic.match( quote( '?[]'(ARR)(100) ) , quote( A(100) )))
    checkTrue(symbolic.match( quote( '?()'(CLL)(100) ) , quote( A(b)(100) )) && CLL==quote(A))
    checkTrue(symbolic.match( quote( '?()'(CLL) ) , quote( A(b) )) && CLL==quote(A))
    checkTrue(!symbolic.match( quote( '?()'(CLL)(100) ) , quote( A(100) )))
    checkTrue(symbolic.match( quote( '?s'(CLL)(100)   ) , quote( A(100) ) ) && CLL == quote(A))
}

test.more.than.one.pat.var = function(){
    
    
    checkTrue(symbolic.match( quote( '?'(X) + '?'(Y) - '?'(X) ),  
        quote( sin(exp(THETA[1])) + ETA[12] - sin(exp(THETA[1])))) && Y==quote(ETA[12]))
    checkTrue(!symbolic.match( quote( '?'(X) + '?'(Y) - '?'(X) ),  
        quote( sin(exp(THETA[1])) + ETA[12] - sin(exp(THETA[2])))))
}

test.simplify.sqrt = function(){
    
    
    checkEquals(simplify.sqrt(quote( sqrt( 4 ))),          2)
    checkEquals(simplify.sqrt(quote( sqrt( (XY*Y1)^2  ))),  quote( XY*Y1) )
    checkEquals(simplify.sqrt(quote( log(sqrt( WT^3 )))),  quote( 1.5 * log(WT)) )
}

test.eval.symbolic = function(){
    
    AA <- parse1('A=10')
    BB <- parse1('B= A + 1 ')
    
    
    checkEquals(eval.symbolic(list(AA, BB)), quote(10 + 1))
    checkEquals(eval.symbolic(quote({ TVCL=THETA[1] ; PWR=THETA[5]; CL=TVCL*(WT/70)^PWR })), 
        quote( THETA[1] * (WT/70)^THETA[5]))
}

test.other.0 = function(){
    
    checkEquals(collect.algebraical.1(quote((a+b)*(a+b))), quote(a^2 + 2*a*b+b^2))
    
    checkEquals(expand.as.sum.of.product(quote(a + exp(b) + c*d + e*sin(f))), 
        list( list(quote(a)), 
            list(quote(exp(b))), 
            list(quote(c), quote(d)), 
            list(quote(e), quote(sin(f))))
    )
    
    checkEquals(factorize.3(quote(  a^3*b*THETA[2] + b *a  )), quote( (1 + a^2*THETA[2]) * a * b))
    checkEquals(factorize.3(quote( THETA[1] + THETA[1]*EPS[5] )), quote( (1+EPS[5])*THETA[1]))
}

# test symbolicRNM
test.rexp.to.string = function(){
    
    
    checkEquals(rexpression.to.string.nonmem( quote( a == b ) ), '(a.EQ.b)' )
    checkEquals(rexpression.to.string.nonmem( quote( THETA[1] ) ), 'THETA(1)' )
    checkEquals(rexpression.to.string.nonmem( quote( a+(b*c) ) ), 'a+(b*c)' )
    checkEquals(rexpression.to.string.nonmem( quote( a+b*c ) ), 'a+b*c' )
    checkEquals(rexpression.to.string.nonmem( quote( (a+b)*c)  ), '(a+b)*c' )
    checkEquals(rexpression.to.string.nonmem( quote( if(a){a1}else if(b){b1} else if(c){c1}  )), 
#        "IFaTHEN\na1\nELSEIFbTHEN\nb1\nELSEIFcTHEN\nc1\nENDIF")
        "IFaTHEN\n    a1\nELSEIFbTHEN\n    b1\nELSEIFcTHEN\n    c1\nENDIF")
    checkEquals(rexpression.to.string.nonmem( quote( a+ -10 ) ), 'a+(-10)' )
    checkEquals(rexpression.to.string.nonmem( quote( a^ -10 ) ), 'a**(-10)' )
}

test.parse.if <- function()
{
    
    # parse.if
    
    # TEST 1
    
    test1 <- try(parse.if(c('IF(A.EQ.B)THEN','C','ENDIF')), silent = TRUE)
    checkEquals(test1[[1]]$val, 'if (A==B)  { C;  }', msg = "parse if == statement")
    
    # TEST 2
    
    txt <-  'IF (NEWWIND.NE.2) AMT1=AMT'
    test2 <- try(parse.if(txt), silent = TRUE)
    checkEquals(test2[[1]]$val, "if (NEWWIND!=2)  AMT1=AMT", msg = "parse if != statement")
    
    # TEST 3
    
    txt <- c('IF (TIME.LT.TLAG) THEN',
            'C=0',
        'ELSE',
            'C=(AMT1/V)*(KA/(KA-(CL/V)))*(EXP(-(CL/V)*(TIME-TLAG))-(EXP(-(KA)*(TIME-TLAG))))',
        'ENDIF'
    )
    test3 <- try(parse.if(txt), silent = TRUE)
    checkEquals( test3[[1]]$val, 
        "if (TIME< TLAG)  { C=0;  } else {  C=(AMT1/V)*(KA/(KA-(CL/V)))*(EXP(-(CL/V)*(TIME-TLAG))-(EXP(-(KA)*(TIME-TLAG)))); }")
    
    # TEST 4 : check val
    
    txt <- 'TVCL=THETA(1)*WT/70'
    checkEquals(parse.if(txt)[[1]]$val, txt, msg = "check val")
}


test.create.graph.from.equations <- function()
{
    
    # graph
    
    # TEST 1
    
    test1 <- try(create.graph.from.equations(c( 'a=x + y', 'b=x1 * sin(y2)', 'c=a+THETA(2)')), silent = TRUE)
    checkEquals(test1$graph$edges[[6]], 
        c('THETA[]','c'))
    
    # TEST 2
    
    pks0 = c(
        'TVCL=THETA(1)*WT/70',
        'TVV =THETA(2)',
        'TVKA=THETA(3)',
        'TVLG=THETA(4)',
        'CL=TVCL*EXP(ETA(1))',
        'V =TVV*EXP(ETA(2))',
        'KA=TVKA*EXP(ETA(3))',
        'ALAG1=TVLG*EXP(ETA(4))',
        'S2=V'
    )

    test2 <- create.graph.from.equations(pks0)
    checkEquals(names(test2), c("eqns", "eqns.rel", "edges.ref.eqn", "graph"), msg = "names of graph")
    G  <- test2$graph
    checkEquals(G$vertex, c("TVCL", "THETA[]", "WT", "TVV", "TVKA", "TVLG", "CL", 
        "ETA[]", "V", "KA", "ALAG1", "S2"), msg = "graph vertices")

    
}

eval.symbolic0 <- function()
{
    
    # TEST 1 : combine equations
    
    pks <- c(
    "TVKA = THETA[1]*WT/70",
    "KA = TVKA * exp(ETA[1])")
    test1 <- eval.symbolic0(pks)
    checkEquals(test1, "KA = TVKA * exp(ETA[1])", msg = "combine equations")
    
    # TEST 2 : empty list
    
    test2 <- eval.symbolic0(list())
    checkEquals(test2, character(0), msg = "empty list")
}
