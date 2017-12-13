# git revision: $Rev: $
# Date of last change: $LastChangedDate: 01/10/2012 $
# Last changed by: $LastChangedBy: ccampbell $
# 
# Original author: ccampbell
# Copyright Mango Solutions, Chippenham, UK
###############################################################################

test.parse1 <- function()
{

    # TEST 1

    txt <- "TVCL = THETA(1) * (WT / 70)**0.75"
    test1 <- parse1(txt)
    checkEquals(test1[[3]][[3]], quote((WT/70)^0.75))

}


test.nonmem.parameters.specialarray.to.r.array <- function()
{    
    
    # nonmem -> r
    
    # TEST 1
    
    test1 <- try(nonmem.parameters.specialarray.to.r.array(quote(100* THETA(1) )), silent = TRUE)
    checkEquals(test1, quote(100*THETA[1]))
    
    # TEST 2
    
    test2 <- try(nonmem.parameters.specialarray.to.r.array(quote(100* DEXP(1) )), silent = TRUE)
    checkEquals(test2, quote(100*exp(1)))
    
    # TEST 3 
    
    test3 <- try(nonmem.parameters.specialarray.to.r.array(parse1('Y = F*EXP(EPS(1)) + EPS(2)')), silent = TRUE)
    # note that this function creates an object for terms in the call
    checkTrue(symbolic.match( CONS('=', quote('?'(x)), quote('?'(y))), test3), msg = "use symbolic match to check whether object has expected structure")
    checkEquals( y , quote( F * exp( EPS[1] ) + EPS[2] ), msg = "was object y created in environment")
    
    # TEST 4 : 
    
    txt <- 'TVCL=THETA(1)*WT/70'
    ptxt <- parse1(txt)
    test4 <- try(nonmem.parameters.specialarray.to.r.array(ptxt[[3]]), silent = TRUE)
    checkEquals(test4, quote(THETA[1] * WT / 70))
    
    # TEST 5 : translate power
    
    txt <- "TVCL = THETA(1) * (WT/70)**0.75"
    ptxt <- parse1(txt)
    test5 <- try(nonmem.parameters.specialarray.to.r.array(ptxt), silent = TRUE)
    checkEquals(test5[[3]], quote(THETA[1] * (WT/70)^0.75), msg = "translate power")
    
    # TEST 5 : DES example
    
    txt = c(
        'DAE=A(1)',
        'DSIZE=A(2)',
        'PD=1-DAE/(AE50+DAE)',
        'DADT(2) = (RIN*PD - DSIZE*KOVER)*DSIZE',
        'DADT(1) = -KPD * A(1)'
    )
    test5 <- unname(sapply(txt, function(x) {
        nonmem.parameters.specialarray.to.r.array(parse1(x)) }))
    checkEquals(test5[[4]][[2]],  quote( DADT[2] ), msg = "check elements of txt for des example")
    checkEquals(test5[[4]][[3]],  quote( (RIN*PD - DSIZE*KOVER)*DSIZE) )
    checkEquals(test5[[5]][[2]],  quote( DADT[1] ))
    checkEquals(test5[[5]][[3]],  quote(  -KPD * A[1] ))
}


test.resolve.to.mapping.chain <- function()
{

    # TEST 1 : string

    txt <- "TVCL = THETA(1) * (WT / 70)**0.75"
    test1 <- resolve.to.mapping.chain(txt)
    checkEquals(test1[[1]]$relations[[2]][[3]][[3]], quote((WT/70)^0.75), msg = "string")
    

}
