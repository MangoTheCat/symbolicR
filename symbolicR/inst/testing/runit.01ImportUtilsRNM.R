# git revision: $Rev: $
# Date of last change: $LastChangedDate: 28/09/2012 $
# Last changed by: $LastChangedBy: ccampbell $
# 
# Original author: ccampbell
# Copyright Mango Solutions, Chippenham, UK
###############################################################################


test.neighbourhood.of.v <- function()
{
    # TEST 1 : Simple Graph left hand neighbour
    
    graph <- list(edges = list(c("B", "A"), c("C", "B")))
    test1 <- neighbourhood.of.v(G = graph, v = "C")
    checkEquals(test1, list(v.income = character(0), v.outcome = "B"), msg = "left hand neighbour")
    
    # TEST 2 : Simple Graph both neighbours
   
    test2 <- neighbourhood.of.v(G = graph, v = "B")
    checkEquals(test2, list(v.income = "C", v.outcome = "A"), msg = "both neighbours")

    pks <- c(
        'TVCL = THETA(1) * WT/70',
        'TVV = THETA(2)',
        'TVKA = THETA(3)',
        'TVLG = THETA(4)',
        'CL = TVCL * EXP(ETA(1))',
        'V = TVV * EXP(ETA(2))',
        'KA = TVKA * EXP(ETA(3))',
        'ALAG1 = TVLG * EXP(ETA(4))',
        'S2 = V')
    
    # TEST 3 : PK Example one vertex 2 left hand neighbours
    
    pkGr <- create.graph.from.equations(lines = pks)$graph
    
    test3 <- neighbourhood.of.v(G = pkGr, v = "KA")
    checkEquals(test3, list(v.income = c("TVKA", "ETA[]"), v.outcome = character(0)), 
        msg = "2 left hand neighbours")
    
    # TEST 4 : PK Example five vertices
    
    test4 <- neighbourhood.of.v(G = pkGr, v = c('CL','V','KA','ALAG1','S2'))
#    # TODO: Check prefered behaviour
#    # what is expected here?
#@answer:
#    should have coming vertices as c('TVCL','TVV','TVKA','TVLG','ETA[]') (the order might be different)
#    because these are the vertex pointing directly into any element of the set {CL,V,KA,ALAG1,S2} and not in the set.
#
    checkTrue(setequal(test4$v.income, c('TVCL','ETA[]','TVV','TVKA','TVLG')))
    checkEquals(test4$v.outcome,character(0))
#
#    checkEquals(test4, list(v.income = c("TVKA", "ETA[]"), v.outcome = character(0)), 
#        msg = "2 left hand neighbours")
}


test.restriction.morphism.of.graph.at.v <- function()
{
    # graph
    
    # TEST 1
    
    pks0 = c(
        'TVCL=THETA(1)*WT/70',
        'TVV =THETA(2)',
        'TVKA=THETA(3)',
        'TVLG=THETA(4)',
        'CL=TVCL*EXP(ETA(1))',
        'V =TVV*EXP(ETA(2))',
        'KA=TVKA*EXP(ETA(3))',
        'ALAG1=TVLG*EXP(ETA(4))',
        'S2=V')

    pkGr <- create.graph.from.equations(pks0)$graph
    
    neighbV <- neighbourhood.of.v(pkGr, 'ALAG1')$v.income
    
    test1 <- restriction.morphism.of.graph.at.v(pkGr, c('ALAG1', neighbV))
    
    checkEquals(test1, c(12, 13), msg = "TEST 1: warf_base")
    
}
