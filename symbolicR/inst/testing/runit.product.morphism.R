
test.product.morphism = function(){
    e = quote( ifelse(RACE == 3, 1 + THETA[5], ifelse(RACE == 2, 1 + THETA[4], 
        ifelse(RACE == 1, 1, NA))) * ifelse(SEX == 2, 1 + THETA[6], 
        ifelse(SEX == 1, 1, NA)) * ifelse(HCTZ == 0, 1 + THETA[7], 
        ifelse(HCTZ == 1, 1, NA)) * THETA[1])


    mor = pattern.product.expression.containing.ifelse(e)
    expected = quote(
        ifelse(RACE == 3, 1 + THETA[2], ifelse(RACE == 2, 1 + THETA[1], 
            ifelse(RACE == 1, 1, NA))) * ifelse(SEX == 2, 1 + THETA[3], 
            ifelse(SEX == 1, 1, NA)) * ifelse(HCTZ == 0, 1 + THETA[4], 
            ifelse(HCTZ == 1, 1, NA)) * THETA[5])

    q1 = instantiate.morphism(mor)
    checkEquals(q1, expected)
}

test.sum.morphism = function(){

    e = quote( ifelse(RACE == 3, 1 + THETA[5], ifelse(RACE == 2, 1 + THETA[4], 
        ifelse(RACE == 1, 1, NA))) + ifelse(SEX == 2, 1 + THETA[6], 
        ifelse(SEX == 1, 1, NA)) + ifelse(HCTZ == 0, 1 + THETA[7], 
        ifelse(HCTZ == 1, 1, NA)) + THETA[1])
    mor = pattern.affine.function.theta.or.ifelse(e)
    q1 = instantiate.morphism(mor)
    checkEquals(q1, quote(
                ifelse(RACE == 3, 1 + THETA[2], ifelse(RACE == 2, 1 + THETA[1], 
                    ifelse(RACE == 1, 1, NA))) + ifelse(SEX == 2, 1 + THETA[3], 
                    ifelse(SEX == 1, 1, NA)) + ifelse(HCTZ == 0, 1 + THETA[4], 
                    ifelse(HCTZ == 1, 1, NA)) + THETA[5] ))


    e2 = quote( ifelse(RACE == 3, 1 + THETA[5], ifelse(RACE == 2, 1 + THETA[4], 
        ifelse(RACE == 1, 1, NA))) * ifelse(SEX == 2, 1 + THETA[6], 
        ifelse(SEX == 1, 1, NA)) + ifelse(HCTZ == 0, 1 + THETA[7], 
        ifelse(HCTZ == 1, 1, NA)) * THETA[1])

    mor2 = pattern.affine.function.theta.or.ifelse(e2)
    q2 = instantiate.morphism(mor2)
}

test.concrete.morphism = function(){
    m1 = morphism.concrete.new('Constant',10)
    checkEquals(instantiate.morphism(m1),10)

    m2 = morphism.concrete.new('Constant',1:2)
    checkEquals(instantiate.morphism(m2),list(1,2))

    m3 = morphism.concrete.new('identicalmapping',dim=1)
    checkEquals(instantiate.morphism(m3), quote(THETA[1]))

    m4 = morphism.concrete.new('identicalmapping',dim=2)
    checkEquals(instantiate.morphism(m4), list(quote(THETA[1]),quote(THETA[2])))

    v0 = list(
        quote( 1+THETA[1]),
        quote( THETA[2] + y * THETA[3]),
        quote( 32 )
    )
    m5 = morphism.concrete.new('AffineMapping', 
        pattern.affine.mapping.of.theta(v0))

    m6 = morphism.concrete.new('Product', dim=3)
    m6.m5 = morphism.composition.2(m6,m5)
    checkEquals( instantiate.morphism(m6.m5, as.list(c(1,2,3))),
        quote( 64 * ( 2 + 3 * y)))

    checkEquals( instantiate.morphism(m6.m5, as.list(c(-1,2,3))), 0)

    m7 = morphism.concrete.new('Sum', dim=3)
    checkEquals(instantiate.morphism(m7, as.list(1:3)), 6)
    
    
}
