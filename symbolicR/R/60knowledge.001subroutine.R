# recored knowledge of NONMEM Subroutine here

#' just a handy macro
Matrix = function(m,...){
    m0 = substitute(m)
    # force to eval it
    if (is.symbol(m0)) m0 = m
    symbolic.CMatrix(tuple.to.list(m0),...)
}

NONMEM.SUBROUTINE.DATABASE = list(
ADVAN1=list(
    compartments = c('CENTRAL','Output'),
    parameters=list(
        basic=list(
            # V might be not used for ODE(kinetic system), however, when we need to compute
            # Concentration, we have to get V, usually using S2=V to cause F to be computed as concentration
            default=list( names = Matrix(TUPLE(K)) ),
            TRANS2=list( names = Matrix(TUPLE(CL, V)),
                         transforms = nonmem.eqns.to.r.eqns(c(
                                        ' K = CL/V ')))
        ),
        scaling=list(names=Matrix(TUPLE(S1,S2)),
                        alias=nonmem.eqns.to.r.eqns(c(
                            'S0 = S2',
                            'SC = S1'))),
        lag    = list(names=Matrix(TUPLE(ALAG1))),
        bioavl = list(names=Matrix(TUPLE(F1))),
        rate   = list(names=Matrix(TUPLE(R1))),
        duration=list(names=Matrix(TUPLE(D1))),
        fraction=list(names=Matrix(TUPLE(F0)),
                      alias=nonmem.eqns.to.r.eqns(c(
                            'FO = F0',
                            'F2 = F0'))))
),
ADVAN2=list(
    compartments = c('DEPOT','CENTRAL','Output'),
    parameters=list(
        basic=list(
            # V might be not used for ODE(kinetic system), however, when we need to compute
            # Concentration, we have to get V, usually using S2=V to cause F to be computed as concentration
            default=list( names = Matrix(TUPLE(K, KA)) ),
            TRANS1=list( names = Matrix(TUPLE(K, KA)) ),
            TRANS2=list( names = Matrix(TUPLE(CL, V, KA)),
                         transforms = nonmem.eqns.to.r.eqns(c(
                                        ' K = CL/V ',
                                        ' KA = KA'
                         ))
            )
        ),
        scaling=list(names=Matrix(TUPLE(S1,S2,S3)),
                        alias=nonmem.eqns.to.r.eqns(c(
                            'S0 = S3',
                            'SC = S2'))),
        lag    = list(names=Matrix(TUPLE(ALAG1,ALAG2))),
        bioavl = list(names=Matrix(TUPLE(F1,F2))),
        rate   = list(names=Matrix(TUPLE(R1,R2))),
        duration=list(names=Matrix(TUPLE(D1,D2))),
        fraction=list(names=Matrix(TUPLE(F0)),
                      alias=nonmem.eqns.to.r.eqns(c(
                            'FO = F0',
                            'F3 = F0'))))
),
ADVAN3=list(
    compartments = c('CENTRAL','Periph','Output'),
    parameters=list(
        basic=list(
            default=list( names = Matrix(TUPLE(K, K12, K21)) ),
            TRANS1=list( names = Matrix(TUPLE(K, K12, K21)) ),
            TRANS3=list( names = Matrix(TUPLE(CL, V, Q, VSS)),
                         transforms = nonmem.eqns.to.r.eqns(c(
                                        ' K = CL/V ',
                                        ' K12 = Q/V ',
                                        ' K21 = Q/(VSS-V)'))
            ),
            TRANS4=list(   names = Matrix(TUPLE(CL, V1, Q, V2) ),
                          transforms = nonmem.eqns.to.r.eqns(c( 
                                        ' K = CL/V1 ',
                                        ' K12 = Q /V1 ',
                                        ' K21 = Q /V2 '))),
            TRANS5=list(   names = Matrix(TUPLE(AOB, ALPHA, BETA) ),
                          transforms = nonmem.eqns.to.r.eqns(c( 
                                        ' K21 = (AOB * BETA + ALPHA)/(AOB + 1) ',
                                        ' K = (ALPHA * BETA) / K21 ',
                                        ' K12 = ALPHA + BETA - K21 -K '))),
            TRANS6=list(   names = Matrix(TUPLE(ALPHA, BETA, K21) ),
                          transforms = nonmem.eqns.to.r.eqns(c( 
                                        ' K = (ALPHA * BETA) / K21 ',
                                        ' K12 = ALPHA + BETA - K21 -K ')))
        ),
        scaling=list(names=Matrix(TUPLE(S1,S2,S3)),
                        alias=nonmem.eqns.to.r.eqns(c(
                            'S0 = S3',
                            'SC = S1'))),
        lag    = list(names=Matrix(TUPLE(ALAG1,ALAG2))),
        bioavl = list(names=Matrix(TUPLE(F1,F2))),
        rate   = list(names=Matrix(TUPLE(R1,R2))),
        duration=list(names=Matrix(TUPLE(D1,D2))),
        fraction=list(names=Matrix(TUPLE(F0)),
                      alias=nonmem.eqns.to.r.eqns(c(
                            'FO = F0',
                            'F3 = F0'))))),
ADVAN4=list(
    compartments = c('DEPOT','CENTRAL','Periph','Output'),
    parameters=list(
        basic=list(
            # V might be not used for ODE(kinetic system), however, when we need to compute
            # Concentration, we have to get V, usually using S2=V to cause F to be computed as concentration
            default=list( names = Matrix(TUPLE(K, K23, K32, KA)) ),
            TRANS1=list( names = Matrix(TUPLE(K, K23, K32, KA)) ),
            TRANS3=list( names = Matrix(TUPLE(CL, V, Q, VSS, KA)),
                         transforms = nonmem.eqns.to.r.eqns(c(
                                        ' K = CL/V ',
                                        ' K23 = Q/V ',
                                        ' K32 = Q/(VSS-V)',
                                        ' KA = KA'
                         ))
            ),
            TRANS4=list(   names = Matrix(TUPLE(CL, V2, Q, V3, KA) ),
                          transforms = nonmem.eqns.to.r.eqns(c( 
                                        ' K = CL/V2 ',
                                        ' K23 = Q /V2 ',
                                        ' K32 = Q /V3 '))),
            TRANS5=list(   names = Matrix(TUPLE(AOB, ALPHA, BETA, KA) ),
                          transforms = nonmem.eqns.to.r.eqns(c( 
                                        ' K32 = (AOB * BETA + ALPHA)/(AOB + 1) ',
                                        ' K = (ALPHA * BETA) / K32 ',
                                        ' K23 = ALPHA + BETA - K32 -K ',
                                        ' KA = KA ')))
        ),
        scaling=list(names=Matrix(TUPLE(S1,S2,S3,S4)),
                        alias=nonmem.eqns.to.r.eqns(c(
                            'S0 = S4',
                            'SC = S2'))),
        lag    = list(names=Matrix(TUPLE(ALAG1,ALAG2,ALAG3))),
        bioavl = list(names=Matrix(TUPLE(F1,F2,F3))),
        rate   = list(names=Matrix(TUPLE(R1,R2,R3))),
        duration=list(names=Matrix(TUPLE(D1,D2,D3))),
        fraction=list(names=Matrix(TUPLE(F0)),
                      alias=nonmem.eqns.to.r.eqns(c(
                            'FO = F0',
                            'F4 = F0'))))
))

#' from equations, find possible PK parameters
#' The reason why we need \code{eqns} here is because alias
#' we actually do not want all possible ones
#' @param eqns equations
#' @param subs PREDPP subroutine name
#' @param trans translation name used
#' @param try_abbr_pk_names should try match ALPH to ALPHA, e.t.c
#' @return list of found PK parameters
find.socp.symbols.by.subroutine = function(eqns, subs='ADVAN4', trans='', try_abbr_pk_names=TRUE) {
    if (!subs %in% names(NONMEM.SUBROUTINE.DATABASE)) {
        stop(sprintf('Subroutine %s not found', subs))
    }
    subr = NONMEM.SUBROUTINE.DATABASE[[ subs ]]
    if (is.null(trans) || trans=='') trans='default'
    trans = sub('^TRAN.*(\\d+)', 'TRANS\\1', trans)
    if (!trans %in% names(subr$parameters$basic) ) {
        if (trans=='TRANS1') trans='default'
        else stop(sprintf('The trans routine %s is not available for %s', trans, subs))
    }
    gg = create.abstractGraph.from.equations(eqns)
    vars = gg$vertex.table
    re.var = list()
    re.type = character(0)
    gg = create.abstractGraph.from.equations(subr$parameters$scaling$alias)
    scaling.alias = gg$vertex.table[ gg$lhs.ind ]
    gg = create.abstractGraph.from.equations(subr$parameters$fraction$alias)
    fraction.alias = gg$vertex.table[ gg$lhs.ind ]
    ##
    if (try_abbr_pk_names) {
        vars.ch = sapply(vars, function(x) paste(deparse(x),collapse=''))
        abbr.str = function(vv) {
            substring(sapply(vv, function(x) paste(deparse(x),collapse='')),1,4)
        }
        # only those LONG names as ALPHA has problem, short ones as scaling and fraction, no problem
        # and ofcouse, ALAG1 can not be ALAG
        basic.para.ch = abbr.str(subr$parameters$basic[[ trans ]]$names)
    }
    ##
    for(i0 in seq_along(vars)) {
    # current, not support ALPH as ALPHA yet, should be exactly same
    # however, we may deal with this in socp.filter code
        i = vars[[ i0 ]]
        if (!is.na(match.symbolic(i, subr$parameters$basic[[ trans ]]$names)) ||
            (try_abbr_pk_names && (vars.ch[ i0 ] %in% basic.para.ch))) {
            re.var = union.symbolic(re.var, i)
            re.type[ length(re.var) ] = 'basic'
            next()
        }
        if (!is.na(match.symbolic(i, subr$parameters$scaling$names))) {
            re.var = union.symbolic(re.var, i)
            re.type[ length(re.var) ] = 'scaling'
            next()
        }
        if (!is.na(match.symbolic(i, scaling.alias))) {
            re.var = union.symbolic(re.var, i)
            re.type[ length(re.var) ] = 'scaling'
            next()
        }
        if (!is.na(match.symbolic(i, subr$parameters$lag$names))) {
            re.var = union.symbolic(re.var, i)
            re.type[ length(re.var) ] = 'lag'
            next()
        }
        if (!is.na(match.symbolic(i, subr$parameters$bioavl$names))) {
            re.var = union.symbolic(re.var, i)
            re.type[ length(re.var) ] = 'bioavl'
            next()
        }
        if (!is.na(match.symbolic(i, subr$parameters$rate$names))) {
            re.var = union.symbolic(re.var, i)
            re.type[ length(re.var) ] = 'rate'
            next()
        }
        if (!is.na(match.symbolic(i, subr$parameters$duration$names))) {
            re.var = union.symbolic(re.var, i)
            re.type[ length(re.var) ] = 'duration'
            next()
        }
        if (!is.na(match.symbolic(i, subr$parameters$fraction$names))) {
            re.var = union.symbolic(re.var, i)
            re.type[ length(re.var) ] = 'fraction'
            next()
        }
        if (!is.na(match.symbolic(i, fraction.alias))) {
            re.var = union.symbolic(re.var, i)
            re.type[ length(re.var) ] = 'fraction'
            next()
        }
    }
    basic.ind = which(re.type=='basic')
    basic.expected.length = length(subr$parameters$basic[[ trans ]]$names)
    if (length(basic.ind) != basic.expected.length) {
        err.type = if (length(basic.ind)<basic.expected.length) 'Less' else 'More'
        warning(sprintf('%s basic PK parameters found than expected, Found [%s], but Expected [%s]',
          err.type,
          paste(sapply(re.var[ basic.ind ], function(x) paste(deparse(x),collapse='')),collapse=','),
          paste(sapply(subr$parameters$basic[[ trans ]]$names, 
            function(x) paste(deparse(x),collapse='')),collapse=',')))
    }
    re.var = symbolic.CMatrix(re.var)
    attr(re.var,'type')=re.type
    re.var
}
