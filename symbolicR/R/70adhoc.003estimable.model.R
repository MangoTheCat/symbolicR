#' Ensure return value of subsetting from a \code{triple.matrix} to be non-\code{NA}
#'
#' If \code{ind0}, the index usr provide is \code{NA}, the return is just the \code{default.val.for.col} default values. 
#' e.g. For column one (lower bound) is minus infinity and upper bound positive infinity and the initial estimation to be a small positive number(NOT zero! NONMEM will complain)
#'
#' If \code{triple.matrix} itself has \code{NA}s, and \code{remove.original.nas} is \code{TRUE}, 
#' those \code{NA}s are also set to \code{default.val.for.col} by their column position.
#'
#' @param ind0 original indices
#' @param triple.matrix typically a Theta in nmmod
#' @param columns select which columns
#' @param default.val.for.col default value for specified column for non-existing indices
#' @param remove.original.nas if there are NAs in the input matrix, we set those to default values
#' @return selected array
get.by.ind.triple.matrix = function(ind0, triple.matrix, columns=1:3, default.val.for.col=c(-Inf,0.0001,Inf), remove.original.nas=TRUE) {
    ind = match(ind0, 1:NROW(triple.matrix))
    ind.na = which(is.na(ind))
    if (length(ind.na)>0) {
        warning(sprintf('These indices [%s] are outside Matrix size.',paste(ind0[ind.na],collapse=',')))
    }
    re = array(0, dim=c(length(ind), length(columns)))
    for(k in seq_along(columns)) {
        re[,k]       = triple.matrix[ind,columns[k]]
        re[ind.na,k] = default.val.for.col[columns[k]]
        # protect against input triple.matrix's NA's
        if (remove.original.nas) {
            ori.ind.na = which(is.na(re[,k]))
            if (length(ori.ind.na)>0) {
                re[ori.ind.na, k ] = default.val.for.col[ columns[k] ]
            }
        }
    }
    if (length(columns)==1) re = as.vector(re)
    re
}

#' From nmmod and nmdata (object returned by \code{RNMImport} package), generate the \code{estimable.model} (Model Spec)
#'
#' This is the main API in this package behaving as a user level API.
#' It's working as an extension of the \code{RNMImport} pacakge. Mango's \code{RNMImport} deal with control file, or NONMEM run,
#' but the PK block(code) and ERROR block(code) in the control file are stored faithfully as plain strings. 
#' This function tries to do deeper information extraction from the return object of \code{RNMImport}.
#' This second phase analysing depends on the following assumptions.
#' Quite a few NONMEM control file are written using the pre-defined PREDPP compartment subroutines, and those subroutines have
#' definite input PK parameters, like \code{CL}, \code{V} and e.t.c. The variant part are how these parameters are constituted,
#' e.g. you may simplely use one parameter for typical value and additional populational random effect. But also you can model 
#' it with multiple parameters and also invole other covariables. 
#' So, even one simple kinetic model as ADVAN2 can derive tens of different concrete models. But there are something in common,
#' and something can be re-used or managed by computer, i.e. the constitution rules. We use \code{estimable.unit} and \code{random.observation} for them.
#' For the math, the \code{estimable.unit} can be a \code{morphism} or \code{random.disturbation}.
#' And based on above, the two models with same kinetic system can be compared, even they have different numbers of \code{THETA}s.
#'
#' Of course, base on above assumptions, we do not support those NONMEM code written in \code{PRED} and those using a generic ODE kinetic system like
#' ADVAN6, e.t.c
#'
#' @note Currently, this function resolve duplicated use of parameters by introducing new variables
#'
#' \code{
#'  V = THETA[1]
#'  CL = THETA[1] + THETA[2] }
#' Then it will think it as
#' \code{
#'  THETAHAT01 = THETA[1]
#'  V = THETAHAT01
#'  CL = THETAHAT01 + THETA[2]}
#' 
#' Also, we assume the \code{nmmod} should have both PK and ERROR part.
#' If the control file, or say, the \code{nmmod} only has PRED code, then it hard to determine what the model
#' means. Now, we stop if only PRED appears.
#'
#' @param nmmod the nonmem model return by \code{importNmMod}
#' @param modData the Input Data returned by \code{RNMImport:::importModelData}
#' @param keep_mu_i_s should keep MU_i in NONMEM statements
#' @param try_abbr_pk_names if try to catch PK parameters using abbreviation like ALPHA -> ALPH
#' @return specification
create.model.spec.from.nmmod.and.nmdata = function(nmmod,modData,keep_mu_i_s=TRUE,try_abbr_pk_names=TRUE) {
    # If only PRED, stop
    if (!is.null(nmmod$PRED))  {
        stop('Analysing nmmod using PRED code is not supported yet.')
    }
    #
    ## INPUT VARIABLES in nmmod$Input
    INPUT.VARIABLES = lapply(toupper(nmmod$Input[,'nmName']), as.symbol)
    ## TABLE VARIABLES nmmod$Tables
    TABLE.VARIABLES = tuple.to.list(
            parse1(sprintf('TUPLE(%s)',paste(nmmod$Table[, 'Columns'], collapse=','))))
    # socp should be those not being a covariate, but need in the output tables
    ind.in.table.not.in.input = is.na(match.symbolic(TABLE.VARIABLES, INPUT.VARIABLES))
    ############################################
    SOCP.VARIABLES.FOR.ERROR.BLOCK = 
        union.symbolic(
                        tuple.to.list(quote(TUPLE(Y,IPRED,IPRE, IRES, IWRES))),
                        TABLE.VARIABLES[ ind.in.table.not.in.input ])
    ## find all appearing ADVAN defined PK parameters
    pks.eqns = nonmem.eqns.to.r.eqns(nmmod$PK)
    err.eqns = nonmem.eqns.to.r.eqns(nmmod$Error)
    ## Try to eliminate possible duplication
    duplication = eliminate.duplicated.appearence.of.rhs.parameters(c(pks.eqns, err.eqns))
    if (any( duplication$rep.var > 0 ) ) {
        N.pks.eqns = length(pks.eqns)
        N.err.eqns = length(err.eqns)
        if (N.pks.eqns>0) {
            pks.eqns = c(
                duplication$intro.eqns$THETA,
                duplication$intro.eqns$ETA,
                duplication$new.eqns[ 1:N.pks.eqns ])
        }
        if (N.err.eqns>0) {
            err.eqns = c(
                duplication$intro.eqns$EPS,
                duplication$new.eqns[ N.pks.eqns + (1:N.err.eqns) ])
        }
    }
    ##
    subs = nmmod$Subroutine
    if (is.null(subs)) {
        stop('No subroutine statement found in nmmod, seems it is not described by PREDPP kinetic system
        Currently, the directly closed formula using $PRED can not be analysed by this function.  ')
    }
    if (length(subs)>1) trans = subs[2] else trans = ''
    advan.pk.para = find.socp.symbols.by.subroutine(pks.eqns, subs=subs[1], trans=trans, 
        try_abbr_pk_names=try_abbr_pk_names)
    ## exclude from error.symbols as input from PK block
    known.excluded.table = union.symbolic(c(SOCP.VARIABLES.FOR.ERROR.BLOCK, 
                                quote(F), INPUT.VARIABLES),
                             advan.pk.para)
    ##
    obs.socp.filter = function(gg=NULL,variables=NULL,eqns=NULL){
        which(!is.na(match.symbolic(variables, SOCP.VARIABLES.FOR.ERROR.BLOCK)))
    }
    ##
    obs.disturbation.patterns = list()
    ## multiple dv, usually indicated by DVID
    dvid.var = NULL
    if (!missing(modData) && 'DVID' %in% toupper(nmmod$Input[,'nmName']) &&
        'DVID' %in% names(modData)) {
        dvid.var = quote(DVID)
        dvid.lvs= sort(unique(modData$DVID))
    }
    if (!is.null(dvid.var)) {
        obs.disturbation.patterns[[ length(obs.disturbation.patterns) + 1 ]] = 
        list(pat=function(defn){
            if(!is.null(pat<-pattern.multiple.dimensional.by.factor(
                                defn$final.e,
                                varname=dvid.var,
                                var.levels=dvid.lvs ))){
                pat$type='multiple.dependent.variable'
                return(pat)
            }
            NULL
        },
        action=function(defn,pat){
                list(critical.varname=quote(Y), type='directly.related.to.nonmem.parameters', defn=pat)
        })
    }
    ## general case
    obs.disturbation.patterns[[ length(obs.disturbation.patterns) + 1 ]]  = list(pat=function(defn) {
            if (symbolic.grep( quote( EPS['?'(epsind)] ), defn$final.e) ) {
                return(list(type='random.observation',defn=defn$final.e))
            }
            NULL
        },
        action=function(defn,pat){
                list(critical.varname=quote(Y), type='directly.related.to.nonmem.parameters', defn=pat)
        })
    ##
    L = analyse.0.0(err.eqns, 
                    socp.filter=obs.socp.filter, 
                    extra.rules=obs.disturbation.patterns, 
                    skip.patterns=TRUE,
                    ignore.non.atomic.coefficient=TRUE
                    )
    if (length(L$other.rels)>0) {
        print(L$other.rels)
        stop('More complex structure in ERROR BLOCK as above detected, stopped')
    }
    # found the starting symbols:
    # should be solely a rhs,
    # should not be a nonmem keyword
    exclude.symbols = function(ll){
        ll = ll[ is.na(match.symbolic(ll, known.excluded.table)) ]
        ll = ll[!sapply(ll, function(x) symbolic.match(quote(A['?a'(ii)]), x))]
        ll = ll[!sapply(ll, function(x) symbolic.match(quote(THETA['?a'(ii)]), x))]
        ll = ll[!sapply(ll, function(x) symbolic.match(quote(ETA['?a'(ii)]), x))]
        ll = ll[!sapply(ll, function(x) symbolic.match(quote(EPS['?a'(ii)]), x))]
        ll
    }
    #
    obsObjs = list()
    error.rhs.symbols = list()
    push.error.rhs.symbols = function(s){
        rhs = exclude.symbols(extract.vertices(s)$rhs.symbol)
        error.rhs.symbols <<- union.symbolic(error.rhs.symbols, rhs)
    }
    symbolic.names = list()
    for(comp in L$analysed.defns) {
        for(varn in comp) {
            foo = varn$defn
            lhs.name = asCharacterSymbolCall(varn$critical.varname)
            symbolic.names[[ length(symbolic.names) + 1 ]] = varn$critical.varname
            if (varn$type %in% c('related.on.previous.variables', 'unknown')) {
                push.error.rhs.symbols(foo)
                obsObjs[[ lhs.name ]] = foo
                next()
            }
            # varn$type should be directly.related.to.nonmem.parameters
            if (foo$type=='random.observation') {
                push.error.rhs.symbols(foo$defn)
                foo = random.disturbation.new.eps(foo$defn)
                if (!is.null(epsind<-attr(foo,'backref'))) {
                    foo$distribution = update.distribution.covariance(
                        foo$distribution, 
                        nmmod$Sigma[ epsind, epsind ,  drop=F])
                }
                obsObjs[[ lhs.name ]] = foo 
                next()
            }
            if (foo$type=='multiple.dependent.variable') {
                for(foo1 in foo$levels) {
                    push.error.rhs.symbols(foo1)
                }
                foo = multiple.dependent.variable.to.disturbation(foo)
                if (!is.null(epsind<-attr(foo,'backref'))) {
                    foo$distribution = update.distribution.covariance(
                        foo$distribution, 
                        nmmod$Sigma[ epsind, epsind ,  drop=F])
                }
                obsObjs[[ lhs.name ]] = foo
                next()
            }
            # THETA or ETA appearing in $ERROR block, quite unusall because it will be computed everytime
            # on an observation record, should be put to PK block
            # however we do met this situation
            if (foo$type %in% c('fixed','normal','log normal')) {
                mean.exp = explain.expression.as.morphism(varn$defn$mean)
                if (is.null(mean.exp) || mean.exp$defn$type=='Constant') mean.exp=varn$defn$mean
                foo = mean.exp
                if (is(foo, 'morphism')) {
                # update domain of morphisms if possible
                    theta.ind = backref.parameter.of.morphism(foo)
                    if (!is.null(theta.ind)) {
                        foo$domain = as.morphism.domain.matrix(get.by.ind.triple.matrix(theta.ind,nmmod$Theta))
                        foo = estimable.Unit.new(foo)
                        attr(foo,'initial.estimation') = get.by.ind.triple.matrix(theta.ind,nmmod$Theta,2)
                    }
                }
                if (!(varn$defn$type=='fixed' || equal.symbol.light(varn$defn$var, 0))) {
                    if (is(foo, 'name') || is(foo,'call') || is.atomic(foo)) {
                    # protect 'name' or 'call' as Constant morphism
                    # before we promote it to an estimable.unit for recording random effect
                        foo = morphism.concrete.new('Constant',foo)
                    }
                    attr(foo,'distribution') = 
                        list(type=varn$defn$type, Var=eval.nonmem.parameter(varn$defn$var, md = nmmod))
                    foo = estimable.Unit.new(foo)
                }
                # verify log normal problem, might be  THETA[i] * exp(ETA[i]) but THETA[i] <0
                foo = estimable.Unit.correct.lognormal(foo)
                obsObjs[[ lhs.name ]] = foo
                next()
            }
            stop(sprintf('Unknow definition type %s ', foo$type))
        }
    }
    attr(obsObjs,'symbolic.names') = symbolic.CMatrix(symbolic.names,nrow=1)
    ##################################
    #### Now begin to look at PKs ####
    ## some table output are not computed in ERROR but PK, should count those
    table.variables.computed.in.pk = TABLE.VARIABLES[ ind.in.table.not.in.input ]
    table.variables.computed.in.pk.ind = is.na(
        match(sapply(table.variables.computed.in.pk, function(x) paste(deparse(x),collapse='')),
        # found those not computed in obsObjs
              names(obsObjs)))
    table.variables.computed.in.pk = table.variables.computed.in.pk[ table.variables.computed.in.pk.ind ]
    socp.pks = union.symbolic(error.rhs.symbols, advan.pk.para)
    socp.pks = union.symbolic(socp.pks, table.variables.computed.in.pk)
    pks.socp.filter = function(gg=NULL,variables=NULL,eqns=NULL){
        ind0 = which(!is.na(match.symbolic(variables, socp.pks)))
    # MU_i should has special meaning in NONMEM, keep them
        if (keep_mu_i_s) {
            v.ch = sapply(variables, function(x) as.character(x)[1])
            ind1 = grep('^MU_\\d+$',v.ch)
            ind0 = union(ind0, ind1)
        }
        ind0
    }
    L1 = analyse.0.0(pks.eqns, 
                     socp.filter=pks.socp.filter,
                     ignore.non.atomic.coefficient=TRUE
                     )
    if (length(L1$other.rels)>0) {
        extra.relations  =  paste(capture.output(print(L1$other.rels)),collapse='')
        warning(sprintf('%s\n%s\nThese relations are ignored in further analysis.',
                'Some extra structure, which is not a PK parameter or 
                 used as TABLE output or used in ERROR block, 
                 is detected in PK BLOCK as following:', extra.relations))
    }
    ## pkObjs
    ## A note about covariance:
    ## covariance between ETA's can't be noticed by the following program if they are not in same expression:
    ## e.g.
    ##  K =  Q*MU_1*ETA[1]  + (1-Q)*MU_2*ETA[3]
    ##  S1 = Q*MU_3*ETA[2]  + (1-Q)*MU_4*ETA[4]
    ## and covariance matrix is a block matrix as following :
    ##  [ 0.04 0.01             ]
    ##  [ 0.01 0.027            ]
    ##  [             0.05 0.01 ]
    ##  [             0.01 0.06 ]
    ## Then, the cov(K,S1) = Q * MU_1*MU_3*cov(ETA[1],ETA[2]) + (1-Q)*MU_2*MU_4*cov(ETA[3],ETA[4])
    ##                     = Q * MU_1*MU_3*0.01 + (1-Q)*MU_2*MU_4*0.01
    ## But the following program can't not recognize this, a basic assumption over this program is that
    ##  every estimable.Unit should be independent with each other
    ## We fix this by check the Omega after creating pkObj, and store the cross estimable unit 
    ##  covariances explicitly there.
    pkObjs = list()
    symbolic.names = list()
    for(comp in L1$analysed.defns){
        for(varn in comp){
            foo = varn$defn
            lhs.name = asCharacterSymbolCall(varn$critical.varname)
            symbolic.names[[ length(symbolic.names) + 1 ]] = varn$critical.varname
            if (varn$type == 'related.on.previous.variables') {
                pkObjs[[ lhs.name ]] = foo
                next()
            }
            if (varn$type == 'unknown'){
                foo = random.disturbation.new.from.theta.eta(foo)
                if (N.fixed.effects(foo) > 0 || N.random.effects(foo)>0) {
                # any non-trivial random disturbation
                    if(!is.null(backref<-attr(foo,'backref'))){
                        if (N.fixed.effects(foo)>0) {
                            theta.ind = backref[foo$FIXED]
                            foo = random.disturbation.update.fixed.domain(foo,
                                as.morphism.domain.matrix(get.by.ind.triple.matrix(theta.ind, nmmod$Theta)))
                            foo = estimable.Unit.new(foo)
                            attr(foo,'initial.estimation') = get.by.ind.triple.matrix(theta.ind, nmmod$Theta, 2)
                        }
                        if (N.random.effects(foo)>0) {
                            eta.ind = backref[foo$RANDOM]
                            # promote the disturbation to a estimable.Unit.new
                            foo = estimable.Unit.new(foo)
                            foo$distribution = update.distribution.covariance(
                                foo$distribution,
                                nmmod$Omega[ eta.ind, eta.ind, drop=F ])
                        }
                    }
                    pkObjs[[ lhs.name ]] = foo
                } else {
                    pkObjs[[ lhs.name ]] = varn$defn
                }
                next()
            }
            # directly.related.to.nonmem.parameters
            mean.exp = explain.expression.as.morphism(varn$defn$mean)
            if (is.null(mean.exp) || mean.exp$defn$type=='Constant') mean.exp=varn$defn$mean
            foo = mean.exp
            if (is(foo, 'morphism')) {
            # update domain of morphisms if possible
                theta.ind = backref.parameter.of.morphism(foo)
                if (!is.null(theta.ind)) {
                    foo$domain = as.morphism.domain.matrix(get.by.ind.triple.matrix(theta.ind, nmmod$Theta))
                    foo = estimable.Unit.new(foo)
                    attr(foo,'initial.estimation') = get.by.ind.triple.matrix(theta.ind, nmmod$Theta, 2)
                }
            }
            if (!(varn$defn$type=='fixed' || equal.symbol.light(varn$defn$var, 0))) {
                if (is(foo, 'name') || is(foo,'call') || is.atomic(foo)) {
                # protect 'name' or 'call' as Constant morphism
                # before we promote it to an estimable.unit for recording random effect
                    foo = morphism.concrete.new('Constant',foo)
                }
                foo = estimable.Unit.new(foo)
                attr(foo,'distribution') = 
                    list(type=varn$defn$type, Var=eval.nonmem.parameter(varn$defn$var, md = nmmod),
                         backref=eval.backref.from.var(varn$defn$var))
            }
            # verify log normal problem, might be  THETA[i] * exp(ETA[i]) but THETA[i] <0
            foo = estimable.Unit.correct.lognormal(foo)
            pkObjs[[ lhs.name ]] = foo
        }
    }
    attr(pkObjs,'symbolic.names') = symbolic.CMatrix(symbolic.names,nrow=1)
    ## Check and repopulate the covariance for Omega part, if there is any
    ## Note here we combine both c(pkObjs, obsObjs) because possible some user may put estimable units in Error BLOCK
    ## Not a good habit, but allowed by NONMEM
    tb = table.random.parts.list.code(c(pkObjs,obsObjs))
    rvs = unique(names(tb$length)[ tb$length!= 0 ])
    n.rvs = length(rvs) 
    PK.COV = list()
    if (n.rvs > 1) {
    # estimable.unit and random.disturbation can catch internal covariance in one unit
    # we only need to add covariance outside the units
        Omega = nmmod$Omega
        for(v1.ind in 1:(n.rvs-1)){
            for(v2.ind in (v1.ind+1):n.rvs){
                v1 = rvs[ v1.ind ]
                v2 = rvs[ v2.ind ]
                for(v1.subind in 1:tb$length[v1]) {
                    for(v2.subind in 1:tb$length[v2]) {
                        ind = c(tb$prestart[v1]+v1.subind, tb$prestart[v2]+v2.subind)
                        bk = tb$backref[ ind ]
                        if (isTRUE(abs(Omega[bk[1],bk[2]]) > .Machine$double.eps)) {
                            PK.COV[[ length(PK.COV) + 1]] = 
                                list(varnames=c(v1,v2),
                                     subind=c(v1.subind,v2.subind),     
                                     val = Omega[ bk[1], bk[2] ])
                        }
                    }
                }
            }
        }
    }
    # Everything is similar for Sigma part, unless the 'length','prestart','lookup','backref' all has '2' as postfix
    rvs = unique(names(tb$length2)[ tb$length2!= 0 ])
    n.rvs = length(rvs) 
    OBS.COV = list()
    if (n.rvs > 1) {
        Sigma = nmmod$Sigma
        for(v1.ind in 1:(n.rvs-1)){
            for(v2.ind in (v1.ind+1):n.rvs){
                v1 = rvs[ v1.ind ]
                v2 = rvs[ v2.ind ]
                for(v1.subind in 1:tb$length2[v1]) {
                    for(v2.subind in 1:tb$length2[v2]) {
                        ind = c(tb$prestart2[v1]+v1.subind, tb$prestart2[v2]+v2.subind)
                        bk = tb$backref2[ ind ]
                        if (isTRUE(abs(Sigma[bk[1],bk[2]]) > .Machine$double.eps)) {
                            OBS.COV[[ length(OBS.COV) + 1]] = 
                                list(varnames=c(v1,v2),
                                     subind=c(v1.subind,v2.subind),     
                                     val = Sigma[ bk[1], bk[2] ])
                        }
                    }
                }
            }
        }
    }
    #------------------------------------------#
    ## start populate the whole object
    # other.para might be indication a PD model
    other.para = error.rhs.symbols[ is.na(match.symbolic(error.rhs.symbols, advan.pk.para)) ]
    if (length(other.para)>0) other.para = symbolic.CMatrix(other.para, nrow=1)
    re = list(
        META.INFO=list(
            PK.PARAMETERS    = symbolic.CMatrix(advan.pk.para,nrow=1),
            OTHER.PARAMETERS = other.para
        ),
        DYNAMICS    = kinetic.system.predpp.new(nmmod$Subroutine),
        PK.PROGRAM  = pkObjs,
        OBS.PROGRAM = obsObjs
    )
    class(re) = 'estimable.model'
    if (length(PK.COV)>0) re$PK.COV = PK.COV
    if (length(OBS.COV)>0) re$OBS.COV = OBS.COV
    re
}

#' For a list of expressions or estimable units or random disturbation,
#' sketch a overall graph of the random parts
#'
#' If any random disturbation is found, 
#' a second \code{length}, \code{prestart}, \code{lookup} and \code{backref} are returned,
#' which are coded as \code{length2}, \code{prestart2}, \code{lookup2} and \code{backref2}.
#'
#' Simple expression or names(symbols) are skipped. Only estimable unit and random disturbations are counted.
#'
#' @param lesu list of estimable.units or random disturbation
#' @return number of random parts, start position in all, from new index to look up name, backref
table.random.parts.list.code = function(lesu){
    nms = names(lesu)
    if (is.null(nms)) nms = as.character(seq_along(lesu))
    ind = integer(0)
    backref = integer(0)
    #
    haseta = F
    haseps = F
    ind0 = integer(0)
    backref0 = integer(0)
    #
    for(nm in nms) {
        eu = lesu[[ nm ]]
        if (is(eu,'estimable.unit')) {
        # THETA or ETA
            haseta = T
            backref = c(backref,
                lookup.random.part.backref.estimable.unit(eu))
            if (is(eu,'random.disturbation')) {
                ind[[ nm ]] = N.random.effects(eu)
                ind0[[nm]] = 0
                next()
            }
            if (is(eu,'morphism') || is(eu,'name') || is(eu,'call') || is.atomic(eu)) {
                if (is.null(attr(eu,'distribution'))) {
                    ind[[nm]]=0
                    ind0[[nm]] = 0
                } else {
                    ind[[ nm ]] = 1
                    ind0[[nm]] = 0
                }
                next()
            }
            stop(sprintf('Do not know how to determine the number of random effects of estimable unit of type: %s',
                paste(class(eu),collapse=',')))
        } else if (is(eu,'random.disturbation')) {
        # POSSIBLE EPS
            haseps = T
            backref0 = c(backref0,
                lookup.random.part.backref.random.disturbation(eu))
            ind0[[ nm ]] = N.random.effects(eu)
            ind[[ nm ]] = 0
            next()
        }
        # ELSE, call or a name
        ind[[ nm ]] = 0
        ind0[[nm ]] = 0
    }
    ind1 = ind[ind!=0]
    ind1 = cumsum(ind1)
    prestart = c(0,ind1[-length(ind1)])
    names(prestart) = names(ind1)
    if (haseps) {
        ind2 = ind0[ind0!=0]
        ind2 = cumsum(ind2)
        prestart2 = c(0, ind2[ - length(ind2) ] )
        names(prestart2) = names(ind2)
    }
    if (haseta) {
        re = list(
            # length of random part of all names
            length = ind, 
            # for those non-zero ones, the prestart part, add to the local index in that estimable unit 
            # to get the over all index
            prestart = prestart,
            # from overall index get the name
            lookup = rep(nms, times=ind),
            # For a give overall index IND :  (IND - prestart[lookup[IND]]) will give the local index
            backref = backref
        )
    } else {
        re = list()
    }
    if (haseps) {
        re$length2 = ind0
        re$prestart2 = prestart2
        re$lookup2 = rep(nms, times=ind0)
        re$backref2 = backref0
    }
    re
}

#' extract index of ETA or EPS 
#' @param e expression
#' @return index
extract.random.index = function(e){
    inds.eta <- inds.eps <- integer(0)
    rules = list(
        list(quote(ETA['?a'(i)]), 
             function(e,d) { 
                inds.eta <<-union(inds.eta, d$i) 
                e} ),
        list(quote(EPS['?a'(i)]), 
             function(e,d) { 
                inds.eps <<-union(inds.eps, d$i) 
                e} )
    )
    walker = symbolic.simplify.gigo(rules)
    walker(e)
    list(eta=inds.eta, eps=inds.eps)
}

#' extract the INDEX of the \code{ETA} entries
#'
#' internal use, Only on \code{ETA}
#' @param e variance expression
#' @return the index
eval.backref.from.var = function(e){
    if (is.atomic(e)) {
        stop('Atomic variance information, backref information is lost.')
    }
    if (symbolic.match(quote(Var(ETA['?a'(II)])), e)) {
        return(II)
    } else if (symbolic.match(quote('?'(coeff)*Var(ETA['?a'(II)])), e)) {
        return(II)
    } 
    stop(sprintf('Variance format for lookup backref not supported for : %s',paste(deparse(e),collapse='')))
}

#' S3 version
print.estimable.model = function(mod){
    indent.n = function(ss,n) {
       if (n==0) return(ss)
       indent = paste(rep(' ',n),collapse='')
       ss = gsub('^', indent, ss)
       gsub('\n', sprintf('\n%s',indent), ss)
    }
    check.meta.information(mod)
    s = print.Matrix0(mod$META.INFO$PK.PARAMETERS)
    s = sprintf('Involved PK parameters: %s',s)
    if (length(mod$META.INFO$OTHER.PARAMETERS)>0) {
        s = sprintf('%s, Other variables(maybe PD related): %s',
            s, print.Matrix0(mod$META.INFO$OTHER.PARAMETERS))
    }
    check.isvalid.basic.pk.parameters(mod)
    s = sprintf('%s\nDynamics Using: %s', s, paste(capture.output(print(mod$DYNAMICS)),collapse='\n'))
    if (length(mod$PK.PROGRAM)>0) {
        s = sprintf('%s\nKinetic definitions:',s)
        for(nm in names(mod$PK.PROGRAM)) {
            foo = mod$PK.PROGRAM[[nm]]
            if (is.atomic(foo)) txt = foo
            else txt = paste(capture.output(print(foo)),collapse='\n')
            txt = indent.n(txt, 8)
            s = sprintf('%s\n %s =\n%s', s, nm, txt)
        }
    }
    compatible.cov = check.compatibility.of.covariance.structure(mod)
    if (!is.null(mod$PK.COV) && length(mod$PK.COV)>0 && compatible.cov[1]) {
    # Note PK.COV can contain covariance even one unit is in PK block
    # and the other estimable unit is in ERROR block
        s = sprintf('%s\nADDITIIONAL COVARIANCE of internal structure AMONG Estimable Units:',s)
        tb = table.random.parts.list.code(mod$PK.PROGRAM)
        for(pair in mod$PK.COV){
            varn = pair$varnames
            varn.prt = varn
            for(k in 1:2) {
                if (tb$length[ varn[k] ]>1) {
                    varn.prt[k] = sprintf('%s.RANDOM[%s]', varn[k], pair$subind[k])
                }
            }
            s = sprintf('%s\n%s', s,
                indent.n(sprintf('<%s,%s> = %s', varn.prt[1],varn.prt[2], pair$val),8))
        }
    }
    if (length(mod$OBS.PROGRAM)>0) {
        s = sprintf('%s\nObservation definitions:',s)
        for(nm in names(mod$OBS.PROGRAM)) {
            foo = mod$OBS.PROGRAM[[nm]]
            if (is.atomic(foo)) txt = foo
            else txt = paste(capture.output(print(foo)),collapse='\n')
            txt = indent.n(txt, 8)
            s = sprintf('%s\n %s =\n%s', s, nm, txt)
        }
    }
    if (!is.null(mod$OBS.COV) && length(mod$OBS.COV)>0 && compatible.cov[2]) {
        s = sprintf('%s\nADDITIIONAL COVARIANCE of internal structure AMONG random disturbations:',s)
        tb = table.random.parts.list.code(mod$OBS.PROGRAM)
        for(pair in mod$OBS.COV){
            varn = pair$varnames
            varn.prt = varn
            for(k in 1:2) {
                if (tb$length2[ varn[k] ]>1) {
                    varn.prt[k] = sprintf('%s.RANDOM[%s]', varn[k], pair$subind[k])
                }
            }
            s = sprintf('%s\n%s', s,
                indent.n(sprintf('<%s,%s> = %s', varn.prt[1],varn.prt[2], pair$val),8))
        }
    }
    cat(s,'\n')
}

#' Compiler which compiles \code{estimable.model} into a NONMEM like specification
#'
#' This is the main working horse of the exporting phase.
#'
#' @param ms model spec by \code{create.model.spec.from.nmmod.and.nmdata}
#' @return a NONMEM specification, which should be quite familar to NONMEM modelers
compile.model.spec.to.nonmem.spec = function(ms) {
    ##
    theta.bounds = list()
    theta.ests   = numeric(0)
    last.theta.ind = 0
    ##
    last.eta.ind = 0
    omega.blocks = list()
    ##
    last.eps.ind = 0
    sigma.blocks = list()
    ##
    PK.EQNS = list()
    OBS.EQNS = list()
    ##
    if (!is.null(ms$PK.COV) && length(ms$PK.COV)>0) {
    # Note, those PK.COV can even contain covariance between one estimable.unit from PK and the other from ERROR
    # So, should carefully also apply pre.defined.eta.index to those estimable.unit in ERROR code
        tb = table.random.parts.list.code(c(ms$PK.PROGRAM,ms$OBS.PROGRAM))
        ## need to extract from the estimable.unit them own
        nms = names(tb$length)[ tb$length>0 ]
        Omega = list()
        for(nm in nms){
            m = ms$PK.PROGRAM[[ nm ]]
            stopifnot(is(m,'estimable.unit'))
            # estimable unit
            if (is(m,'random.disturbation')){
                Omega[[ length(Omega) + 1 ]] = Cov.random.disturbation(m)
            } else if( is(m,'morphism') || is(m,'name') || is(m,'call')) {
                Omega[[ length(Omega) + 1 ]] = attr(m,'distribution')$Var
            }
        }
        Omega = matrix.direct.sum.list(Omega)
        ## fill in the additional covariance
        for(pair in ms$PK.COV){
            v1 = pair$varnames[1]
            v2 = pair$varnames[2]
            v1.subind = pair$subind[1]
            v2.subind = pair$subind[2]
            v1.ind = tb$prestart[v1] + v1.subind
            v2.ind = tb$prestart[v2] + v2.subind
            Omega[v1.ind,v2.ind] <- Omega[v2.ind,v1.ind] <- pair$val
        }
        # carefully choose the sequence of ETA's so they will become a block matrix
        perm = matrix.blockrize.permutation(Omega)
        pre.defined.eta.index = inverse.permutation(perm)
        has.pk.cov = TRUE
        Omega = Omega[perm,perm]
        if(any(eigen(Omega)$values<=0)) {
            warning('The compiled Covariance matrix OMEGA is not positive definite. The control file may not be able to run under NONMEM Estimation.')
        }
    } else {
        has.pk.cov = FALSE
    }
    ##
    for(nm in names(ms$PK.PROGRAM)) {
        m = ms$PK.PROGRAM[[ nm ]]
        if (is(m,'estimable.unit')) {
            if (is(m,'morphism')) {
                dim.of.m = degree.of.freedom(m)
                paras = NULL
                if (dim.of.m > 0) {
                    theta.ind = last.theta.ind + 1:dim.of.m
                    theta.bounds = c( theta.bounds, as.low.upper.list.morphism.domain(m$domain) )
                    last.theta.ind = last.theta.ind + dim.of.m
                    paras = lapply(theta.ind, function(x) CONS('[',quote(THETA),as.numeric(x)))
                    if (!is.null(attr(m, 'initial.estimation'))) {
                        theta.ests = c( theta.ests, attr(m,'initial.estimation'))   
                    } else { #padding as 0
                        theta.ests = c( theta.ests, rep(0, dim.of.m))
                    }
                }
                e = simplify.2(instantiate.morphism(m, paras))
                if (!is.null(distr<-attr(m,'distribution'))){
                    # should model the PARAMETER nm as random effect
                    rd = convert.distrib.description.to.random.disturbation(distr$type, distr$Var)
                    dim.rd = N.random.effects(rd)
                    if (dim.rd>0) {
                        if (!has.pk.cov) {
                            omega.blocks[[ length(omega.blocks) + 1]] = Cov.random.disturbation(rd)
                        }
                        eta.ind = last.eta.ind + 1:dim.rd
                        last.eta.ind = last.eta.ind + dim.rd
                        if (has.pk.cov) {
                            eta.ind = pre.defined.eta.index[ eta.ind ]
                        }
                        r.paras = lapply(eta.ind, function(x) CONS('[', quote(ETA), as.numeric(x)))
                        e = disturbation.apply(rd, fixed=e, random=r.paras)
                        e = simplify.2(e)
                    }
                }
            } else if (is(m,'random.disturbation')) {
                dim.fixed = N.fixed.effects(m)
                f.para = NULL
                if (dim.fixed > 0) {
                    theta.ind = last.theta.ind + 1:dim.fixed
                    theta.bounds = c( theta.bounds, 
                        as.low.upper.list.morphism.domain(random.disturbation.extract.fixed.domain(m)))
                    last.theta.ind = last.theta.ind + dim.fixed
                    f.para = lapply(theta.ind, function(x) CONS('[',quote(THETA),as.numeric(x)))
                    if (!is.null(attr(m, 'initial.estimation'))) {
                        theta.ests = c( theta.ests, attr(m,'initial.estimation'))   
                    } else { #padding as 0
                        theta.ests = c( theta.ests, rep(0, dim.fixed))
                    }
                }
                dim.rd = N.random.effects(m)
                r.paras = NULL
                if (dim.rd>0) {
                    if (!has.pk.cov){
                        omega.blocks[[ length(omega.blocks) + 1]] = Cov.random.disturbation(m)
                    }
                    eta.ind = last.eta.ind + 1:dim.rd
                    last.eta.ind = last.eta.ind + dim.rd
                    if (has.pk.cov) {
                        eta.ind = pre.defined.eta.index[ eta.ind ]
                    }
                    r.paras = lapply(eta.ind, function(x) CONS('[', quote(ETA), as.numeric(x)))
                }
                e = simplify.2(disturbation.apply(m, fixed=f.para, random=r.paras))
            } else if (is(m,'name') && !is.null(distr<-attr(m,'distribution'))) {
                    e = m
                    rd = convert.distrib.description.to.random.disturbation(distr$type, distr$Var)
                    dim.rd = N.random.effects(rd)
                    if (dim.rd>0) {
                        if (!has.pk.cov){
                            omega.blocks[[ length(omega.blocks) + 1]] = Cov.random.disturbation(rd)
                        }
                        eta.ind = last.eta.ind + 1:dim.rd
                        last.eta.ind = last.eta.ind + dim.rd
                        if (has.pk.cov) {
                            eta.ind = pre.defined.eta.index[ eta.ind ]
                        }
                        r.paras = lapply(eta.ind, function(x) CONS('[', quote(ETA), as.numeric(x)))
                        e = disturbation.apply(rd, fixed=e, random=r.paras)
                        e = simplify.2(e)
                    }
            } else {
                stop('Can not export', paste(deparse(m),collapse=''))
            }
        } else {
            e = m
        }
        # make it pretty
        # exp(THETA[i] + ETA[1]) kept as is, but
        # exp(log(THETA[i] + ETA[1])) => THETA[i] * exp(ETA[1])
        tmp = simplify.collect.log.exp(e)
        if (symbolic.match( quote( exp('?'(anything1)) * '?'(anything2) ), tmp)) tmp = e
        PK.EQNS[[ length(PK.EQNS) + 1 ]] = CONS('=', as.symbol(nm), tmp)
    }
    ##
    if (!is.null(ms$OBS.COV) && length(ms$OBS.COV)>0) {
    # carefully choose the sequence of ETA's so they will become a block matrix
        tb = table.random.parts.list.code(c(ms$PK.PROGRAM,ms$OBS.PROGRAM))
        ## need to extract from the estimable.unit them own
        nms = names(tb$length2)[ tb$length2>0 ]
        Sigma = list()
        for(nm in nms){
            m = ms$OBS.PROGRAM[[ nm ]]
            stopifnot(is(m,'random.observation'))
            if (is(m,'random.disturbation')){
                Sigma[[ length(Sigma) + 1 ]] = Cov.random.disturbation(m)
            }
        }
        Sigma = matrix.direct.sum.list(Sigma)
        ## fill in the additional covariance
        for(pair in ms$OBS.COV){
            v1 = pair$varnames[1]
            v2 = pair$varnames[2]
            v1.subind = pair$subind[1]
            v2.subind = pair$subind[2]
            v1.ind = tb$prestart2[v1] + v1.subind
            v2.ind = tb$prestart2[v2] + v2.subind
            Sigma[v1.ind,v2.ind] <- Sigma[v2.ind,v1.ind] <- pair$val
        }
        perm = matrix.blockrize.permutation(Sigma)
        pre.defined.eps.index = inverse.permutation(perm)
        has.obs.cov = TRUE
        Sigma = Sigma[perm,perm]
        if(any(eigen(Sigma)$values<=0)) {
            warning('The compiled Covariance matrix Sigma is not positive definite. The control file may not be able to run under NONMEM Estimation.')
        }
    } else {
        has.obs.cov = FALSE
    }
    ##
    for(nm in names(ms$OBS.PROGRAM)){
        m = ms$OBS.PROGRAM[[ nm ]]
        if (is(m,'estimable.unit')) {
            warning('Using a estimable unit for a fixed effect in the ERROR block is not usual.
                     If what you really want is to use a existing THETA value, rather than 
                     introduce new THETA or ETA parameters to estimate here.')
            if (is(m,'morphism')) {
                dim.of.m = degree.of.freedom(m)
                paras = NULL
                if (dim.of.m > 0) {
                    theta.ind = last.theta.ind + 1:dim.of.m
                    theta.bounds = c( theta.bounds, as.low.upper.list.morphism.domain(m$domain) )
                    last.theta.ind = last.theta.ind + dim.of.m
                    paras = lapply(theta.ind, function(x) CONS('[',quote(THETA),as.numeric(x)))
                    if (!is.null(attr(m, 'initial.estimation'))) {
                        theta.ests = c( theta.ests, attr(m,'initial.estimation'))   
                    } else { #padding as 0
                        theta.ests = c( theta.ests, rep(0, dim.of.m))
                    }
                }
                e = simplify.2(instantiate.morphism(m, paras))
                if (!is.null(distr<-attr(m,'distribution'))){
                    # should model the PARAMETER nm as random effect
                    rd = convert.distrib.description.to.random.disturbation(distr$type, distr$Var)
                    dim.rd = N.random.effects(rd)
                    if (dim.rd>0) {
                        if (!has.pk.cov) {
                            omega.blocks[[ length(omega.blocks) + 1]] = Cov.random.disturbation(rd)
                        }
                        eta.ind = last.eta.ind + 1:dim.rd
                        last.eta.ind = last.eta.ind + dim.rd
                        if (has.pk.cov) {
                            eta.ind = pre.defined.eta.index[ eta.ind ]
                        }
                        r.paras = lapply(eta.ind, function(x) CONS('[', quote(ETA), as.numeric(x)))
                        e = disturbation.apply(rd, fixed=e, random=r.paras)
                        e = simplify.2(e)
                    }
                }
            } else if (is(m,'random.disturbation')) {
                dim.fixed = N.fixed.effects(m)
                f.para = NULL
                if (dim.fixed > 0) {
                    theta.ind = last.theta.ind + 1:dim.fixed
                    theta.bounds = c( theta.bounds, 
                        as.low.upper.list.morphism.domain(random.disturbation.extract.fixed.domain(m)))
                    last.theta.ind = last.theta.ind + dim.fixed
                    f.para = lapply(theta.ind, function(x) CONS('[',quote(THETA),as.numeric(x)))
                    if (!is.null(attr(m, 'initial.estimation'))) {
                        theta.ests = c( theta.ests, attr(m,'initial.estimation'))   
                    } else { #padding as 0
                        theta.ests = c( theta.ests, rep(0, dim.fixed))
                    }
                }
                dim.rd = N.random.effects(m)
                r.paras = NULL
                if (dim.rd>0) {
                    if (!has.pk.cov) {
                        omega.blocks[[ length(omega.blocks) + 1]] = Cov.random.disturbation(m)
                    }
                    eta.ind = last.eta.ind + 1:dim.rd
                    last.eta.ind = last.eta.ind + dim.rd
                    if (has.pk.cov) {
                        eta.ind = pre.defined.eta.index[ eta.ind ]
                    }
                    r.paras = lapply(eta.ind, function(x) CONS('[', quote(ETA), as.numeric(x)))
                }
                e = simplify.2(disturbation.apply(m, fixed=f.para, random=r.paras))
            } else {
                stop('Can not export', paste(deparse(m),collapse='') )
            }
        } else if (is(m,'multiple.dependent.variable')) {
            r.paras = NULL
            if (is(m,'random.disturbation') && (dim.eps <- N.random.effects(m))>0) {
                eps.ind = last.eps.ind + 1:dim.eps
                last.eps.ind = last.eps.ind + dim.eps
                if (!has.obs.cov) {
                    sigma.blocks[[ length(sigma.blocks) + 1]] = Cov.random.disturbation(m)
                }
                if (has.obs.cov){
                    eps.ind = pre.defined.eps.index[ eps.ind ]
                }
                r.paras = lapply(eps.ind, function(x) CONS('[', quote(EPS), as.numeric(x)))
            }
            e = simplify.2(multiple.dependent.variable.apply(m, fixed=NULL, random=r.paras))
        } else if (is(m,'random.observation')) {
            r.paras = NULL
            if (is(m,'random.disturbation') && (dim.eps <- N.random.effects(m))>0) {
                eps.ind = last.eps.ind + 1:dim.eps
                last.eps.ind = last.eps.ind + dim.eps
                if (!has.obs.cov) {
                    sigma.blocks[[ length(sigma.blocks) + 1]] = Cov.random.disturbation(m)
                }
                if (has.obs.cov){
                    eps.ind = pre.defined.eps.index[ eps.ind ]
                }
                r.paras = lapply(eps.ind, function(x) CONS('[', quote(EPS), as.numeric(x)))
            }
            e = simplify.2(disturbation.apply(m, fixed=NULL, random=r.paras))
        } else {
            e = m
        }
        OBS.EQNS[[ length(OBS.EQNS) + 1 ]] = CONS('=', as.symbol(nm), e)
    }
    theta = matrix(0, last.theta.ind, 3)
    theta[,1] = sapply(theta.bounds, function(x) x[1])
    theta[,3] = sapply(theta.bounds, function(x) x[2])
    theta[,2] = theta.ests
    ### New added Code for output pretty names
    tab = list.all.fixed.parameters.estimable.model(ms)
    if (NROW(tab)>0) rownames(theta)=sprintf('%s[%s]', tab$varname, tab$subind)
    tab = list.all.random.parameters.estimable.model(ms)
    omega.names.tab = tab[tab$indexed.by=='subject',]
    omega.names = sprintf('%s[%s]', omega.names.tab$varname, omega.names.tab$subind)
    if (has.pk.cov) omega.names = omega.names[ inverse.permutation(pre.defined.eta.index) ]
    sigma.names.tab = tab[tab$indexed.by=='observation',]
    sigma.names = sprintf('%s[%s]', sigma.names.tab$varname, sigma.names.tab$subind)
    if (has.obs.cov) sigma.names = sigma.names[ inverse.permutation(pre.defined.eps.index) ]
    ###
    re =list(SUBR = as.nonmem.subroutine.dynamics(ms$DYNAMICS),
             PK   =PK.EQNS, 
             ERROR=OBS.EQNS,
             THETA=theta, 
             OMEGA=if (has.pk.cov) matrix.decompose.to.list(Omega) else omega.blocks, 
             SIGMA=if (has.obs.cov) matrix.decompose.to.list(Sigma) else sigma.blocks
    )
    # usually in looking back from a NMRUN, usually this part do not need to be print
    if (has.pk.cov) {
        re$pre.defined.eta.index = pre.defined.eta.index
    }
    if (has.obs.cov) {
        re$pre.defined.eps.index = pre.defined.eps.index
    }
    ### New added Code for output pretty names
    if (length(re$OMEGA)>0) {
        re$OMEGA = inject.names.into.list.of.matrix(re$OMEGA, omega.names)
    }
    if (length(re$SIGMA)>0) {
        re$SIGMA = inject.names.into.list.of.matrix(re$SIGMA, sigma.names)
    }
    ###
    class(re) = 'nonmem.model.specification'
    re
}

#' inverse permutation
#' @param perm
#' @return inverse of perm
inverse.permutation = function(perm){
    re = seq_along(perm)
    re[perm] = re
    re
}

#' S3 printer
print.nonmem.model.specification = function(ms){
    indent.n = function(ss,n) {
       if (n==0) return(ss)
       indent = paste(rep(' ',n),collapse='')
       ss = gsub('^', indent, ss)
       gsub('\n', sprintf('\n%s',indent), ss)
    }
    if (length(ms$SUBR)>0) {
        cat(sprintf('SUBR %s', ms$SUBR), '\n')
    }
    if (length(ms$PK)>0){
        # PK BLOCK
        tmp = CONS('=',quote(PK),
                    as.call(c(as.symbol('{'),ms$PK)))
        txt = paste(capture.output(print(tmp)),collapse='\n')
        cat(txt,'\n')
    }
    if (length(ms$ERROR)>0) {
        # ERROR BLOCK
        tmp = CONS('=',quote(ERROR),
                    as.call(c(as.symbol('{'),ms$ERROR)))
        txt = paste(capture.output(print(tmp)),collapse='\n')
        cat(txt,'\n')
    }
    if (length(ms$THETA)>0){
        # THETA formater
        cat('THETA','\n')
        tab = as.table(ms$THETA)
        colnames(tab)=c('Lower','Est','Upper')
        if (is.null(rownames(tab))) {
            rownames(tab)=sprintf('THETA[%s]',1:NROW(tab))
        }
        txt = paste(capture.output(print(tab)),collapse='\n')
        cat(indent.n(txt,4), '\n')
    }
    if (length(ms$OMEGA)>0) {
        # OMEGA formater
        cat('OMEGA','\n')
        txt = print.diagonal.block.matrix0(ms$OMEGA, prefix='ETA')
        cat(indent.n(txt,4), '\n')
    }
    if (length(ms$SIGMA)>0) {
        # SIGMA
        cat('SIGMA','\n')
        txt = print.diagonal.block.matrix0(ms$SIGMA, prefix='EPS')
        cat(indent.n(txt,4), '\n')
    }
    invisible()
}

print.diagonal.block.matrix0 = function(ml, prefix='ETA'){
    extra.nms = overall.dim.names.matrix.list(ml, prefix=prefix)
    m = matrix.direct.sum.list(ml)   
    m[m==0]=''
    m1 = cbind('[',m,']')
    tab = as.table(m1)
    if (is.null(extra.nms)) {
        nm = sprintf('%s[%s]',prefix, 1:NCOL(m))
    } else {
        nm = extra.nms
    }
    colnames(tab) = c('',nm,'')
    rownames(tab) = nm
    txt = capture.output(print(tab))
    paste(txt,collapse='\n')
}
