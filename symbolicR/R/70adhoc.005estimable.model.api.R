#' Get the index(position) of a parameter of a estimable unit when it's part of a model
#'
#' One estimable unit can be meaningful on it own.
#' But we usually see it as an embeded object in an estimable model.
#' When we list all the parameters (corresponding to those THETAs, ETAs) in a big table, it's useful to know their positions.
#'
#' @param esm estimable model
#' @param v1 estimable.unit 1
#' @param v2 estimable.unit 2
#' @param v1.subind one estimable.unit may have more than one random varible in it, this indicate which variable are operated on
#' @param v2.subind similar to \code{v1.subind}
#' @return embeded index and the corresponding backref (if available, or just NAs)
embeded.ind.of.two.estimable.unit = function(esm,v1,v2,v1.subind=1,v2.subind=1) {
    tb = table.random.parts.list.code(c(esm$PK.PROGRAM,esm$OBS.PROGRAM))
    ind = c(tb$prestart[[ v1 ]] + v1.subind, tb$prestart[[ v2 ]] + v2.subind)
    bk = tb$backref[ ind ]
    list(ind=ind,backref=bk)
}

#' Check if the additional convariance is compatible with the model definition
#'
#' We allow user to change the estimable model,
#' usually, the distribution type or covariance matrix of parameters of \code{estimable.model} can be changed.
#' Even some \code{estimable.unit}s or \code{random.observation}s can be added into or removed from the \code{estimable.model}.
#' The change may cause conflications between the \code{PK.PROGRAM} or \code{OBS.PROGRAM} and
#'  the additional convariance structures defined in \code{PK.COV} or \code{OBS.COV}.
#' This function will check the situation and print a warning.
#' @param esm estimable.model
#' @param need.print if print diagonize information
#' @return vector of TRUE or FALSE indication if compatible for PK.COV or OBS.COV
check.compatibility.of.covariance.structure = function(esm, need.print=TRUE) {
    tb = table.random.parts.list.code(c(esm$PK.PROGRAM,esm$OBS.PROGRAM))
    pk.incompatible = FALSE
    obs.incompatible = FALSE
    if (length(esm$PK.COV)>0){
        rd = tb$length[ tb$length > 0 ]
        nms = names(rd)
        for(pair in esm$PK.COV) {
            for(k in seq_along(pair$varnames)) {
                if (!(pair$varnames[k] %in% nms)) {
                    pk.incompatible = TRUE
                    if (need.print) {
                        cat(sprintf('[Error] %s is not a estimable.unit in the model\n', pair$varnames[k]))
                    }
                } else if (pair$subind[k] > rd[ pair$varnames[k] ]) {
                    pk.incompatible = TRUE
                    if (need.print) {
                        cat(sprintf('[Error] subindex %s is larger than number of random effect of
                         estimable.unit %s\n', pair$subind[k], pair$varnames[k]))
                    }
                }
            }
        }
    }
    if (length(esm$OBS.COV)>0){
        rd = tb$length2[ tb$length2 > 0 ]
        nms = names(rd)
        for(pair in esm$OBS.COV) {
            for(k in seq_along(pair$varnames)) {
                if (!(pair$varnames[k] %in% nms)) {
                    obs.incompatible = TRUE
                    if (need.print) {
                        cat(sprintf('[Error] %s is not a random.observation in the model\n', pair$varnames[k]))
                    }
                } else if (pair$subind[k] > rd[ pair$varnames[k] ]) {
                    obs.incompatible = TRUE
                    if (need.print) {
                        cat(sprintf('[Error] subindex %s is larger than number of random effect of
                         random.observation %s\n', pair$subind[k], pair$varnames[k]))
                    }
                }
            }
        }
    }
    !c(pk.incompatible,obs.incompatible)
}

#' Check if the meta info is compatible with the \code{PK.PROGRAM} or \code{OBS.PROGRAM} code in the model
#'
#' The story is, after deleting some \code{estimable.unit} which is essential for the model definition, 
#' it will cause problem. We will check before convert it into a control file.
#' @param esm estimable.model
#' @param need.print print diagonize information
#' @return TRUE or FALSE
check.meta.information = function(esm, need.print=TRUE) {
    tmp = lapply(esm$META.INFO$PK.PARAMETERS, function(x) paste(deparse(x),collapse=''))
    ind = which(!(tmp %in% names(esm$PK.PROGRAM)))
    pk.ok = TRUE
    if (length(ind)>0) {
        pk.ok = FALSE
        if (need.print) {
            cat(sprintf('[Error] The PK parameters [%s] are declared in META.INFO, but not appearing in PK.PROGRAM\n',
                paste(tmp[ind],collapse=',')))

        }
    }
    others.ok = TRUE
    if (length(esm$META.INFO$OTHER.PARAMETERS)>0){
        tmp = lapply(esm$META.INFO$OTHER.PARAMETERS, function(x) paste(deparse(x),collapse=''))
        ind = which(!(tmp %in% c(names(esm$PK.PROGRAM),names(esm$OBS.PROGRAM))))
        if (length(ind)>0) {
            others.ok = FALSE
            if (need.print) {
                cat(sprintf('[Error] The Some Other parameters [%s] are declared in META.INFO, 
                but not appearing in either PK.PROGRAM or OBS.PROGRAM\n',
                    paste(tmp[ind],collapse=',')))

            }
        }
    }
    c(pk.ok, others.ok)
}

#' check if all basic pk parameter needed by subroutine and trans appearing in the PK.PROGRAM
#' @param esm estimable.model
#' @param need.print print diagonize information
#' @return TRUE or FALSE
check.isvalid.basic.pk.parameters = function(esm, need.print=TRUE){
    pks = names(esm$PK.PROGRAM)
    nms = get.basic.pk.parameters(esm$DYNAMICS)
    subs = attr(nms,'subs')
    trans = attr(nms,'trans')
    nms.ind = which(! nms %in% pks)
    ok = TRUE
    if (length(nms.ind)>0) {
        ok = FALSE
        if (need.print) {
            cat(sprintf('[Error] Subroutine [%s,%s] in the Model need [%s] basic pk parameters
                         while only [%s] appearing in the PK.PROGRAM, lacking [%s].\n',
                            subs,trans,
                            paste(nms,collapse=','),
                            paste(pks[pks %in% nms],collapse=','),
                            paste(nms[nms.ind], collapse=',')
         ))
        }
    }
    ok
}

#' For dynamics description, return the necessary basic pk parameters
#'
#' The basic pk parameters are necessary for NONMEM PREDPP subroutines to run.
#' Any lack of those parameters will make NONMEM complain.
#' @param dyn dynamics
#' @return basic pk parameters
get.basic.pk.parameters = function(dyn) {
    subs0 = get.subs.trans.dynamics(dyn)
    subs  = subs0[1]
    trans = subs0[2]
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
    nms = sapply(subr$parameters$basic[[ trans ]]$names, function(x) paste(deparse(x),collapse=''))
    attr(nms,'subs') = subs
    attr(nms,'trans') = trans
    nms
}

#' add covariance relations for inner random variable of \code{estimable.unit}s in an \code{estimable.model}
#'
#' @note This function is a scaler version, all arguments should be scalar.
#' @param esm estimable model
#' @param v1 estimable.unit 1
#' @param v2 estimable.unit 2
#' @param x a new value, should be one number
#' @param v1.subind one estimable.unit may have more than one random varible in it, this indicate which variable are operated on
#' @param v2.subind similar to \code{v1.subind}
#' @param do.check check if the input \code{v1,v2,v1.subind,v2.subind} are valid
#' @return the updated estimable.model
add.covariance.of.estimable.unit.estimable.model = function(esm,v1,v2,x,v1.subind=1,v2.subind=1,do.check=TRUE) {
    if (do.check) {
        tb = table.random.parts.list.code(c(esm$PK.PROGRAM,esm$OBS.PROGRAM))
        if (!(v1 %in% names(tb$length))) {
            stop(sprintf('v1: [%s] is not random part of the model PK',v1))
        }
        if (v1.subind > tb$length[ v1 ]) {
            stop(sprintf('v1 index %s out of bound',v1.subind))
        }
        if (!(v2 %in% names(tb$length))){
            stop(sprintf('v2: [%s] is not random part of the model PK',v2))
        }  
        if (v2.subind > tb$length[ v2 ]) {
            stop(sprintf('v2 index %s out of bound',v2.subind))
        }
    }
    if (v1==v2) {
        stop(sprintf('Please update the internal covariance of %s directly using update.cov.estimable.unit',v1))
    }
    cov.l = esm$PK.COV
    if (is.null(cov.l)) cov.l = list()
    ind = NULL
    subind = c(v1.subind,v2.subind)
    for(i in seq_along(cov.l)) {
        if (setequal(c(v1,v2), cov.l[[i]]$varnames) && 
            all( cov.l[[i]]$subind== subind)) {
            ind = i
            break()
        }
    }
    if (is.null(ind)) {
        cov.l[[ length(cov.l)+1 ]] = 
            list(varnames=c(v1,v2),
                 subind=subind,
                 val = x)
    }else{
        cov.l[[ ind ]] = 
            list(varnames=c(v1,v2),
                 subind=subind,
                 val = x)
    }
    esm$PK.COV=cov.l
    esm
}

#' remove covariance relations for inner random variable of \code{estimable.unit}s in an \code{estimable.model}
#'
#' @note scaler version, each argument has to be a scalar.
#' @param esm estimable model
#' @param v1 estimable.unit 1
#' @param v2 estimable.unit 2
#' @param v1.subind one estimable.unit may have more than one random varible in it, this indicate which variable are operated on
#' @param v2.subind similar to \code{v1.subind}
#' @return the updated estimable.model
remove.covariance.of.estimable.unit.estimable.model = function(esm,v1,v2,v1.subind=1,v2.subind=1) {
    if (v1==v2) {
        stop(sprintf('Please remove the internal covariance of %s directly using update.cov.estimable.unit',v1))
    }
    cov.l = esm$PK.COV
    if (is.null(cov.l)) return(esm)
    ind = NULL
    subind = c(v1.subind,v2.subind)
    for(i in seq_along(cov.l)) {
        if (setequal(c(v1,v2), cov.l[[i]]$varnames) && 
            all( cov.l[[i]]$subind== subind)) {
            ind = i
            break()
        }
    }
    if (!is.null(ind)) {
        esm$PK.COV = cov.l[ - ind ]
    }
    esm
}

#' remove covariance relations for inner random variable of random.observation in an estimable.model
#'
#' @note scaler version, each argument has to be a scalar
#' @param esm estimable model
#' @param v1 estimable.unit 1
#' @param v2 estimable.unit 2
#' @param v1.subind one estimable.unit may have more than one random varible in it, this indicate which variable are operated on
#' @param v2.subind similar to \code{v1.subind}
#' @return the updated estimable.model
remove.covariance.of.random.observation.estimable.model = function(esm,v1,v2,v1.subind=1,v2.subind=1) {
    if (v1==v2) {
        stop(sprintf('Please remove the internal covariance of %s directly using update.cov.random.observation',v1))
    }
    cov.l = esm$OBS.COV
    if (is.null(cov.l)) return(esm)
    ind = NULL
    subind = c(v1.subind,v2.subind)
    for(i in seq_along(cov.l)) {
        if (setequal(c(v1,v2), cov.l[[i]]$varnames) && 
            all( cov.l[[i]]$subind== subind)) {
            ind = i
            break()
        }
    }
    if (!is.null(ind)) {
        esm$OBS.COV = cov.l[ - ind ]
    }
    esm
}

#' Batch removing the random effects in an \code{estimable.model}
#'
#' Remember the function \code{\link{remove.random.effect.estimable.unit}} can specialize an \code{estimable.unit} into simpler one without random effect.
#' But it's boring to do it in a \code{for} loop by usr, to update any possible \code{estimable.unit} in a model.
#' So we provide this API for the purpose. The usage is simple, put the names of the \code{estimable.unit} in argument \code{varnames},
#  and all the specified ones is changed. You can also omit \code{varnames} and set \code{all} to \code{TRUE} to remove all random effects.
#' @param esm estimable.model
#' @param varnames variable names
#' @return new model
remove.random.effects.estimable.model = function(esm, varnames, all=FALSE) {
    if (missing(varnames) && !all) {
        stop("no input for which random effect to be removed, and not remove all")
    }
    pks = names(esm$PK.PROGRAM)
    obs = names(esm$OBS.PROGRAM)
    if (!all) {
        pks = pks[ pks %in% varnames ]
        obs = obs[ obs %in% varnames ]
    }
    for (vn in pks) {
        foo = esm$PK.PROGRAM[[vn]]
        if (is(foo, 'estimable.unit')) {
            esm$PK.PROGRAM[[ vn ]] = remove.random.effect.estimable.unit(foo)
        } 
    }
    if (length(esm$PK.COV)>0) {
        rm.ind = which(sapply(esm$PK.COV, function(x) any(x$varnames %in% pks)))
        if (length(rm.ind)>0) {
            esm$PK.COV = esm$PK.COV[ - rm.ind ]
        }
    }
    for (vn in obs) {
        foo = esm$OBS.PROGRAM[[vn]]
        if (is(foo, 'estimable.unit')) {
            esm$OBS.PROGRAM[[ vn ]] = remove.random.effect.estimable.unit(foo)
        } else if (is(foo, 'random.disturbation')) {
            esm$OBS.PROGRAM[[ vn ]] = remove.random.effect.random.disturbation(foo)
        }
    }
    if (length(esm$OBS.COV)>0) {
        rm.ind = which(sapply(esm$OBS.COV, function(x) any(x$varnames %in% obs)))
        if (length(rm.ind)>0) {
            esm$OBS.COV = esm$OBS.COV[ - rm.ind ]
        }
    }
    esm
}

#' add covariance relations for inner random observations 
#'
#' @note scaler version, each argument has to be a scalar.
#' @param esm estimable model
#' @param v1 estimable.unit 1
#' @param v2 estimable.unit 2
#' @param x a new value, should be one number
#' @param v1.subind one estimable.unit may have more than one random varible in it, this indicate which variable are operated on
#' @param v2.subind similar to \code{v1.subind}
#' @param do.check check if the input \code{v1,v2,v1.subind,v2.subind} are valid
#' @return the updated estimable.model
add.covariance.of.random.observation.estimable.model = function(esm,v1,v2,x,v1.subind=1,v2.subind=1,do.check=TRUE) {
    if (do.check) {
        tb = table.random.parts.list.code(esm$OBS.PROGRAM)
        if (!(v1 %in% names(tb$length2))) {
            stop(sprintf('v1: [%s] is not random part of the model PK',v1))
        }
        if (v1.subind > tb$length2[ v1 ]) {
            stop(sprintf('v1 index %s out of bound',v1.subind))
        }
        if (!(v2 %in% names(tb$length2))){
            stop(sprintf('v2: [%s] is not random part of the model PK',v2))
        }  
        if (v2.subind > tb$length2[ v2 ]) {
            stop(sprintf('v2 index %s out of bound',v2.subind))
        }
    }
    if (v1==v2) {
        stop(sprintf('Please update the internal covariance of %s directly using update.cov.random.observation',v1))
    }
    cov.l = esm$OBS.COV
    if (is.null(cov.l)) cov.l = list()
    ind = NULL
    subind = c(v1.subind,v2.subind)
    for(i in seq_along(cov.l)) {
        if (setequal(c(v1,v2), cov.l[[i]]$varnames) && 
            all( cov.l[[i]]$subind== subind)) {
            ind = i
            break()
        }
    }
    if (is.null(ind)) {
        cov.l[[ length(cov.l)+1 ]] = 
            list(varnames=c(v1,v2),
                 subind=subind,
                 val = x)
    }else{
        cov.l[[ ind ]] = 
            list(varnames=c(v1,v2),
                 subind=subind,
                 val = x)
    }
    esm$OBS.COV=cov.l
    esm
}

#' update the initial.estimations for all estimable unit
#'
#' For NONMEM users, it's common to consider all the parameters need to be estimated in a whole, no matter if 
#' \code{THETA[1]} and \code{THETA[2]} are both for modeling the first PK parameter \code{CL}.
#' So, we provide this API, user can use a vector of values \code{thetas} and used sequentially as the initial estimation for
#' all \code{estimable.unit}s
#' @param esm estimable.model
#' @param thetas a vector
#' @return updated estimable model
update.all.initial.estimation.estimable.unit.in.estimable.model = function(esm, thetas){
    if (length(thetas)<1) return(esm)
    # recycle thetas, usefull if you just want to set all to zero
    ind.theta = 0
    N.theta = length(thetas)
    if (length(esm$PK.PROGRAM)>0) {
        for(i in seq_along(esm$PK.PROGRAM)) {
            m = esm$PK.PROGRAM[[i]]
            if (!is(m,'estimable.unit')) next()
            # estimable.unit might be call/name/atomic also, only for the RANDOM purpose
            N = 0
            if (is(m,'morphism')) {
                N = degree.of.freedom(m)
            } else if (is(m,'random.disturbation')) {
                N = N.fixed.effects(m)
            }
            if (N==0) next()
            vec = thetas[ 1 + (ind.theta + 0:(N-1)) %% N.theta ]
            ind.theta = ind.theta + N
            attr(m,'initial.estimation') = vec 
            esm$PK.PROGRAM[[i]] = m
        }
    }
    if (length(esm$OBS.PROGRAM)>0) {
        for(i in seq_along(esm$OBS.PROGRAM)) {
            m = esm$OBS.PROGRAM[[i]]
            if (!is(m,'estimable.unit')) next()
            # estimable.unit might be call/name/atomic also, only for the RANDOM purpose
            N = 0
            if (is(m,'morphism')) {
                N = degree.of.freedom(m)
            } else if (is(m,'random.disturbation')) {
                N = N.fixed.effects(m)
            }
            if (N==0) next()
            vec = thetas[ 1 + (ind.theta + 0:(N-1)) %% N.theta ]
            ind.theta = ind.theta + N
            attr(m,'initial.estimation') = vec 
            esm$OBS.PROGRAM[[i]] = m
        }
    }
    esm
}

#' fill in the estimable unit or random observations 
#'
#' All inner random parameters are views as big random vector. And the \code{cov.} is viewed as
#' the covariance of the vector to update the model.
#' @param esm estimable.model
#' @param cov. covariance matrix
#' @param indexed.by can be 'subject' or 'observation'
#' @return new estimable.model
update.cov.estimable.model = function(esm, cov. , indexed.by = 'subject') {
    tab = list.all.random.parameters.estimable.model(esm)
    tab = tab[ tab$indexed.by == indexed.by, ]
    if (!is.matrix(cov.)) cov. = as.matrix(cov.)
    if (any(eigen(cov.)$values<=0)) {
        warning('cov. is not positive definite.')
    }
    if (NROW(cov.) != NROW(tab)) {
        stop(sprintf('The random effects indexed by %s size mismatch with input covariance matrix.', indexed.by))
    }
    # internal part of one estimable.unit or random.observation
    rvs = unique(tab$varname)
    for(vn in rvs) {
        # assume cov. should have index
        ind = which(vn == tab$varname)
        place = tab$defined.place[ind[1]]
        ## two class
        foo = esm[[place]][[vn]]
        if (is(foo,'estimable.unit')) {
            esm[[ place ]][[ vn ]] = update.cov.estimable.unit( foo , cov.[ind,ind] )
        } else if (is(foo,'random.observation')) {
            esm[[ place ]][[ vn ]] = update.cov.random.observation( foo , cov.[ind,ind] )
        } else {
            stop(sprintf('can not update cov. of class [%s].', paste(class(foo),collapse=',')))
        }
    }
    # extra covariance among different estimable.unit
    n.rvs = length(rvs)
    # Note if you put estimable.unit in OBS.PROGRAM, the cov still stored in PK.COV
    # so, place only affect by to which index the random effect is indexed
    if (indexed.by == 'subject') {
        place = 'PK.COV' 
    } else {
        place = 'OBS.COV'
    }
    cov.l = esm[[ place ]]
    if (is.null(cov.l)) cov.l = list()
    if (n.rvs > 1) {
        for(v1.ind in 1:(n.rvs-1)){
            for(v2.ind in (v1.ind+1):n.rvs){
                v1 = rvs[ v1.ind ]
                v2 = rvs[ v2.ind ]
                v1.inds = which(tab$varname==v1)
                v2.inds = which(tab$varname==v2)
                v1.subinds = tab[v1.inds,'subind']
                v2.subinds = tab[v2.inds,'subind']
                for(v1.subind in v1.subinds) {
                    for(v2.subind in v2.subinds) {
                        global.ind = c(v1.inds[ v1.subind ], v2.inds[ v2.subind ] )
                        new.val = cov.[global.ind[1],global.ind[2]]
                        subind = c(v1.subind, v2.subind)
                        cov.l.ind = NULL
                        for(i in seq_along(cov.l)) {
                            if (setequal(c(v1,v2), cov.l[[i]]$varnames) && 
                                all( cov.l[[i]]$subind==subind)) {
                                cov.l.ind = i
                                break()
                            }
                        }
                        if (is.null(cov.l.ind)) {  
                            if (isTRUE(abs(new.val) > .Machine$double.eps)) {
                                cov.l[[ length(cov.l)+1 ]] = 
                                    list(varnames=c(v1,v2),
                                         subind=subind,
                                         val = new.val)
                            }
                        } else {
                            cov.l[[ cov.l.ind ]] = 
                                list(varnames=c(v1,v2),
                                     subind=subind,
                                     val = new.val)
                        }
                    }
                }
            }
        }
    }
    # remove those 0 ones
    cov.l = cov.l[ unlist(sapply(cov.l, function(x) isTRUE(abs(x$val)>.Machine$double.eps))) ]
    esm[[ place ]]  = cov.l
    if (length(cov.l)==0) {
    # remove it if does not have entries
        esm[[ place ]] = NULL
    }
    esm
}

#' list all parameters of fixed part in all \code{estimable.unit} from a list
#' @param le list of code
#' @return data.frame
list.all.fixed.parameters.list = function(le) {
    nms = names(le)
    re = character(0)
    subind = integer(0)
    for(i in seq_along(le)) {
        eu = le[[i]]
        if (!is(eu,'estimable.unit')) next()
        if (is.atomic(eu) || is(eu,'name') || is(eu,'call')) next()
        if (is(eu,'morphism')) {
            if (eu$defn$type=='Constant') next()
            N = degree.of.freedom(eu)
            if (N==0) next()
            re = c(re, rep(nms[i], N))
            subind = c(subind, 1:N)
            next()
        }
        if (is(eu,'random.disturbation')) {
            N = N.fixed.effects(eu)
            if (N==0) next()
            re = c(re, rep(nms[i], N))
            subind = c(subind, 1:N)
            next()
        }
        stop(sprintf('do not know how to get parameters of estimable.unit of class [%s]', 
            paste(class(eu),collapse=',')))
    }
    if (length(re)==0) {
        return(
            data.frame(INDEX=integer(0), varname=character(0), subind=integer(0), stringsAsFactors=FALSE))
    }
    data.frame(INDEX=1:length(re), varname=re, subind=subind, stringsAsFactors=FALSE)
}

#' list all fixed parameters in a model spec
#' @param esm estimable.model
#' @param with.triple output lower est upper of fixed parameters
#' @return data.frame
list.all.fixed.parameters.estimable.model = function(esm, with.triple=FALSE) {
    L = c(esm$PK.PROGRAM, esm$OBS.PROGRAM)
    tab = list.all.fixed.parameters.list(L)
    nms = names(esm$PK.PROGRAM)
    tab$defined.place = ifelse(tab$varname %in% nms, 'PK.PROGRAM', 'OBS.PROGRAM')
    if (with.triple) {
        varnames = unique(tab$varname)
        foo = NULL
        for(vn in varnames) {
            where.defn = tab[ which(tab$varname==vn)[1] , 'defined.place' ]
            m = esm[[ where.defn ]][[ vn ]]
            foo = rbind(foo, get.triple.estimable.unit(m))
        }
        tab$Lower = foo[,1]
        tab$Est   = foo[,2]
        tab$Upper = foo[,3]
    }
    tab
}

#' list all random parameters in a model spec
#' @param esm estimable.model
#' @return data.frame
list.all.random.parameters.estimable.model = function(esm){
    tb = table.random.parts.list.code(c(esm$PK.PROGRAM,esm$OBS.PROGRAM))
    etas = tb$length
    etas = etas[ etas > 0 ]
    re = NULL
    if (length(etas)>0) {
        re = data.frame(INDEX=1:sum(etas), 
                        varname=tb$lookup, 
                        subind=unlist(sapply(unname(etas), function(x) 1:x)),
                        indexed.by = 'subject',
                        stringsAsFactors=FALSE
                        )
    }
    epss = tb$length2
    epss = epss[ epss > 0 ]
    if (length(epss)>0) {
        re1 = data.frame(INDEX=1:sum(epss), 
                    varname=tb$lookup2, 
                    subind=unlist(sapply(unname(epss), function(x) 1:x)),
                    indexed.by = 'observation',
                    stringsAsFactors=FALSE
                    )
        if (is.null(re)) {
            re = re1
        } else {
            re = rbind.data.frame(re, re1)
        }
    }
    if (is.null(re)) {
    # ensure return a data.fram evenif non-random effect
        re = data.frame(INDEX=integer(0),
                        varname=character(0),
                        subind=integer(0),
                        indexed.by=character(0),
                        defined.place=character(0),
                        stringsAsFactors=FALSE)
        return(re)
    }
    nms = names(esm$PK.PROGRAM)
    re$defined.place = ifelse(re$varname %in% nms, 'PK.PROGRAM', 'OBS.PROGRAM')
    re
}

#' return covariance matrix of all random parameters for the model spec
#' @param esm estimable.model
#' @param indexed.by subject or observation random effect
#' @return named covariance matrix
get.cov.estimable.model = function(esm, indexed.by='subject') {
    tab = list.all.random.parameters.estimable.model(esm)
    tab = tab[tab$indexed.by==indexed.by, ]
    cov. = list()
    for(vn in unique(tab$varname)) {
        where.defn = tab[ which(tab$varname==vn)[1], 'defined.place' ]
        m = esm[[ where.defn ]][[ vn ]]
        if (is(m,'estimable.unit')) {
            if (is(m,'random.disturbation')) {
                cov.[[length(cov.)+1]] = Cov.random.disturbation(m)
            } else {
                cov.[[length(cov.)+1]] = attr(m,'distribution')$Var
            }
        } else if (is(m,'random.observation')) {
            cov.[[length(cov.)+1]] = Cov.random.disturbation(m)
        } else {
            stop(sprintf('The [%s] of class [%s] should not have random part.', vn, paste(class(m),collapse=',')))
        }
    }
    cov. = matrix.direct.sum.list( cov. )
    if (indexed.by == 'subject' ) {
        place = 'PK.COV'
    } else {
        place = 'OBS.COV'
    }
    if (!is.null(esm[[ place ]])) {
        for(pair in esm[[ place ]] ) {
            v1 = pair$varnames[1]
            v2 = pair$varnames[2]
            v1.subind = pair$subind[1]
            v2.subind = pair$subind[2]
            v1.ind = which(tab$varname==v1 & tab$subind==v1.subind)
            v2.ind = which(tab$varname==v2 & tab$subind==v2.subind)
            cov.[v1.ind, v2.ind] <- cov.[v2.ind,v1.ind] <- pair$val
        }
    }
    nms = sprintf('%s[%s]',tab$varname,tab$subind)
    colnames(cov.) = nms
    rownames(cov.) = nms
    cov.
}

#' list all paramters in one model spec
#' @param esm estimable.model
#' @return data.frame
list.all.parameters.estimable.model = function(esm){
    fixed = list.all.fixed.parameters.estimable.model(esm, with.triple=FALSE)
    random = list.all.random.parameters.estimable.model(esm)
    if (NROW(fixed)>0) {
        fixed$indexed.by = NA
        fixed$type       = 'fixed'
    }
    if (NROW(random)>0){
        random$type      = 'random'
    }
    rbind.data.frame(fixed,random)
}

#' match a simple assignment pattern
#' @param e expression
#' @return lhs and rhs
pattern.expression.simple.assignment = function(e) {
    if (!symbolic.match( CONS('=',quote('?'(any.lhs)),quote('?'(any.rhs))), e)) {
        stop('The Code for Inserting should be an assignment.')
    }
    if (!(is.symbol(any.lhs)) || 
        (is.call(any.lhs) && any.lhs[[1]]=='[' && length(any.lhs)>2 && is.symbol(any.lhs[[2]]))) {
        stop('The lhs of assignment should be a simplie variable or an array.')
    }
    c(any.lhs,any.rhs)   
}

#' insert code into a given list
#'
#' Insert one line code into list. The code is an assignment, which the \code{lhs} and \code{rhs} should be given explicitly.
#' The \code{rhs} can be an atomic or a symbol or an expression or supported object like \code{estimable.unit} and \code{random.observation}.
#'
#' If both \code{lhs} and \code{rhs} are missing, you have to provide an assignment expression in \code{e} and the function
#' can split it into lhs and rhs.
#' @param le list of code
#' @param lhs lhs symbol
#' @param rhs rhs definition, can be estimable unit
#' @param e assignment expression
#' @param index where to insert, should be number
#' @return new list
insert.one.line.code.list = function(le, lhs, rhs, e, index ) {
    N = length(le)
    if (missing(index)) index = N + 1
    if (index<1 || index>N+1) {
        stop('Index to insert is out of range.')
    }
    if (missing(lhs) || missing(rhs)) {
        if (missing(e)) {
            stop('If not lhs and rhs give, one assignment statement is expected.')
        }
        e1 = pattern.expression.simple.assignment(e)
        lhs = e1[[1]]
        rhs = e1[[2]]
    }
    nms = names(le)
    if (is.null(nms)) nms = as.character(1:N)
    symbolic.names = attr(le,'symbolic.names')
    if (is.null(symbolic.names)) symbolic.names = lapply(nms,as.symbol)
    lhs.name = asCharacterSymbolCall(lhs)
    new.code = vector(N+1, mode='list')
    if (index==1) {
        new.code[[1]]=rhs
        new.code[2:(N+1)] = le
        symbolic.names  = c(lhs,       symbolic.names   )
        nms             = c(lhs.name,  nms            )
    } else if (index == N + 1 ) {
        new.code[1:N]   = le
        new.code[[N+1]] = rhs
        symbolic.names  = c(symbolic.names,          lhs)
        nms             = c(nms,              lhs.name)  
    } else {
        new.code[1:(index-1)] = le[1:(index-1)]
        new.code[[index]]     = rhs
        new.code[(index+1):(N+1)] = le[index:N]
        symbolic.names = c(symbolic.names[1:(index-1)],lhs,       symbolic.names[index:N])
        nms            = c(nms[1:(index-1)],         lhs.name,  nms[index:N]         )
    }
    attr(new.code,'symbolic.names') = symbolic.CMatrix(symbolic.names,nrow=1)
    names(new.code) = nms
    new.code
}

#' insert code into model spec
#' 
#' Handy character representation of place.
#' \code{^} refer to the beginning of code, \code{$} means appending after the code
#' \code{^CL} or \code{CL} means insert the new code before \code{CL} in the existing code
#' \code{CL$} means after that code
#'
#' @param esm estimable.model
#' @param lhs the lhs of inserted code, should be quoted expression
#' @param rhs definition of lhs, can be a expression or estimable.unit
#' @param which.part insert to which part, PK or Observation PROGRAM
#' @param where.ind index, can be numeric, or character
#' @return updated model
insert.one.line.code.estimable.model = function(esm,lhs,rhs,which.part = 'PK.PROGRAM', where.ind = 1) {
    if (! which.part %in% names(esm)) {
        stop(sprintf('The model spec doesnot contain [%s].', which.part))
    }
    if (is.character(where.ind)) {
        if (where.ind == '^') {
            where.ind = 1
        } else if (where.ind == '$') {
            where.ind = length(esm[[ which.part ]]) + 1
        } else {
            offset = 0
            if (substring(where.ind,1,1)=='^') {
            # insert before, default behavior
                where.ind = substring(where.ind,2)
                offset = 0 
            } else if (substring(where.ind,nchar(where.ind),nchar(where.ind))=='$') {
            # insert after, e.g. 'CL$' to insert the new one after 'CL'
                where.ind = substring(where.ind,1,nchar(where.ind)-1)
                offset = 1 
            }
            where.ind = match(where.ind, names(esm[[ which.part ]])) + offset
            if (is.na(where.ind)) where.ind = length(esm[[ which.part ]]) + 1
        }
    }
    esm[[ which.part ]] = insert.one.line.code.list(esm[[ which.part ]], lhs, rhs , index = where.ind)
    esm
}

#' edit code
#'
#' pay attention, not all code are related to a parameter/estimable.unit/random.observation
#' some are just expression for computation
#'
#' @param esm estimable.unit
#' @param lhs quoted, or character
#' @param rhs new definition
#' @param which.part which code, necessary
#' @return new mode
edit.one.line.code.estimable.model = function(esm, lhs, rhs, which.part='PK.PROGRAM' ) {
    if (is.character(lhs)) {
        lhs.ind = match(lhs, names(esm[[ which.part ]] ))
    } else if (is.symbol(lhs) || is.call(lhs)) {
        lhs.ind = match.symbolic(lhs, attr(esm[[ which.part ]] , 'symbolic.names'))
    } else if (is.numeric(lhs)) {
        lhs.ind = round(lhs)
        if (lhs.ind < 1 || lhr.ind > length(esm[[ which.part]])) {
            stop('lhs ind out of range')
        }
    } else {
        stop('lhs should be integer or character or expression')
    }
    if (is.na(lhs.ind)) {
        stop('lhs not found in the model')
    }
    # might affect PK.COV or OBS.COV
    vn = names(esm[[which.part]])[lhs.ind]
    tab = list.all.random.parameters.estimable.model(esm)
    # if vn is a random part, should remove all relation entries in PK.COV or OBS.COV
    # as treating the new added one are independent with old ones
    if (vn %in% tab$varname) {
        for(ind.by in unique(tab[ vn==tab$varname, 'indexed.by' ])){
            if (ind.by=='subject') {
                place = 'PK.COV'
            } else {
                place = 'OBS.COV'
            }
            if (length(esm[[ place ]]) > 0 ) {
                kept.ind = which(sapply(esm[[place]], function(x) !(vn %in% x$varnames)))
                if (length(kept.ind)==0) {
                    esm[[ place ]] = NULL
                } else {
                    esm[[ place ]] = esm[[place]][kept.ind]
                }
            }
        }
    }
    if (is.null(rhs)) {
        # remove it
        esm[[ which.part ]] = esm[[ which.part ]][ -lhs.ind ]
        if (!check.isvalid.basic.pk.parameters(esm,need.print=FALSE)) {
            warning('The new model is not complete caused by lacking necessary definition for basic PK paramters.')
        }
    } else {
        esm[[ which.part ]][[ lhs.ind ]] = rhs
    }
    esm
}

#' remove one line code in the model spec
#' @param esm estimable.model
#' @param lhs integer or character or quoted expression
#' @return updated estimable.model
remove.one.line.code.estimable.model = function(esm, lhs, which.part='PK.PROGRAM') {
    edit.one.line.code.estimable.model(esm,lhs,rhs=NULL,which.part)
}

#' specialize \code{estimable.unit} into simpler one
#'
#' The background of this design is following:
#' NONMEM has convenient \code{FIX} notation for remove one parameter from to be estimabled.
#' Our infrastructure does not allow this, but specialization plays a same meaning.
#' Specilize some parameters to a given value means degenerate the original model into simpler one(with less parameters).
#'
#' @param esm estimable.model
#' @param varnames which variable need to be specialize
#' @param all boolean specialize all estimable.unit to their initial guess
#' @return updated model
specialize.estimable.units.estimable.model = function(esm, varnames, all = FALSE) {
    if (missing(varnames) && !all) {
        stop('not all, and no varnames given for specifizing')
    }
    if (all) {
        varnames = unique(list.all.fixed.parameters.estimable.model(esm)$varname)
    }
    ind = match(varnames,names(esm$PK.PROGRAM))
    ind = ind[complete.cases(ind)]
    for(i in ind) {
        esm$PK.PROGRAM[[ i ]] = specialize.estimable.unit(esm$PK.PROGRAM[[ i ]])
    }
    ind = match(varnames,names(esm$OBS.PROGRAM))
    ind = ind[complete.cases(ind)]
    for(i in ind) {
        esm$OBS.PROGRAM[[ i ]] = specialize.estimable.unit(esm$OBS.PROGRAM[[ i ]])
    }
    esm
}

#' compare two estimable model
#'
#' The result is a table showing the type and typical value and initial estimations of basic PK parameters.
#' 
#' If two model are of different kinetic model (using different PREDPP subroutine), we do not real compare them, but output the subroutine they used.
#' @param m1 model spec 1
#' @param m2 model spec 2
#' @param model.names names used in output table
#' @return summary table
compare.estimable.model = function(m1, m2, model.names=c('Model1','Model2')){
    dyn1 = as.nonmem.subroutine.dynamics(m1$DYNAMICS)
    dyn2 = as.nonmem.subroutine.dynamics(m2$DYNAMICS)
    
    if (dyn1 != dyn2) {
        tab = data.frame(Model=model.names, DYNAMICS=c(dyn1,dyn2))
        return(tab)
    }
    # same dynamics, compare the definitions of basic parameters
    tab1 = cbind(Model=model.names[1],summary.estimable.model(m1))
    tab2 = cbind(Model=model.names[2],summary.estimable.model(m2))
    tab = NULL
    for(i in 1:NROW(tab1)) {
        tab = rbind(tab, tab1[i,])
        tab = rbind(tab, tab2[i,])
    }
    tab
}

#' summary one estimable model in a data.frame
#'
#' By default, the result returns only those information for basic PK parameters.
#' The option \code{all} can switch on the output for all variable appearing in the \code{PK.PROGRAM} and \code{OBS.PROGRAM}
#' @param esm estimable.model
#' @param all if all variable are output
#' @return data.frame
summary.estimable.model = function( esm, all=FALSE ) {
    tab0 = NULL
    if (all) {
        for(vn in names(esm$PK.PROGRAM)) {
            tab = summary.estimable.unit(esm$PK.PROGRAM[[ vn ]])
            tab0 = rbind(tab0,tab)
        }
        for(vn in names(esm$OBS.PROGRAM)) {
            tab = summary.estimable.unit(esm$OBS.PROGRAM[[ vn ]])
            tab0 = rbind(tab0,tab)
        }
        nms = c(names(esm$PK.PROGRAM),names(esm$OBS.PROGRAM))
    } else {
        nms = get.basic.pk.parameters(esm$DYNAMICS)
        for(vn in nms) {
            tab = summary.estimable.unit(esm$PK.PROGRAM[[ vn ]])
            tab0 = rbind(tab0,tab)
        }
    }
    cbind(Variable=nms,tab0)
}
