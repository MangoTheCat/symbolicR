# Estimable unit has implicit meaning of ETA's , not EPS

#' lookup possible information of ETA[i]
#'
#' If available, find which the random part corresponding to original ETAs
#'
#' @param esu est estimable unit
#' @return NULL, or NAs or index
lookup.random.part.backref.estimable.unit = function(esu){
    if (is(esu,'random.disturbation')) {
        return(lookup.random.part.backref.random.disturbation(esu))
    }
    # estimable unit can be a morphism, random disturbation, call, name, or atomic
    d = attr(esu,'distribution')
    if (is.null(d)) return(NULL)
    bk = d$backref
    if (is.null(d)) return(NA)
    bk
}

#' correct log normal case when the value is defined by \code{x -> log(x)}
#' 
#' @param e estimable unit
#' @return fixed version
estimable.Unit.correct.lognormal = function(e){
    if(!is(e,'morphism') || 
        is.null(attr(e,'distribution')) || 
        attr(e,'distribution')$type!='log normal'
    ) return(e)
    if ( e$defn$type!='log.change' ) return(e)
    boundary = as.low.upper.list.morphism.domain(e)[[1]]
    if (boundary[1] > 0 ) return(e)
    # usually, we assume RealNumber for it, so if the initial.estimation is positive, we are not willing to "fix" it
    if (boundary[2] >0 && attr(e,'initial.estimation')>0) return(e)
    fixed.syms = lapply(1:degree.of.freedom(e), function(x) CONS('[',quote(FIXED.dummy),as.numeric(x)))
    fixed = simplify.collect.log.exp(CONS('exp',instantiate.morphism(e, fixed.syms)))
    m = random.disturbation.new(CONS('*', fixed, quote( exp(RANDOM[1]) )), 
            fixed = fixed.syms,
            distribution = distribution.common.new(cov.=attr(e,'distribution')$Var))
    m$defn$domain = CONS('*',e$domain,quote(RealNumber^1))
    attr(m,'initial.estimation') = attr(e,'initial.estimation')
    # put back both the fixed part and the random part
    backref = c(backref.parameter.of.morphism(e), 
        lookup.random.part.backref.estimable.unit(e))
    attr(m,'backref') = backref
    class(m) = c('estimable.unit', class(m))
    m
}

#' Promotion a expression to estimable unit
#'
#' promote a morphism or a expression to a estimable.unit
#' However, 'name' is special in R, if an class attribute is added to a name, R will be forced to evaluate it,
#' need to avoid this.
#'
#' @param m morphism as estimable unit
#' @return m with class attribute updated
estimable.Unit.new = function(m) {
    # R will treat 'name' quite different from 'call'
    # a name will be evaluated by R when passed to a function if an extra class attribute is added
    # and other attributes will be lost
    if (is(m,'name')) return(m)
    if (!is(m,'estimable.unit')) class(m) = c('estimable.unit', class(m))
    m
}

#' S3 version formater
print.estimable.unit = function(m) {
    indent.n = function(ss,n) {
       if (n==0) return(ss)
       indent = paste(rep(' ',n),collapse='')
       ss = gsub('^', indent, ss)
       gsub('\n', sprintf('\n%s',indent), ss)
    }
    if (is(m,'morphism')){
        if (m$defn$type=='Constant') {
        # we don't want to print a Constant as * |---> constant
            m0 = paste(deparse(instantiate.morphism(m)),collapse='')
        } else {
            m0 = print.morphism0(m)
        }
        m0 = indent.n(m0, 4)
        m0 = sprintf('Estimable Unit ~ \n%s', m0)
        if(!is.null(tmp <- attr(m, 'initial.estimation'))) {
            m0 = sprintf('%s\nInitial guess of parameters: (%s)',m0, paste(tmp, collapse=','))
        }
    } else if (is(m,'random.disturbation')) {
        m0 = paste(capture.output(print.random.disturbation(m)),collapse='\n')
        m0 = indent.n(m0, 4)
        m0 = sprintf('Estimable Unit ~ \n%s', m0)
        if(!is.null(tmp <- attr(m, 'initial.estimation'))) {
            m0 = sprintf('%s\nInitial guess of FIXED parameters: (%s)',m0, paste(tmp, collapse=','))
        }
    } else if (is(m,'name') && !is.null(attr(m,'distribution'))) {
        m0 = paste(deparse(m), collapse=',')
        m0 = indent.n(m0, 4)
        m0 = sprintf('Estimable Unit ~ \n%s', m0)
    } else if (is(m,'call') && !is.null(attr(m,'distribution'))) {
        m0 = paste(deparse(m), collapse=',')
        m0 = indent.n(m0, 4)
        m0 = sprintf('Estimable Unit ~ \n%s', m0)
    } else if (is.atomic(m)) {
        m0 = sprintf('Estimable Unit ~ \n%s', indent.n(m,4))
    } else {
        stop('Can not print ', paste(deparse(m),collapse=''),' as estimable unit.\n')
    }
    if(!is.null(tmp <- attr(m, 'distribution'))){
        m0 = sprintf('%s\n%s',m0,
                sprintf('random effect of type %s, with variance: %s', 
                    tmp$type, paste(deparse(tmp$Var),collapse='')))
    }
    cat(m0, '\n')
}

#' update the FIXed parameters for estimable unit
#'
#' This function will check validity of the input value.
#' @param m estimable unit
#' @return updated estimable unit
update.initial.estimation.estimable.unit = function(m, new.val) {
    stopifnot(is(m,'estimable.unit'))
    if (is(m,'morphism')) {
        N = degree.of.freedom(m)
        if ((m$defn$type!='Constant') && length(new.val) != N) {
            stop(sprintf('The new.val has length %s while the morphism expect %s values.', 
                length(new.val), N ))
        }
        if (!does.morphism.domain.contain(m,new.val)) {
            warning('The new value is not in the domain of the morphism. 
            We suggest update the domain accordingly.')
        }
        attr(m,'initial.estimation') = new.val
        return(m)
    } else if (is(m,'random.disturbation')) {
        N = N.fixed.effects(m)
        if (length(new.val) != N) {
            stop(sprintf('The new.val has length %s while the random disturbation expect %s 
                number of Fixed parameter(s).', length(new.val), N ))
        }
        domain = random.disturbation.extract.fixed.domain(m)
        if (!does.morphism.domain.contain(domain,new.val)) {
            warning('The new value is not in the domain of fix effect of the disturbation. 
            We suggest update the domain of fixed effect accordingly.')
        }
        attr(m,'initial.estimation') = new.val
        return(m)
    }
    stop('Seems this estimable unit does not contain any initial estimation for Fix Part.')
}

#' update the fixed domain of a estimable.unit if there is any 
#'
#' Estimable unit can be implemented by either morphism or random disturbation. 
#' But the word domain here always related to the fixed part. For morphism, random effect is just an attribute,
#' for random disturbation, domain means the domain of fixed parameters.
#' @param m estimable unit
#' @param domain domain expression
#' @return updated estimable unit
update.domain.estimable.unit = function(m , domain) {
    stopifnot(is(m,'estimable.unit'))
    if (!is.null(vec<-attr(m,'initial.estimation'))){
        if (!does.morphism.domain.contain(domain, vec)) {
            warning('The new domain does not contain current initial estimation. 
            We suggest update the initial estimation accordingly.')
        }
    }
    if (is(m,'morphism')){
        if (m$defn$type=='Constant') {
            m$domain = domain
        } else if ((dm<-degree.of.freedom(m)) == (dd<-dimension.of.domain(domain))) {
            m$domain = domain
        } else {
            stop(sprintf('Previous dimension is %s, while the new domain has dimension %s', dm, dd ))
        }
        return(m)
    } else if (is(m,'random.disturbation')) {
        dm = N.fixed.effects(m)
        dd = dimension.of.domain(domain)
        if (dm!=dd){
            stop(sprintf('Previous dimension of fixed part is %s, while the new domain has dimension %s', dm, dd ))
        }
        return(random.disturbation.update.fixed.domain(m, domain))
    }
    stop('Seems this estimable unit does not contain any FIX effects, hence no domain can be updated.')
}

#' convert estimable unit to a simpler one without random effect
#'
#' Atomic or Constant morphism are simplified into expreesion or atomic, 
#  may not remain an estimable unit.
#'
#' @param m estimable.unit
#' @return estimable.unit with random effect removed, or simply an expression
remove.random.effect.estimable.unit = function(m) {
    stopifnot(is(m,'estimable.unit'))
    if (is(m,'morphism')) {
        # For log normal, we treat it specially
        if (!is.null(distr<-attr(m,'distribution'))) {
            if (distr$type=='log normal') {
                if (is(m,'morphism')) {
                    m1 = morphism.composition.2(morphism.concrete.new('exponential.change'), m)
                    m1 = estimable.Unit.new(m1)
                    attr(m1,'distribution') = NULL
                    attr(m1,'backref') = NULL
                    attr(m1,'initial.estimation') = attr(m,'initial.estimation')
                } else if (is.atomic(m)) {
                    m1 = exp(m)
                } else if (is.symbol(m) || is.call(m)) {
                    m1 = simplify.collect.log.exp(CONS('exp',m))
                }
                return(m1)
            }
        }
        attr(m,'distribution') = NULL
        if (m$defn$type=='Constant') m = instantiate.morphism(m)
        return(m)
    }
    if (is(m,'random.disturbation')) {
        m1 = remove.random.effect.random.disturbation(m)
        # only the morphism need to be promote into estimable unit
        if (is(m1,'morphism')) {
            m1 = estimable.Unit.new(m1)
            # If it's estimable.unit
            # should also keep the the initial.estimation
            attr(m1,'initial.estimation') = attr(m,'initial.estimation')
        }
        return(m1)
    }
    # May also be atomic
    if (is.atomic(m)) {
        attr(m,'distrubation') = NULL
        # it should not remain a estimable unit
        m = unclass(m)
        return(m)
    }
    stop(sprintf('Can not remove the distrubation of class [%s]', paste(class(m),collapse=',')))
}

#' remove the random part of a random disturbation
#'
#' The origianl random disturbation may degenerate into a morphism if there is non-trivial fixed parameters.
#' If no fixed parameter exists, then the original disturbation may even degenerated into an expression or a symbol.
#'
#' @param m random.disturbation
#' @return disturbation or morphism
remove.random.effect.random.disturbation = function(m) {
    stopifnot(is(m,'random.disturbation'))
    if (N.random.effects(m)>0) {
        rd = as.symbolic.vector(rep(0, N.random.effects(m)))
    } else {
        rd = NULL
    }
    if (N.fixed.effects(m)==0) {
        # then it degenerate to expression
        re = simplify.2(disturbation.apply(m, random=rd))
    } else {
        # degenerate into morphism
        fixed = lapply(1:N.fixed.effects(m), function(ii) CONS('[', quote(Rtmpo5mJgn), as.numeric(ii)))
        re = morphism.unapply(simplify.2(disturbation.apply(m, fixed=fixed, random=rd)), fixed)
        # should keep the domain
        re$domain = random.disturbation.extract.fixed.domain(m)
    }
    if (is(m,'multiple.dependent.variable')) {
        if (is(re,'Rmatrix')) {
            re = asCMatrix.RMatrix(re)
        }
        if (is(re,'CMatrix')) {
            re = list(levels=re)
        }
        attr(re,'selecting.variable') = attr(m,'selecting.variable')
        attr(re,'levels') = attr(m,'levels')
        class(re) = union('multiple.dependent.variable', class(re) )
    }
    re
}

#' change the distribution attribute in a estimable unit
#' @param m estimable.unit
#' @param distr distribution, type and Var
#' @return updated estimable.unit
update.distribution.estimable.unit = function(m, distr) {
    stopifnot(is(m,'estimable.unit'))
    if(is(m,'morphism') || is.atomic(m)) {
        attr(m,'distribution') = distr
        return(m)
    }
    if(is(m,'random.disturbation')) {
        if (!is(distr, 'distribution')) 
            distr = convert.distrib.description.to.random.disturbation(type=distr$type, Var=distr$Var)
        m$distribution = distr
        return(m)
    }
    stop(sprintf('Not able to update the distribution of the estimable unit of class [%s].',
        paste(class(m),collapse=',')))
}

#' update the covariance matrix of a estimable.unit
#'
#' Similar to \code{\link{update.distribution.estimable.unit}}, but
#' more common user story is we only want to change the covariance.
#' 
#' If the previous estimable unit has no random effect, we simple add one multi-normal distribution from give covariance matrix.
#'
#' @param m estimable.unit
#' @param Cov Variance or Covariance
#' @return updated estimable.unit
update.cov.estimable.unit = function(m, Cov) {
    stopifnot(is(m,'estimable.unit'))
    if (any(eigen(Cov)$values <=0)) {
        warning('The covariance matrix is not positive definite.')
    }
    if(is(m,'morphism') || is.atomic(m)) {
        old.distr = attr(m,'distribution')
        if (is.null(old.distr)) {
            distr = convert.distrib.description.to.random.disturbation(type='normal', Var = Cov)
        } else {
            distr = old.distr
            distr$Var = Cov
        }
        m = update.distribution.estimable.unit(m, distr)
        return(m)
    }
    if(is(m,'random.disturbation')) {
        old.distr = attr(m,'distribution')
        if (is.null(old.distr)) {
            distr = distribution.common.new(dim = NROW(Cov), cov.=Cov)
        } else {
            distr = update.distribution.covariance(old.distr, Cov)
        }
        if (N.random.effects(m)!=NROW(Cov)) {
            stop('Number of Random effect and order of Cov matrix dismatch')
        }
        m = update.distribution.estimable.unit(m, distr)
        return(m)
    }
    stop('Not able to update the distribution of the estimable unit.')
}

#' update the distribution for a random.observation object
#'
#' By design, we record expression with \code{EPS[i]} as random observation.
#'
#' @param m random observation 
#' @param distr distribution object or list of type and Var
#' @return updated random observation
update.distribution.random.observation = function(m , distr) {
    stopifnot(is(m,'random.observation'))
    if (is(m,'random.disturbation')) {
        if (!is(distr, 'distribution')) 
            distr = convert.distrib.description.to.random.disturbation(type=distr$type, Var=distr$Var)
        m$distribution = distr
        return(m)
    }
    stop('Not able to update the distribution of the random observation.')
}

#' update the covariance of distribution for a random.observation object
#' 
#' More common situation is only change the covariance, not the whole type.
#' @param m random observation 
#' @param Cov  covariance matrix
#' @return updated random observation
update.cov.random.observation = function(m, Cov) {
    stopifnot(is(m,'random.observation'))
    if (any(eigen(Cov)$values <=0)) {
        warning('The covariance matrix is not positive definite.')
    }
    if(is(m,'random.disturbation')) {
        old.distr = attr(m,'distribution')
        if (is.null(old.distr)) {
            distr = distribution.common.new(dim = NROW(Cov), cov.=Cov)
        } else {
            distr = update.distribution.covariance(old.distr, Cov)
        }
        if (N.random.effects(m)!=NROW(Cov)) {
            stop('Number of Random effect and order of Cov matrix dismatch')
        }
        m = update.distribution.random.observation(m, distr)
        return(m)
    }
    stop('Not able to update the covariance of the random observation.')
}

#' get the initial estimation from a estimable unit
#' @param m estimable unit
#' @return value of initial estimation
get.initial.estimation.estimable.unit = function(m) {
    stopifnot(is(m,'estimable.unit'))
    attr(m,'initial.estimation')
}

#' extract the domain part of a estimable.unit
#'
#' Only fixed part has domain.
#' @param m
#' @return domain expression
get.domain.estimable.unit = function(m) {
    stopifnot(is(m,'estimable.unit'))
    if (is(m,'morphism')) {
        return(m$domain)
    } else if (is(m,'random.disturbation')) {
        domain = random.disturbation.extract.fixed.domain(m)
        return(domain)
    }
    stop('This estimable unit does not have fixed part (hence no domain).')
}

#' the lower est upper triple 
#'
#' Usually we use a triple (three column data.frame) for storing the lower bound, estimation value and upper bound.
#' After extracting the comain of estimable unit, we convert it to the triple representation.
#' @param m estimable.unit
#' @return triple
get.triple.estimable.unit = function(m) {
    stopifnot(is(m,'estimable.unit'))
    est = attr(m,'initial.estimation')
    if (is(m,'morphism')) {
        domain = m$domain
    } else if (is(m,'random.disturbation')) {
        domain = random.disturbation.extract.fixed.domain(m)
    } else {
        stop('This estimable unit does not have fixed part (hence no domain).')
    }
    domain = as.low.upper.list.morphism.domain(domain)
    data.frame(Lower = sapply(domain, el, 1), Est = est, Upper = sapply(domain, el, 2))
}

#' Specialize a estimable unit
#'
#' Change estimable.unit to expression or simpler estimable.unit by specializing parameter values
#' or just substitute the initial estimations into the morphism or random disturbation.
#'
#' @param m estimable.unit
#' @param val specialize on which argument
#' @param ind optional index
#' @param simplify if simplify to simple atomic or expression, only possible give full value for vals
#' @return new estimable.unit
specialize.estimable.unit = function(m, val=get.initial.estimation.estimable.unit(m), 
                                    ind=seq_along(as.symbolic.vector(val)), simplify=FALSE) {
    if(!is(m,'estimable.unit') || length(val)<1) return(m)
    if(is(m,'morphism')) {
        N = degree.of.freedom(m)
        if (simplify && (length(ind)==N || N==0)) {
            if (N==0) {
                tmp = instantiate.morphism(m)
            } else {
                # protect it against atomic input
                params = as.symbolic.vector(val)
                params[ind] = val
                tmp = instantiate.morphism(m, params)
            }
            if (!is.null(distr<-attr(m,'distribution'))) {
                if (distr$type=='log normal') {
                    tmp = simplify.collect.log.exp(CONS('exp', tmp))
                }
            }
            return(tmp)
        }
        if (N==0) return(m)
        if (length(ind) > N) {
            warning('More indices supplied than number of fixed parameters')
            ind = ind[1:N]
        }
        m1 = estimable.Unit.new(morphism.curry.point(m, val, ind))
        if (!is.null(attr(m,'initial.estimation'))) {
            attr(m1,'initial.estimation') = attr(m, 'initial.estimation')[-(1:length(val))]
            if (length(attr(m1,'initial.estimation'))==0) attr(m1,'initial.estimation')=NULL
        }
        attr(m1,'distribution') = attr(m, 'distribution')
        if (m1$defn$type=='Constant' && 
            is.null(attr(m1,'initial.estimation')) && 
            is.null(attr(m1,'distribution'))){
            # simplify it to expression
            m1 = instantiate.morphism(m1)
        }
        return(m1)
    } else if (is(m,'random.disturbation')) {
        N = N.fixed.effects(m)
        if (simplify && (length(ind)==N || N==0)) {
            if (N==0) {
                tmp = disturbation.apply(m, random=0)
            } else {
                params = as.symbolic.vector(val)
                params[ind] = val
                tmp = disturbation.apply(m, fixed=params, random=0)
            }
            return(tmp)
        }
        if (N==0) return(m)
        if (length(ind) > N) {
            warning('More indices supplied than number of fixed parameters')
            ind = ind[1:N]
        }
        m1 = morphism.curry.point(m$defn, val, ind)
        rd = estimable.Unit.new(m)
        rd$defn  = m1
        rd$FIXED = m$FIXED[ -ind ]
        # backref is not longer valid, remove it
        attr(rd,'backref') = NULL
        if (!is.null(attr(m,'initial.estimation'))) {
            attr(rd,'initial.estimation') = attr(m, 'initial.estimation')[-(1:length(val))]
            if (length(attr(rd,'initial.estimation'))==0) attr(rd,'initial.estimation')=NULL
        }
        if (N.random.effects(rd)==0 && N.fixed.effects(rd)==0) {
            rd = disturbation.apply(rd)
        }
        return(rd)
    }
    if (simplify) {
        return(m)
    }
    warning('Not curryable estimable.unit, return as is')
    m
}

#' specialize a random disturbation
#' @param m random.disturbation
#' @param val values for fixed paramters
#' @param ind optional index 
#' @param simplify when all val given, simpify to typical value
#' @return random.disturbation or expression
specialize.random.disturbation = function(m, val, ind=seq_along(as.symbolic.vector(val)),
                                        simplify = FALSE) {
    if (!is(m,'random.disturbation')) return(m)
    if (missing(val)) {
        stop('To specialize one random disturbation, you need to provide val arguments.')
    }
    N = N.fixed.effects(m)
    if (simplify && (length(ind)==N || N==0)) {
        if (N==0) {
            tmp = disturbation.apply(m, random=0)
        } else {
            params = val
            params[ind] = val
            tmp = disturbation.apply(m, fixed=params, random=0)
        }
        return(tmp)
    }
    if (N==0) return(m)
    if (length(ind) > N) {
        warning('More indices supplied than number of fixed parameters')
        ind = ind[1:N]
    }
    m1 = morphism.curry.point(m$defn, val, ind)
    rd = m
    rd$defn  = m1
    rd$FIXED = m$FIXED[ -ind ]
    # backref is not longer valid, remove it
    attr(rd,'backref') = NULL
    if (N.random.effects(rd)==0 && N.fixed.effects(rd)==0) {
        rd = disturbation.apply(rd)
    }
    return(rd)
}

#' description for a estimable unit
#'
#' similar to print, but results a data.frame
#' if input is not an estimable unit, return 0-rows data.frame
#'
#' @param m estimable unit
#' @return data.frame
summary.estimable.unit = function(m) {
    types = class(m)
    morphism.type = NA
    if (!is(m,'estimable.unit')) {
        type = types[1]
        defined.by = paste(setdiff(types,type),collapse=',')
        N.fixed = NA
        fixed.domain = NA
        initial.estimation = NA
        variance = NA
        if (is.call(m)) {
            typical.val = paste(deparse(m),collapse='')
        } else if (is.atomic(m)) {
            typical.val = paste(m,collapse=',')
        } else if (is.symbol(m)) {
            typical.val = paste(deparse(m),collapse=',')
        } else {
            typical.val = NA
        }
    } else {
        # Export the result data.frame
        type = 'estimable.unit'
        if (is(m, 'morphism')) {
            defined.by = 'morphism'
            morphism.type = m$defn$type
            N.fixed = degree.of.freedom(m)
            fixed.domain = paste(deparse(m$domain),collapse='')
            initial.estimation = paste(attr(m,'initial.estimation'),collapse=',')
            if (is.null(attr(m,'distribution'))) {
                variance = 'NONE'
            } else {
                variance = sprintf('%s(Var=%s)', 
                    attr(m,'distribution')$type,
                    attr(m,'distribution')$Var)
            }
            typical.val = paste(deparse(specialize.estimable.unit(m,simplify=TRUE)),collapse='')
        } else if (is(m, 'random.disturbation')) {
            defined.by = 'random.disturbation'
            N.fixed = N.fixed.effects(m)
            if (N.fixed>0) {
                fixed.domain = paste(deparse(random.disturbation.extract.fixed.domain(m)),collapse='')
            } else {
                fixed.domain = NA
            }
            initial.estimation = paste(attr(m,'initial.estimation'),collapse=',')
            variance = print.distribution.common0(m$distribution)
            typical.val = paste(deparse(specialize.random.disturbation(m, simplify=TRUE)),collapse='')
        } else {
        # atomic or call or symbol
            defined.by = paste(setdiff(types, type),collapse=',')
            N.fixed = 0
            fixed.domain = NA
            initial.estimation = NA
            if (is.null(attr(m,'distribution'))) {
                variance = 'NONE'
            } else {
                variance = sprintf('%s (Var=%s)', 
                    attr(m,'distribution')$type,
                    attr(m,'distribution')$Var)
            }
            if (is.call(m)) {
                typical.val = paste(deparse(m),collapse='')
            } else if (is.atomic(m)) {
                typical.val = paste(m,collapse=',')
            } else if (is.symbol(m)) {
                typical.val = paste(deparse(m),collapse=',')
            } else {
                typical.val = NA
            }
        }
    }
    data.frame(
               Typical.val = typical.val,
               Estimable= (type == 'estimable.unit'), 
               IsMorphism=(defined.by=='morphism'),
#               Morphism.type=morphism.type,
               N.fixed=N.fixed,
               Domain.of.fixed.parameter=fixed.domain,
               Initial.Est=initial.estimation,
               Random.effects = variance
   )
}

#' compare two estimable units
#' @param m1 first  estimable unit
#' @param m2 second estimable unit
#' @return data.frame
compare.estimable.unit = function(m1,m2){
    tab1 = summary.estimable.unit(m1)
    tab2 = summary.estimable.unit(m2)
    rbind(tab1,tab2)
}
