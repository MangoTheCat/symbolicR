## Canonical Estimable units
DATABASE.of.estimable.unit = list(
list(keywords=c('normal','basic'),
     constructor=function(init=0.01, lower=-Inf, upper=Inf, Var=1){
        m = morphism.concrete.new('identicalmapping', dim=1)
        m$domain=morphism.domain.one.dimensional(lower,upper)
        m = estimable.Unit.new(m)
        attr(m,'initial.estimation') = init
        attr(m,'distribution') = list(type='normal', Var=Var)
        m
     }),
list(keywords=c('log normal','basic'),
     constructor=function(init=0.01, lower=0.00001, upper=Inf, Var=1){
        m = morphism.concrete.new('log.change')
        m$domain=morphism.domain.one.dimensional(lower,upper)
        m = estimable.Unit.new(m)
        attr(m,'initial.estimation') = init
        attr(m,'distribution') = list(type='log normal', Var=Var)
        m
     }),
list(keywords=c('log normal','weigth'),
     constructor=function(init=0.01, lower=0.00001, upper=Inf, Var=1){
        m = morphism.unapply(quote(log(x * (WT/70))), quote(x))
        m$domain=morphism.domain.one.dimensional(lower,upper)
        m = estimable.Unit.new(m)
        attr(m,'initial.estimation') = init
        attr(m,'distribution') = list(type='log normal', Var=Var)
        m
     }),
list(keywords=c('normal','weigth'),
     constructor=function(init=0.01, lower=0.00001, upper=Inf, Var=1){
        m = morphism.unapply(quote( x * (WT/70)), quote(x))
        m$domain=morphism.domain.one.dimensional(lower,upper)
        m = estimable.Unit.new(m)
        attr(m,'initial.estimation') = init
        attr(m,'distribution') = list(type='normal', Var=Var)
        m
     }),
list(keywords=c('log normal','exponential'),
     constructor=function(init=0.01, lower=0.00001, upper=Inf, Var=1){
        m = morphism.unapply(quote(log(x * (WT/70)^0.75)), quote(x))
        m$domain=morphism.domain.one.dimensional(lower,upper)
        m = estimable.Unit.new(m)
        attr(m,'initial.estimation') = init
        attr(m,'distribution') = list(type='log normal', Var=Var)
        m
     }),
list(keywords=c('normal','exponential'),
     constructor=function(init=0.01, lower=0.00001, upper=Inf, Var=1){
        m = morphism.unapply(quote(x * (WT/70)^0.75), quote(x))
        m$domain=morphism.domain.one.dimensional(lower,upper)
        m = estimable.Unit.new(m)
        attr(m,'initial.estimation') = init
        attr(m,'distribution') = list(type='normal', Var=Var)
        m
     }),
list(keywords=c('fixed','affine','unconditioned'),
     constructor=function(A,b=NULL,est.bound.mat=NULL){
        m = convert.matrix.to.morphism(A,b)
        est = numeric(degree.of.freedom(m))
        if (!is.null(est.bound.mat)) {
            m$domain=as.morphism.domain.matrix(est.bound.mat)
            if ('Est' %in% colnames(est.bound.mat)) {
                est = est.bound.mat[,'Est']
            } else { # usually 2 is good
                est = est.bound.mat[,2]
            }
        }
        m = estimable.Unit.new(m)
        attr(m,'initial.estimation') = est
        m
     }),
list(keywords=c('fixed','affine','conditioned'),
     constructor=function(A, cond.var,cond.level,b=NULL,est.bound.mat=NULL){
        m = convert.matrix.to.morphism(A,b)
        est = numeric(degree.of.freedom(m))
        if (!is.null(est.bound.mat)) {
            m$domain=as.morphism.domain.matrix(est.bound.mat)
            if ('Est' %in% colnames(est.bound.mat)) {
                est = est.bound.mat[,'Est']
            } else { # usually 2 is good
                est = est.bound.mat[,2]
            }
        }
        m = combine.affinemorphism.factor(m, cond.var, as.character(cond.level))
        m = estimable.Unit.new(m)
        attr(m,'initial.estimation') = est
        m
     })
)

Canonical.Estimable.Unit.new = function(keywords, ...) {
    kws = lapply(DATABASE.of.estimable.unit, function(x) x$keywords)
    score = score.keywords.match(keywords, kws)
    ind = which.max(score)
    DATABASE.of.estimable.unit[[ind]]$constructor(...)
}
