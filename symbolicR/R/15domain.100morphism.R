# explain a expression as morphism
morphism.concrete.database = list(
list(type='identicalmapping',
     constructor=function(dim, ...){
            dim = as.integer(dim)
            stopifnot(dim>0)
            list(domain=CONS('^',quote(RealNumber) , as.numeric(dim)),
                  image =CONS('^',quote(RealNumber) , as.numeric(dim)),
                  defn  =list(type='identicalmapping', ... )) },
     instantiator=function(m, parameter.pool){
     # we following the conversion that a morphism to 1 dimensional space should not be a list format
            if ((dims<-degree.of.freedom(m))==1) return(parameter.pool[[1]])
            parameter.pool[1:dims]
     }),
list(type='exponential.change',
     constructor=function(...){
            list(domain=quote(RealNumber^1),
              image =quote(PositiveRealNumber^1),
              defn  =list(type='exponential.change', ...)) },
     instantiator=function(m, parameter.pool) {
            CONS('exp', parameter.pool[[1]])
     }),
list(type='log.change',
     constructor=function(...){
            list(domain=quote(PositiveRealNumber^1),
              image =quote(RealNumber^1),
              defn  =list(type='log.change', ...)) },
     instantiator=function(m, parameter.pool) {
              CONS('log', parameter.pool[[1]])
     }),
list(type='mirror.reflection',
     constructor=function(...){
            list(domain=quote(RealNumber^1),
              image =quote(RealNumber^1),
              defn  =list(type='mirror.reflection', ...)) },
     instantiator=function(m, parameter.pool) {
              simplify(CONS('-', parameter.pool[[1]]))
     }),
list(type='inverse.logit.change',
     constructor=function(...){
            list(domain=quote(RealNumber^1),
                  image =quote(UnitInterval^1),
                  defn  =list(type='inverse.logit.change', ...))},
     instantiator=function(m, parameter.pool){
            CONS('/',1 , CONS('+', 1, CONS('exp', CONS('-', parameter.pool[[1]]))))
     }),
list(type='logit.change',
     constructor=function(...){
            list(domain=quote(RealNumber^1),
                  image =quote(UnitInterval^1),
                  defn  =list(type='logit.change', ...))},
     instantiator=function(m, parameter.pool){
            CONS('log', CONS('/',parameter.pool[[1]] , CONS('-', 1, parameter.pool[[1]])))
     }),
list(type='Constant',
     constructor=function(val,...) {
         domain = list(...)$domain
         val = as.symbolic.vector(val)
         N.val = length(val)
         if (N.val==1) val=val[[1]]
         re = list(domain=quote(AnyDomain^n),
                image=CONS('^',quote(RealNumber), as.numeric(N.val)),
                defn =list(type='Constant', defn=val))
         if (!is.null(domain)) re$domain = domain
         re },
     instantiator=function(m, parameter.pool=NULL){
         dims = dimension.of.image(m)
         if (is.atomic(m$defn$defn) && length(m$defn$defn)>1) {
             re = as.list(m$defn$defn)
         } else {
             if(dims==1 && is.list(m$defn$defn))
              re = m$defn$defn[[1]]
             else
              re = m$defn$defn   
         }
         re
     }),
list(type='Product', # prod
     constructor=function(dim,...){
         list(domain=CONS('^',quote(RealNumber),dim),
                image=quote(RealNumber^1),
                defn =list(type='Product',...)) },
     instantiator=function(m, parameter.pool) {
         expression.from.list.of.product( parameter.pool )
     }),
list(type='Sum', # sum
     constructor=function(dim,...){
         list(domain=CONS('^',quote(RealNumber),dim),
                image=quote(RealNumber^1),
                defn =list(type='Sum',...)) },
     instantiator=function(m, parameter.pool) {
         expression.from.list.of.sum( parameter.pool )
     }),
list(type='AffineMapping',
     constructor=function(affine.mapping,...){
         list(domain=CONS('^',quote(RealNumber), as.numeric(length(affine.mapping$A[[1]]))),
              image =CONS('^',quote(RealNumber), as.numeric(length(affine.mapping$A))),
              defn  =list(type='AffineMapping', defn=affine.mapping, ... )) },
     instantiator=function(m, parameter.pool) {
         affine = m$defn$defn
         re = symbolic.affine.mapping(affine$A, parameter.pool, affine$b)
         if (is.list(re) && dimension.of.image(m)==1) re=re[[1]]
         re
     }),
list(type='affine.mapping.factor.effect',
     constructor=function(mapping,...){
         list(domain=CONS('^',quote(RealNumber), as.numeric(length(mapping$affine.mapping$A[[1]]))),
              image =quote(RealNumber^1),
              defn  =list(type='affine.mapping.factor.effect', 
                          defn=mapping,...))},
     instantiator=function(m, parameter.pool) {
         af = m$defn$defn$affine.mapping
         y = m$defn$defn
         v1 = vector(length(y$level.name),mode='list')
         for(i in seq_along(y$level.name)){
             # convert levels back to numeric if possible
             tmp = suppressWarnings(as.numeric(y$level.name[i]))
             if (is.na(tmp)) tmp = y$level.name[i]
             v1[[i]]=CONS('==',as.symbol(y$var.name), tmp)
         }
         v2 = symbolic.affine.mapping(af$A, parameter.pool, af$b)
#        a little different between following two, first will return 0 if no condition success
#        second will return NA is no condition success
#         symbolic.inner.product(v1,v2)
         symbolic.inner.product.logical(v1,v2)
     }),
list(type='Product.of.Morphisms',
     constructor=function(list.of.morphism, ...) {
         domain=list.of.morphism[[1]]$domain
         if (length(list.of.morphism)>1) {
             for(i in 2:length(list.of.morphism)) {
                 domain = simplify.cartesian.product(CONS('*', domain, list.of.morphism[[i]]$domain))
             }
         }
         list(domain=domain,
              image =quote(RealNumber^1),
              defn  =list(type='Product.of.Morphisms', 
                          defn=list.of.morphism,...))},
     instantiator=function(m, parameter.pool){
         lms = m$defn$defn
         ind = 0
         l = vector(length(lms),mode='list')
         for(i in seq_along(lms)) {
             dg = degree.of.freedom(lms[[i]])
             l[[i]]=instantiate.morphism(lms[[i]], parameter.pool[ (ind+1):(ind+dg) ])
             ind = ind + dg
         }
         simplify.2(expression.from.list.of.product(l))
     }),
list(type='Sum.of.Morphisms',
     constructor=function(list.of.morphism, ...) {
         domain=list.of.morphism[[1]]$domain
         if (length(list.of.morphism)>1) {
             for(i in 2:length(list.of.morphism)) {
                 domain = simplify.cartesian.product(CONS('*', domain, list.of.morphism[[i]]$domain))
             }
         }
         list(domain=domain,
              image =quote(RealNumber^1),
              defn  =list(type='Sum.of.Morphisms', 
                          defn=list.of.morphism,...))},
     instantiator=function(m, parameter.pool){
         lms = m$defn$defn
         ind = 0
         l = vector(length(lms),mode='list')
         for(i in seq_along(lms)) {
             dg = degree.of.freedom(lms[[i]])
             l[[i]]=instantiate.morphism(lms[[i]], parameter.pool[ (ind+1):(ind+dg) ])
             ind = ind + dg
         }
         simplify.2(expression.from.list.of.sum(l))
     }),
list(type='general.pointwise', # general pointwise mapping
     constructor=function(domain,image,pts.from,pts.to,...) {
         list(domain=domain,
                image=image,
                defn =list(type='general.pointwise',
                           defn=list(pts.from=pts.from,pts.to=pts.to),
                           ...)) },
     instantiator=function(m, parameter.pool) {
        mapping = m$defn$defn
        exps = mapping$pts.to
        if (!is.list(exps)) exps = list(exps)
        args = mapping$pts.from
        if (!is.list(args)) args = list(args)
        if (!is.list(parameter.pool)) parameter.pool = list(parameter.pool)
        if (length(parameter.pool)!=length(args)) stop('Parameter arity not equal to definition')
        for(i0 in seq_along(exps)) {
            for(i1 in seq_along(args)) {
                exps[[i0]] = symbolic.gsub( args[[i1]], parameter.pool[[i1]], exps[[i0]] )
            }
        }
        if (dimension.of.image(m)==1) exps = exps[[1]]
        exps
     })
)

#' morphism constructor
#' @param type type of morphism
#' @return the morphism object
morphism.concrete.new = function(type, ...) {
    ind = match(type, sapply(morphism.concrete.database, function(x) x$type))
    if (is.na(ind)) {
        stop(sprintf('Do not know how to construct %s type morphism', type))
    } else {
        re = morphism.concrete.database[[ ind ]]$constructor(...)
    }
    class(re) = c('morphism',class(re))
    re
}

#' morphism instantiator
#' @param m morphism
#' @param parameter.pool to which symbols we instantiate it
#' @return an expression or a set of expressions
morphism.concrete.instantiate = function(m, parameter.pool) {
    stopifnot(is(m,'morphism'))
    ind = match(m$defn$type, sapply(morphism.concrete.database, function(x) x$type))
    if (is.na(ind)) {
        return(NULL)
    } else {
        re = morphism.concrete.database[[ ind ]]$instantiator(m, parameter.pool)
    }
    re
}

#' for one or more expressions
#' find in \code{available.symbol.table} an not used one, which can represent the formal arguments
#' @param es one or more expressions
#' @param available.symbol.table characters which looks better for printing
#' @return first found symbol, or NULL if not found
generate.unused.symbol.in.expressions = function(es, 
    available.symbol.table=c('x','y','z','w','a','alpha','beta','gamma')) {
    if (!is.list(es)) {
        es = list(es)
    }
    for(sym0 in available.symbol.table){
        sym = as.symbol(sym0)
        has.sym = FALSE
        for(e0 in es) {
            if (symbolic.grep(CONS('[',sym, quote('?a'(any.ind))), e0)) {
                has.sym = TRUE
                break()
            }
        }
        if (!has.sym) return(sym)
    }
    NULL
}

#' \code{unapply( x+y+z , (x,y))} will result the following morphism:
#' \code{(x,y) -> x+y +z}
#' @param e expression
#' @param tuple the arguments
#' @param ... other arguments passed into general.pointwise morphism constructor
#' @return morphism
morphism.unapply = function( e, tuple, ... ) {
    if (is.call(tuple) && tuple[[1]]=='TUPLE') {
        tuple = as.list(tuple)[-1]
    }
    if (!is.list(tuple)) {
        domain = CONS('^', quote(RealNumber), 1)
        args = list(tuple)
    } else {
        if (length(tuple)==0) {
        # Constant
            return(morphism.concrete.new('Constant', e))
        }
        domain = CONS('^',quote(RealNumber), as.numeric(length(tuple)))
        args = tuple
    }
    if (!is.list(e)) {
        image = CONS('^',quote(RealNumber), 1)
        eqns = list(e)
    } else {
        image = CONS('^',quote(RealNumber),as.numeric(length(e)))
        eqns = e
    }
    # make it pretty as x[1],x[2],...,x[n]
    if (!is.null(new.sym<-generate.unused.symbol.in.expressions(eqns))) {
        for(i0 in seq_along(eqns)) {
            for(i1 in seq_along(args)) {
                eqns[[i0]] = symbolic.gsub( args[[i1]], CONS('[',new.sym,as.numeric(i1)), eqns[[i0]] )
            }
        }
        #
        if (!is.list(tuple)) {
            tuple = CONS('[',new.sym,1)
        } else {
            tuple = lapply(seq_along(args), function(ii) CONS('[',new.sym,as.numeric(ii)))
        }
        #
        if (!is.list(e)) {
            e = eqns[[1]]
        } else {
            e = eqns
        }
    }
    #
    morphism.concrete.new('general.pointwise', domain=domain, image=image, pts.from=tuple, pts.to=e, ... )
}

#' simplify compose two morphisms
#' do not check if compatible
#' @param m1 morphism1
#' @param m2 morphism2
#' @return a new morphism as \code{ m1 o m2 }, i.e. first \code{m2} then followed by \code{m1}
morphism.composition.2 = function(m1,m2) {
   # special deal with exp(log(*)) or log(exp(*))
   if ((m1$defn$type=='exponential.change' && m2$defn$type=='log.change') ||
       (m1$defn$type=='log.change' && m2$defn$type=='exponential.change')) {
         re = m2
         re$defn$type='identicalmapping'
         return(re)
   }
   re = list(domain =m2$domain,
              image  =m1$image,
              defn   =list(type='morphism.composition',
                           defn=list(m1,m2)))
   class(re) = c('morphism', class(re))
   re
}

#' other handy interface for Product of morphisms
#' @param ml morphism list
#' @return morphism being the chained up morphisms in \code{ml}
morphism.composition.list = function(ml) {
   N = length(ml)
   if (N==0) return(NULL)
   for(i in 1:N) {
       if (is.atomic(ml[[i]]) || is.symbol(ml[[i]]) || ((is.language(ml[[i]]) && is.call(ml[[i]])))) {
            ml[[i]] = morphism.concrete.new('Constant', val=ml[[i]])
       }
       stopifnot(is(ml[[i]],'morphism'))
       if (ml[[i]]$defn$type=='Constant') {
           N=i
           break()
       }
   }
   if (N==1) return(ml[[1]])
   re = list(domain =ml[[N]]$domain,
              image =ml[[1]]$image,
              defn  =list(type='morphism.composition',
                           defn=ml))
   class(re) = c('morphism', class(re))
   re
}

#' another handy interface for undetermined number of morphisms composition
#' @param ... morphisms
#' @return morphism as the composition result
morphism.composition.n = function(...) {
   morphism.composition.list(list(...))
}

#' extract the dimension of domain
#' @param m morphism
#' @return dimension number
degree.of.freedom = function(m) {
    # of course the Constant mapping can from any space, however, it doesn't cost degree of freedom
    # it's reasonable to set it as 0
    if (m$defn$type=='Constant') return(0)
    if (m$defn$type=='morphism.composition') {
    # composition morphism 's degree.of.freedom is defined by the last one, assuming all chain are compatible
    # if we directly read from domain, Constant will make trouble
        return(Recall(m$defn$defn[[ length(m$defn$defn) ]] ))
    }
    ambienceDomain = calculate.ambience.domain(m$domain)
    ambienceDomain[[3]]
}

#' return the dimension of the image
#' @param m morphism
#' @return dimension of the ambience space of the image
dimension.of.image = function(m) {
   calculate.ambience.domain( m$image )[[3]]
}

#' generate formal arguments, which not existing in the symbols
#' @param m morphism
#' @param fsymbol symbol used for arguments, e.g. THETA or x
#' @return list of indeterminators
generate.formal.arguments.morphism = function(m, fsymbol=quote(x)) {
    N = degree.of.freedom(m)
    if (N==0) return(NULL)
    lapply(1:N, function(ii) CONS('[', fsymbol, as.numeric(ii)))
}

#' for a given morphism export as one or more equations
#' just for applying on parameters
#' but not for currying
#' You can instantiate using THETAs as parameters, and then \code{explain.expression.as.morphism} to analyse
#' its structure
#' @param m morphism
#' @param parameter.pool e.g. \code{THETA[i]} which will be used for constructing the output mapping
#' @return equations
instantiate.morphism = function(m, parameter.pool=NULL) {
    if (m$defn$type=='Constant') {
        # no need for parameter.pool
        morphism.concrete.instantiate(m)
    }
    if (is.null(parameter.pool)) {
        parameter.pool = generate.formal.arguments.morphism(m, fsymbol=quote(THETA))
    } else if (!is.list(parameter.pool) && 
        (is.atomic(parameter.pool) || is.language(parameter.pool))) parameter.pool = list(parameter.pool)
    if (length(parameter.pool)!=degree.of.freedom(m)) {
        stop('number of parameter.pool and dimension of domain mismatch')
    }
    if (!is.null(re<-morphism.concrete.instantiate(m , parameter.pool))) {
        return(re)
    }
    if (!is.null(re<-morphism.groupeffect.instantiator(m, parameter.pool))) {
        return(re)
    }
    if (m$defn$type=='morphism.composition') {
    # for m1 %o% m2 %o% m3
    # first instantiate m3 and use is as parameter.pool to instantiate m2
    # and then m1
        N.m = length(m$defn$defn)
        if (N.m == 0) stop('morphism composition has nothing in it')
        if (N.m == 1) return(Recall(m$defn$defn[[1]], parameter.pool))
        last.mor = Recall(m$defn$defn[[N.m]], parameter.pool) 
        for(ind in (N.m-1):1) {
            last.mor = Recall(m$defn$defn[[ind]], last.mor)
        }
        return(last.mor)
    }
    stop(sprintf('Unknown type: %s', m$defn$type ))
}

#' FORMATTER
#' simple pretty printer for morphism
#' @param m morphism
#' @return string 
print.morphism0 = function(m, WW=60) {
    domain = paste(deparse(m$domain), collapse='')
    image = paste(deparse(m$image) , collapse='')
    if (m$defn$type=='Constant') {
        pars = NULL
        str.pars = '*'
    } else {
        if (degree.of.freedom(m)==1) pars = 'x'
        else pars = sprintf('x%s', 1:degree.of.freedom(m))
        str.pars = paste(pars,collapse=',')
        if (length(pars)>1) str.pars = sprintf('(%s)', str.pars)
    }
    e.s = instantiate.morphism(m, parameter.pool=lapply(pars, as.symbol))
    if(is.list(e.s)){
       str.img = sprintf('(%s)',paste(sapply(e.s, function(x) paste(deparse(simplify.pretty.01(x)),collapse='')) ,collapse=','))
    } else {
       str.img = paste(deparse(simplify.pretty.01(e.s)),collapse='')
    }
    str.img = gsub(' {3,}',' ', str.img)
    if (nchar(str.img)>WW) str.img = sprintf('%s ...',substring(str.img, 1, WW))
    stype = m$defn$type
    if (stype=='groupeffect') stype = sprintf('%s:%s', stype, m$defn$defn$grptype)
    s1 = sprintf('%s -----> %s  [%s]', domain, image, stype)
    s2 = sprintf('%s |--> %s', str.pars, str.img)
    sep0 = ifelse(nchar(s1)+nchar(s2)>80, '\n','')
    sprintf('%s  %s%s',s1,sep0,s2)
}

#' will deal with Composition special
#' printing additional information of the components morphism
#' affine mapping will also show the matrix
#' the structure will be shown more clearer
#' @param m morphism
#' @return NULL
print.morphism = function(m) {
   indent.n = function(ss,n) {
       if (n==0) return(ss)
       indent = paste(rep(' ',n),collapse='')
       ss = gsub('^', indent, ss)
       gsub('\n', sprintf('\n%s',indent), ss)
   }
   comp.types = c('morphism.composition','Sum.of.Morphisms','Product.of.Morphisms')
   summary.of.type = c('Composition','Sum','Product')
   walker = function(m1,id0=0) {
       ss = print.morphism0(m1)
       ss = indent.n(ss, id0)
       id1 = id0 + 4
       if (!is.na(tp1 <- match(m1$defn$type,comp.types))) {
           N = length(m$defn$defn) 
           if (N>1) {
               ss = sprintf('%s\n%s',ss,
                        indent.n(sprintf('%s of %s morphisms', summary.of.type[tp1] ,
                            N),id0+2))
               for(i in 1:N) {
                   sub.mor = Recall(m$defn$defn[[i]], id1)
                   # add index
                   sub.mor = sub('^(\\s+)',sprintf('\\1%s) ',i), sub.mor)
                   ss = sprintf('%s\n%s', ss, sub.mor)
               }
           }
       }
       if (m1$defn$type=='AffineMapping') {
           s2 = indent.n(print.affine.mapping0(m1$defn$defn),id0+10)
           ss = sprintf('%s\n%s', ss, s2)
       }
       if (m1$defn$type=='affine.mapping.factor.effect') {
           s2 = indent.n(print.affine.mapping.factor.effect0(m1$defn$defn),id0+10)
           ss = sprintf('%s\n%s', ss, s2)
       }
       ss
   }
   re=walker(m)
   cat(re,'\n')
}
