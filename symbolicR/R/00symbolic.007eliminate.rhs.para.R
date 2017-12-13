#' eliminate duplicated appearence of NONMEM parameters
#'
#' multiple appearence of THETA[i] or ETA[i] or EPS[i] as rhs in equations may cause problem
#' because of our assumption any appearence of THETA[] ETA[] EPS[] will create new estimable.unit 
#' or random observations
#' one simple solution is, if find multiple rhs THETA[1] then introduce a symbol parameter and replace all
#' appearence of THETA[1] by the introduced symbol
#'
#' @param eqns list of equations
#' @return equivalent equations
eliminate.duplicated.appearence.of.rhs.parameters = function(eqns) {
    first.dup.ind = function(gg, pat=quote(THETA)){
        for(i in seq_along(gg$vertex.table)) {
            if (symbolic.match(CONS('[',pat,quote('?a'(parind))), gg$vertex.table[[i]]) ) {
                if (sum(sapply(gg$eqns.lhs.rhs, function(x) i %in% x$rhs))>1) return(i)
            }
        }
        NULL
    }
    newvar = numeric(3)
    names(newvar) = rev(c('THETA','ETA','EPS'))
    intro.code = list(EPS=list(),ETA=list(),THETA=list())
    for(pat in names(newvar)) {
        while(TRUE) {
            gg = create.abstractGraph.from.equations(eqns)
            ind = first.dup.ind(gg, as.symbol(pat))
            if (is.null(ind)) break()
            # introduce new symbols
            newvar[ pat ] = newvar[pat]+1
            # symbols
            oldsym = gg$vertex.table[[ ind ]]
            newsym = as.symbol(sprintf('%sHAT%02d', pat, newvar[pat]))
            intro.code[[pat]][[ length(intro.code[[pat]])+1 ]] = CONS('=', newsym, oldsym)
            for(i in seq_along(eqns)) {
                if (ind %in% gg$eqns.lhs.rhs[[i]]$rhs) {
                    eqns[[i]] = symbolic.gsub(oldsym,newsym,eqns[[i]])
                }
            }
        }
    }
    list(new.eqns = eqns, intro.eqns=intro.code, rep.var = newvar)
}
