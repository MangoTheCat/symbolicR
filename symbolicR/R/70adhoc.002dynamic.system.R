#' description as kinetic.system
#' @param ncompartment number of compartments
#' @param defn definition of dynamics, can be kinetic system or even general ODE
#' @return skeleton
dynamic.system.new = function(ncompartment=3,defn){
    re = list(
            ncompartment=ncompartment,
            defn=defn
    )
    class(re) = union('dynamic.system', class(re))
    re
}

#' S3
print.dynamic.system = function(da) {
    print(da$defn)
}

#' S3
print.kinetic.system.predpp = function(da){
    if (da$type=='PREDPP') {
        s = sprintf('NONMEM PREDPP, %s', da$subs)
        if (!is.null(da$trans)) {
            s = sprintf('%s using %s',s, da$trans)
        }
        cat(s, '\n')
    }else{
        print(da)
    }
}

#' create a predpp known kinetic system
#' @param s description from NONMEM model
#' @return kinetic system
kinetic.system.predpp.new = function(s) {
    subs = s[1]
    if (length(s)>1) trans = s[2]
    else trans =''
    trans = sub('^TRAN.*(\\d+)', 'TRANS\\1', trans)
    if (is.null(NONMEM.SUBROUTINE.DATABASE[[ subs ]])) {
        stop(sprintf('Can not find PREDPP subroutine %s ', subs))
    }
    ncomp = length(NONMEM.SUBROUTINE.DATABASE[[ subs ]]$compartments)
    defn = list(type='PREDPP', subs=subs, trans=trans)
    class(defn) = 'kinetic.system.predpp'
    re = dynamic.system.new(
            ncompartment=ncomp,
            defn = defn )
    re
}

#' convert description in kinetic to NONMEM PREDPP subroutine
#' @param d kinetic description
#' @return string used in NONMEM control file SUBR block
as.nonmem.subroutine.dynamics = function(d){
    stopifnot(is(d$defn, 'kinetic.system.predpp'))
    sprintf('%s %s',d$defn$subs,d$defn$trans)
}

#' convert description in kinetic to NONMEM PREDPP subroutine, as characters
#' @param d kinetic description
#' @return strings of subs and trans
get.subs.trans.dynamics = function(d){
    stopifnot(is(d$defn, 'kinetic.system.predpp'))
    c(d$defn$subs, d$defn$trans)
}
