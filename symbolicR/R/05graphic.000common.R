# Graphic relation functions                  #
###############################################

#' topo.sort
#' Topological sort according to dependency
#'
#' Assume \code{G} is directed graph as following: \cr
#'  \code{G=list(V=V,E=E)}
#' Do a topological sort of vertex of \code{G}
#'
#' @param G graph
#' @param rev logical \code{TRUE} then reverse direction of all edges
#' @return a list has two names \code{L} and \code{stage.number} \cr 
#'      where \code{L} is the sorted vertex and \code{N} are the corresponding stage number for the vertex in \code{L}
#' @author jjxie
topo.sort=function(G,rev=F){
    from = 1
    to   = 2
    if (rev) {
        from=2
        to = 1
    }
    # S = no incoming
    stage.number = integer(0)
    S = setdiff(G$vertex, sapply(G$edges, function(x)x[to]))
    stage.number[ as.character(S) ] = 1
    L = character(0)
    E = G$edges
    while( length(S) > 0 ){
        n = S[1]
        S = S[-1]
        L = c(L,n)
        E1 = integer(0)
        N.E = length(E)
        if (N.E==0){
            L=c(L,S)
            break
        }
        for(i in 1:N.E){
            k = E[[i]]
            if (k[from] != n) {
                E1 = c(E1,i)
                next
            }
            m = k[to]
            if (any(sapply(1:N.E, function(x) E[[x]][from]!=n && E[[x]][to]==m))) next()
            stage.number[ as.character(m) ] = stage.number[ as.character(n) ] + 1
            S = c(S,m)
        }
        E = E[E1]
    }
    #
    if (length(E)>0) stop('Loop found!')
    list(L=L, N=stage.number)
}

#' restriction.morphism.of.graph.at.v
#'
#' sub graph of \code{G} restricted by vertex subset \code{v}, \cr
#' only the morphism is essential
#'
#' @param G : \code{G} is a directed graph, which is represented as list of vertex and edges; \code{edges} are list of vertex pairs <from, to>
#' @param v1 : to which vertex the Graph is restricted at
#' @return the indices of the edges appearing in the graph
#' @author jjxie
restriction.morphism.of.graph.at.v = function(G,v1){
    e = G$edges
    if(length(e)<1) return(integer(0))
    ind = sapply(1:length(e), function(x) all(e[[x]] %in% v1))
    which(ind)
}

#' neighbourhood.of.v
#'
#' For a vertex set \code{v1}, find the neighbourhood of \code{v1} in the directed Graph \code{G} \cr
#' return \code{v.income} and \code{v.outcome}, i.e. the incoming vertex to \code{v1} and output vertex from \code{v1}
#'
#' @param G : \code{G} is a directed graph, which is represented as list of vertex and edges; edges are list of vertex pairs <from, to>
#' @param v1 : vertex set
#' @return list of incoming and outcoming list
#' @author jjxie
neighbourhood.of.v = function(G, v1){
    from = 1
    to   = 2
    e = G$edges
    v.outcome <- v.income <- character(0)
    for(ei in e){
        if (ei[from] %in% v1 && !(ei[to] %in% v1)) {
            v.outcome = union(v.outcome, ei[to])
            next
        }
        if (ei[to] %in% v1 && !(ei[from] %in% v1)) {
            v.income = union(v.income, ei[from])
            next
        }
    }
    list(v.income=v.income, v.outcome=v.outcome)
}

#' connected.component
#' 
#' @title get the connected component
#' @param G graph
#' @param v vertices set from which looking
#' @param direction 0 denotes incoming and 1 denotes outcoming, and c(0,1) for both directions
#' @return set of vertex
#' @author jjxie
connected.component = function(G, v, direction=0){
    if (length(G$edges)==0) return(v)
    vs = v
    while(T){
        ind1 = which(sapply(G$edges, function(x) any(x[2-direction] %in% vs)))
        if (length(ind1)>0){
            v1 = sapply(G$edges[ind1], function(x) x[1+direction])
            v2 = union(v1,vs)
            if (setequal(v2,vs)) break
        } else {
        # isolated vs (usually just v itselt)
            return(vs)
        }
        vs = v2
    }
    vs
}

#' divide graph into connected.components
#' this will go wrong if more vertex than edges
#' you have to put self loop if you want those vertex being component
#' @param G graph
#' @return decomposition of vertex
divide.graph = function(G){
    if (length(G$vertex)==0) return(NULL)
    L=list()
    nL=0
    while(T) {
        v0 = connected.component(G,G$vertex[1],c(0,1))
        L[[ nL<-nL+1 ]] = v0
        v1 = setdiff(G$vertex, v0)
        if (length(v1)==0) break()
        E = restriction.morphism.of.graph.at.v(G, v1)
        if (length(E)==0) {
            L = c(L, as.list(v1))
            break()
        }
        G = list(vertex=v1,edges=G$edges[E])
    }
    return(L)
}
