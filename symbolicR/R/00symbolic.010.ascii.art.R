#' display an expression as a tree
#' 
#' Only for debug and learning purpose
#'
#' @param e0 expression
#' @param display \code{TRUE} for print the result out
#' @param do.not.expand We may avoid expanding the call starting by '?'. (If you do not want to expand '[', this can be set here)
#' @return string for the tree
symbolic.draw.tree = function(e0, display=TRUE, do.not.expand = c('?','?a','?s','?c','?()','?[]','??','?na')) {
    charToBlock = function(ce) {
        list(val = unlist(strsplit(ce,'')),
                dim = c(1, nchar(ce)),
                anchor = rep(max(round(nchar(ce)/2),1),2))
    }
    embeding = function(b, val, new.pos, val.size) {
        d1 = b$dim[1]
        d2 = b$dim[2]
        for(i in 1:d1) {
            for(j in 1:d2) {
                pos = (i-1)*d2+j
                newpos = (i-2+new.pos[1])*val.size[2] + (j-1+new.pos[2])
                val[newpos] = b$val[pos]
            }
        }
        val
    }
    # +----> y
    # |
    # |
    # V
    # x
    varrow = function(l=1) {
        list(val=rep('|',l), dim=c(l,1), anchor=c(1,1))
    }
    plus.sym = function(){
        list(val='+', dim=c(1,1), anchor=c(1,1))
    }
    hlink = function(l){
        list(val=rep('-',l), dim=c(1,l), anchor=c(ceiling(l/2),ceiling(l/2)))
    }
    vbind = function(b1, b2) {
        new.dim.1 = b1$dim[1] + b2$dim[1] + 1
        new.dim.2 = max( b1$anchor[2], b2$anchor[1]) + 
            max( b1$dim[2] - b1$anchor[2], b2$dim[2] - b2$anchor[1] )
        new.dim = c(new.dim.1,new.dim.2)
        new.val = rep(' ', new.dim.1 * new.dim.2)
        new.anchor = c(0,0)
        b1.x = 1
        b2.x = b1$dim[1] + 2
        delta.y = b1$anchor[2] - b2$anchor[1]
        if (delta.y >= 0) {
            b1.y = 1
            b2.y = 1 + delta.y
            new.anchor = c(b1$anchor[1], b2$anchor[2]+delta.y)
        } else {
            b1.y = 1 - delta.y
            b2.y = 1
            new.anchor = c(b1$anchor[1] - delta.y, b2$anchor[2])
        }
        new.val = embeding(b1, new.val, c(b1.x, b1.y), new.dim)
        new.val = embeding(b2, new.val, c(b2.x, b2.y), new.dim)
        new.val = embeding(varrow(1), new.val, 
                        c(b1.x+b1$dim[1], b1.y + b1$anchor[2] - 1), 
                        new.dim)
        list(val=new.val, dim=new.dim, anchor=new.anchor)
    }
    hbind = function(bl) {
        n = length(bl)
        new.dim.1 = max(sapply(bl, function(x) x$dim[1])) + 2
        new.dim.2 = sum(sapply(bl, function(x) x$dim[2])) + n - 1
        new.val = rep(' ', new.dim.1 * new.dim.2)
        new.dim = c(new.dim.1,new.dim.2)
        cy = 1
        last.a = a0 = NA
        for(i in 1:n) {
            new.val = embeding(bl[[i]],new.val,c(3, cy), new.dim)
            a = cy+bl[[i]]$anchor[1]-1
            # |
            new.val = embeding(varrow(1),new.val, c(2,a), new.dim)
            new.val = embeding(plus.sym(),new.val, c(1,a), new.dim)
            if (i>1) {
            # ---
                new.val = embeding(hlink(a-last.a-1),new.val, c(1,last.a+1), new.dim)
            } else {
                a0 = a
            }
            last.a = a
            cy = cy + bl[[i]]$dim[2] + 1
        }
        mid = a0 + round((a-a0)/2)
        # put + at middle
        new.val = embeding(plus.sym(),new.val, c(1,mid), new.dim)
        list(val=new.val, dim=new.dim, anchor=c(mid,NA))
    }
    pretty.single = function(e) {
        #
        if (is.symbol(e)) {
            ce = sprintf('`%s`',as.character(e))
        } else if (is.character(e)) {
            ce = dQuote(as.character(e))
        } else if (is.numeric(e)) {
            ce = as.character(e)
        }
        ce
    }
    walker = function(e) {
        if (is.atomic(e) || is.symbol(e)) {
            return(charToBlock(pretty.single(e)))
        }
        # assert it's a call
        stopifnot(is.call(e))
        ##
        root = e[[1]]
        if (! (is.symbol(root) || is.atomic(root)) ) {
            stop('The function is too complex - itself is another call.')
        }
        if (as.character(root) %in% do.not.expand) {
            ce = paste(deparse(e),collapse='')
            ce = gsub('\\s{2,}',' ',ce)
            return(charToBlock(ce))
        }
        root = charToBlock(pretty.single(root))
        #
        arity = length(e) - 1
        if (arity == 0) {
        # no argument
            return(root)
        } else if (arity == 1) {
        # one argument
            arg1 = Recall(e[[2]])
            return(vbind(root,arg1))
        } else {
        # more than 2 arguments
            args = vector(arity, mode='list')
            for(i in 1:arity) {
                args[[i]] = Recall(e[[ i + 1 ]])
            }
            return(vbind(root,hbind(args)))
        }
    }
    ascii = walker(e0)
    print.block = function(b){
        d1 = b$dim[1]
        d2 = b$dim[2]
        lines = character(d1)
        for(i in 1:d1) {
            lines[i] = paste(b$val[(i-1)*d2 + (1:d2)],collapse='')
        }
        txt = paste(lines,collapse='\n')
        txt
    }
    re = print.block(ascii)
    if (display) cat(re,'\n')
    invisible(re)
}
