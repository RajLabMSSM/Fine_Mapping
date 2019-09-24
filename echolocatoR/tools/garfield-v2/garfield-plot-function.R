# The function in this file is based on the radial.plot function
# from the R package plotrix. The original code is released under
# the GPL v2 license, which has been modified here and all changes
# to the original code are released under the GPL v2 license with
# copyright notice for these included below.


# Modified code license:

# GARFIELD - GWAS analysis of regulatory or functional information enrichment with LD correction.
# Copyright (C) 2014 Wellcome Trust Sanger Institute / EMBL - European Bioinformatics Institute
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#print ("reading circus_plot_script.R")

shadowtext <- function(x, y=NULL, labels, col='white', bg='black',bg2='white',
	theta= seq(0, 2*pi, length.out=32), r=0.02, ... ) {
	

	xy <- xy.coords(x,y)
	xo <- r*strwidth('A')
	yo <- r*strheight('A')
	for (i in theta) {
		text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, 
		  labels, col=bg, ... )
	}
	#xy <- xy.coords(x,y)
	#xo <- r/10*strwidth('A')
	#yo <- r/10*strheight('A')
	#for (i in theta) {
#		text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, 
#		  labels, col=bg2, ... )
#	}
	text(xy$x, xy$y, labels, col=col, ... )
}


garfield.plot<-function (lengths, radial.pos = NULL, labels = NA, breaks= NA, label.pos = NULL, 
    radlab = FALSE, start = 0, clockwise = FALSE, rp.type = "r", 
    label.prop = 1.05, main = "", xlab = "", ylab = "", line.col = par("fg"), 
    lty = par("lty"), lwd = par("lwd"), mar = c(2, 2, 3, 2), 
    show.grid = TRUE, show.grid.labels = 4, show.radial.grid = TRUE, 
    rad.col = "gray", grid.col = "gray", grid.bg = "transparent", 
    grid.left = FALSE, grid.unit = NULL, point.symbols = 1, point.col = par("fg"), 
    show.centroid = FALSE, radial.lim = NULL, radial.labels = NULL, 
    poly.col = NA, add = FALSE, ann.col=1,ann.pch=15,ann.cols=1,compact=TRUE,...) 
{	
    lengths.ann = rep(max(lengths,na.rm=T),ncol(lengths))
    radial.lim=c(0,max(lengths,na.rm=T))
    
    #if (is.null(radial.lim)) 
    #    radial.lim <- range(lengths)
    length.dim <- dim(lengths)
    if (is.null(length.dim)) {
        npoints <- length(lengths)
        nsets <- 1
        lengths <- matrix(lengths, nrow = 1)
    }
    else {
        npoints <- length.dim[2]
        nsets <- length.dim[1]
        lengths <- as.matrix(lengths)
    }
    lengths <- lengths - radial.lim[1]
    lengths[lengths < 0] <- NA
    if (is.null(radial.pos[1])) 
        radial.pos <- seq(0, pi * (2 - 2 * (rp.type != "l")/npoints), 
            length.out = npoints)
    radial.pos.dim <- dim(radial.pos)
    if (is.null(radial.pos.dim)) 
        radial.pos <- matrix(rep(radial.pos, nsets), nrow = nsets, 
            byrow = TRUE)
    else radial.pos <- as.matrix(radial.pos)
    if (rp.type == "l") {
        clockwise <- TRUE
        start <- pi/2
    }
    if (clockwise) 
        radial.pos <- -radial.pos
    if (start) {
        radial.pos <- radial.pos + start
    }
    if (show.grid) {
        if (length(radial.lim) < 3) 
            grid.pos <- pretty(radial.lim)
        else grid.pos <- radial.lim
        if (grid.pos[1] < radial.lim[1]) 
            grid.pos <- grid.pos[-1]
        maxlength <- max(grid.pos - radial.lim[1])
        angles <- seq(0, 1.96 * pi, by = 0.04 * pi)
    }
    else {
        grid.pos <- NA
        maxlength <- diff(radial.lim)
    }
    oldpar <- par("xpd", "mar", "pty")
    if (!add) {
        par(mar = mar, pty = "s")
        plot(c(-maxlength, maxlength), c(-maxlength, maxlength), 
            type = "n", axes = FALSE, main = main, xlab = xlab, 
            ylab = ylab)
        if (show.grid) {
            for (i in seq(length(grid.pos), 1, by = -1)) {
                xpos <- cos(angles) * (grid.pos[i] - radial.lim[1])
                ypos <- sin(angles) * (grid.pos[i] - radial.lim[1])
                polygon(xpos, ypos, border = grid.col, col = grid.bg)
            }
        }
    }
    par(xpd = TRUE)
    if (length(line.col) < nsets) 
        line.col <- 1:nsets
    if (length(rp.type) < nsets) 
        rp.type <- rep(rp.type, length.out = nsets)
    if (length(point.symbols) < nsets) 
        point.symbols <- rep(point.symbols, length.out = nsets)
    if (length(point.col) < nsets) 
        point.col <- rep(point.col, length.out = nsets)
    if (length(poly.col) < nsets) 
        poly.col <- rep(poly.col, length.out = nsets)
    if (length(lty) < nsets) 
        lty <- rep(lty, length.out = nsets)
    if (length(lwd) < nsets) 
        lwd <- rep(lwd, length.out = nsets)

    for (i in 1:nsets) {
        if (nsets > 1) {
            linecol <- line.col[i]
            polycol <- poly.col[i]
            pointcol <- point.col[i]
            pointsymbols <- point.symbols[i]
            ltype <- lty[i]
            lwidth <- lwd[i]
        }
        else {
            linecol <- line.col
            polycol <- poly.col
            pointcol <- point.col
            pointsymbols <- point.symbols
            ltype <- lty
            lwidth <- lwd
        }
        rptype <- unlist(strsplit(rp.type[i], ""))
        if (match("s", rptype, 0)) {
            if (is.null(pointsymbols)) 
                pointsymbols <- i
            if (is.null(pointcol)) 
                pointcol <- i
        }
        xpos <- cos(radial.pos[i, ]) * lengths[i, ]
        ypos <- sin(radial.pos[i, ]) * lengths[i, ]
        if (match("r", rptype, 0)) 
            segments(0, 0, xpos, ypos, col = linecol, lty = ltype, 
                lwd = lwidth, ...)
        if (match("p", rptype, 0)) 
            polygon(xpos, ypos, border = linecol, col = polycol, 
                lty = ltype, lwd = lwidth, ...)
        if (match("s", rptype, 0)) 
            points(xpos, ypos, pch = pointsymbols, col = pointcol, 
                ...)
        if (match("l", rptype, 0)) 
            lines(xpos, ypos, lty = ltype, lwd = lwidth, col = linecol, 
                ...)
        if (show.centroid) 
            if (match("p", rptype, 0)) {
                nvertices <- length(xpos)
                polygonarea <- xpos[nvertices] * ypos[1] - xpos[1] * 
                  ypos[nvertices]
                for (vertex in 1:(nvertices - 1)) polygonarea <- polygonarea + 
                  xpos[vertex] * ypos[vertex + 1] - xpos[vertex + 
                  1] * ypos[vertex]
                polygonarea <- polygonarea/2
                centroidx <- (xpos[nvertices] + xpos[1]) * (xpos[nvertices] * 
                  ypos[1] - xpos[1] * ypos[nvertices])
                centroidy <- (ypos[nvertices] + ypos[1]) * (xpos[nvertices] * 
                  ypos[1] - xpos[1] * ypos[nvertices])
                for (vertex in 1:(nvertices - 1)) {
                  centroidx <- centroidx + (xpos[vertex] + xpos[vertex + 
                    1]) * (xpos[vertex] * ypos[vertex + 1] - 
                    xpos[vertex + 1] * ypos[vertex])
                  centroidy <- centroidy + (ypos[vertex] + ypos[vertex + 
                    1]) * (xpos[vertex] * ypos[vertex + 1] - 
                    xpos[vertex + 1] * ypos[vertex])
                }
                points(centroidx/(6 * polygonarea), centroidy/(6 * 
                  polygonarea), col = point.col[i], pch = point.symbols[i], 
                  cex = 2, ...)
            }
            else points(mean(xpos), mean(ypos), col = pointcol, 
                pch = pointsymbols, cex = 2, ...)

    }
    xpos.ann <- cos(radial.pos[1, ]) * maxlength*1.02
	 #ypos.ann <- sin(radial.pos[1, ]) * lengths.ann*1.09    
    ypos.ann <- sin(radial.pos[1, ]) * maxlength*1.02

    points(xpos.ann, ypos.ann, pch = ann.pch, col = ann.col,...)
    for (ij in 1:nrow(ann.cols)){
        xpos.ann2 <- cos(radial.pos[1, ]) * (maxlength*(0.992-0.015*(ij-1)))#1.015)
        ypos.ann2 <- sin(radial.pos[1, ]) * (maxlength*(0.992-0.015*(ij-1)))#1.015)
        points(xpos.ann2, ypos.ann2, pch = 19, col = ann.cols[ij,], cex=0.5,...) 
    }
    if (!add) {
        if (is.na(labels[1])) {
            label.pos <- seq(0, 1.8 * pi, length = 9)
            labels <- as.character(round(label.pos, 2))
 	    lablen <- length(labels)
	    breaks.pos <- seq(0, pi * (2 - 2/length(breaks)), length.out = length(breaks))
        }
        if (is.null(label.pos[1])) {
            lablen <- length(labels)
            label.pos <- seq(0, pi * (2 - 2/lablen), length.out = lablen)
	    breaks.pos <- seq(0, pi * (2 - 2/length(breaks)), length.out = length(breaks))
        }
        if (clockwise) 
            label.pos <- -label.pos
        if (start) 
            label.pos <- label.pos + start
        xpos <- cos(breaks.pos-pi/lablen) * maxlength
        ypos <- sin(breaks.pos-pi/lablen) * maxlength
	
   	if (show.radial.grid) 
    		segments(0, 0, xpos[which(!duplicated(breaks))], ypos[which(!duplicated(breaks))], col = "gray",lty=4)        
		  xpos <- cos(label.pos) * maxlength * label.prop
        ypos <- sin(label.pos) * maxlength * label.prop
        label.adj <- round(abs(1 - cos(label.pos))/2-10^{-10})
        if (radlab) {
	labn = names(table(labels))
	nlab = as.numeric(table(labels))/sum(as.numeric(table(labels)))*length(as.numeric(table(labels)))
	labelsn = nlab[match(labels,labn)]
            for (label in 1:length(labels)) {
		if (!is.na(labels[label])){
			flr=ceiling(median(which(labels == labels[label])))
		} else {
			flr=ceiling(median(which(is.na(labels))))
		}
		if ((flr==label & compact==TRUE) | compact==FALSE){
                labelsrt <- (180 * label.pos[label]/pi) + 180 * 
                  (label.pos[label] > pi/2 && label.pos[label] < 
                    3 * pi/2)
                ##text(xpos[label], ypos[label], labels[label],cex = 0.45, srt = labelsrt, adj = label.adj[label])
                #text(xpos[label], ypos[label], labels[label],cex = (0.35+0.3*labelsn[label]), srt = labelsrt, adj = label.adj[label]+lengths.ann[1]*1.07/1500,col="gray30")
		#text(xpos[label], ypos[label], labels[label],cex = (0.35+0.3*labelsn[label]), srt = labelsrt, adj = label.adj[label],col=ann.col[label])
		shadowtext(xpos[label], ypos[label], labels[label],cex = (0.45+0.3*labelsn[label]), srt = labelsrt, adj = label.adj[label],bg="white",col=1)
#		shadowtext(xpos[label], ypos[label], labels[label],cex = (0.25+0.3*labelsn[label]), srt = labelsrt, adj = label.adj[label],bg="black",col=ann.col[label])
		
		}
            }
        }
        else {
            for (label in 1:length(labels)) {
                text(xpos[label], ypos[label], labels[label], 
                  cex = par("cex.axis"), adj = label.adj[label])
            }
        }
        if (show.grid.labels) {
            if (show.grid.labels%%2) {
                ypos <- grid.pos - radial.lim[1]
                xpos <- rep(0, length(grid.pos))
                if (show.grid.labels == 1) 
                  ypos <- -ypos
            }
            else {
                xpos <- grid.pos - radial.lim[1]
                ypos <- rep(0, length(grid.pos))
                if (show.grid.labels == 2) 
                  xpos <- -xpos
            }
            if (is.null(radial.labels)) 
                radial.labels = as.character(grid.pos)
            if (!is.null(grid.unit)) 
                radial.labels[length(grid.pos)] <- paste(radial.labels[length(grid.pos)], 
                  grid.unit)
           text(xpos+0.003*max(xpos), ypos+0.003*max(ypos), radial.labels, cex = par("cex.lab"),col="white")
	   text(xpos, ypos, radial.labels, cex = par("cex.lab"))
	
        }
    }
    invisible(oldpar)
}

