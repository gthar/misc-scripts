#!/usr/bin/env Rscript

###############################################################################
# Input parameters

# input files
xfile <- "newCG_twist.dat"
yfile <- "newCG_tilt.dat"
xefile <- "CGe_twist.dat"
yefile <- "CGe_tilt.dat"

# output file
outfile <- "GC_tilt_roll.tiff"

# 2D area range
xlim <- c(0, 55)
ylim <- c(-25, 25)
zlim <- c(0, 0.004)

# axis labels
xlab <- "tilt"
ylab <- "roll"

# number of bins in the contour
nbins <- 4

# placement of the segment legend relative to xlim and ylim
legend.xprop <- 0.15
legend.yprop <- 0.2

# offset of the number labels to the position of the line, relative to ylim
lab.offset.prop <- 0.015

###############################################################################
# imports

library(ggplot2)
library(MASS)
library(plyr)
library(reshape2)

###############################################################################
# Some function definitions


getLimProp <- function (lim, prop, at.start=TRUE)
{   # Given a limit range and a proportion, return a sub-limit at the start
    # or at the end of the original limit
    segment.len <- (lim[2] - lim[1]) * prop
    if (at.start) {
        c(lim[1], lim[1] + segment.len)
    } else {
        c(lim[2] - segment.len, lim[2])
    }
}

getSegmentLegendDf <- function (xlim, ylim, xprop, yprop, n,
                                xplace="right", yplace="bottom")
{   # Create a data.frame with the data for segments to be
    # represented as an ad hoc legend for the contour
    if (xplace == "right") {
        subxlim <- getLimProp(xlim, xprop, at.start=FALSE)
    } else if (xplace == "left") {
        subxlim <- getLimProp(xlim, xprop, at.start=TRUE)
    }
    if (yplace == "bottom") {
        subylim <- getLimProp(ylim, xprop, at.start=TRUE)
    } else if (yplace == "top") {
        subylim <- getLimProp(ylim, xprop, at.start=FALSE)
    }

    total.len <- subylim[2] - subylim[1]
    sect.len <- total.len/n

    ypos <- subylim[1] + 0:(n-1)*sect.len

    data.frame(x1=subxlim[1],
               x2=subxlim[2],
               y=ypos)
}

###############################################################################
# Read the input files

# read the values in the input files
vals <- lapply(list(s=list(x=xfile,
                           y=yfile),
                    e=list(x=xefile,
                           y=yefile)),
               lapply,
               scan)

# arrange the input files values into data.frames
dfs <- lapply(vals,
              with,
              data.frame(x, y))

###############################################################################
# calculate 2D density over a grid to colour points and to draw contours

dens <- with(dfs$s, kde2d(x, y, n=100, lims=c(xlim, ylim)))

m <- dens$z

colnames(m) <- dens$y
rownames(m) <- dens$x

dens.df <- melt(m)
names(dens.df) <- c("x", "y", "z")

###############################################################################

compose <- function (...) {
    comp2 <- function(f, g) {
        force(f)
        force(g)
        function (...) f(g(...))
    }
    Reduce(comp2, list(...))
}

getClosests <- function (xs, is)
    sapply(xs, compose(which.min, abs, `-`), is)

is <- getClosests(dfs$s$x, dens$x)
js <- getClosests(dfs$s$y, dens$y)

dfs$s$density <- mapply(function (i, j) dens$z[i, j], is, js)

###############################################################################
# let's find where to draw labels with numbers on the plot

###############################################################################

p <- ggplot(dens.df, aes(x, y, z=z)) + stat_contour(bins=nbins)
val.df <- ggplot_build(p)$data[[1]]

segments <- getSegmentLegendDf(xlim,
                               ylim,
                               legend.xprop,
                               legend.yprop,
                               nbins,
                               xplace="right",
                               yplace="bottom")

segments$txt <- round(sort(unique(val.df$level),
                           decreasing=FALSE),
                      4)
segments$labx <- with(segments, (x1+x2) / 2)
segments$laby <- segments$y + (ylim[2] - ylim[1]) * lab.offset.prop

###############################################################################
# Make the plot

myplot <- ggplot() +
    geom_point(data=dfs$s,
               mapping=aes(x=x, y=y, color=density)) +
    stat_contour(data=dens.df,
                 mapping=aes(x=x, y=y, z=z),
                 colour="#FFFFFF",
                 alpha=0.5,
                 bins=nbins) +
    geom_point(data=dfs$e, aes(x, y), size=2, alpha=0.7) +
    scale_colour_gradientn(colours=c("#0000FF",
                                     "#00FFFF",
                                     "#FFFF00",
                                     "#FF0000"),
                           limits=zlim) +
    geom_segment(data=segments,
                 (aes(x=x1, y=y, xend=x2, yend=y)),
                 alpha=0.5) +
    geom_text(data=segments,
              aes(x=labx, y=laby, label=txt)) +
    xlab(xlab) +
    ylab(ylab) +
    coord_cartesian(xlim=xlim, ylim=ylim) +
    theme(axis.text   = element_text(size=20),
          axis.title  = element_text(size=20),
          legend.text = element_text(size=20)) +
    theme_bw()

###############################################################################
# And save it

myplot
#ggsave(filename=outfile, plot=myplot, height=5, width=6)

###############################################################################
