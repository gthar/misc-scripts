#!/usr/bin/env Rscript

###############################################################################
# Input parameters

# input files
xfile <- "newAT_twist.dat"
yfile <- "newAT_tilt.dat"
xefile <- "ATe_twist.dat"
yefile <- "ATe_tilt.dat"

# output file
outfile <- "2d_densplot.pdf"

# 2D area range
xlim <- c(0, 60)
ylim <- c(-20, 20)

# axis labels
xlab <- "twist"
ylab <- "tilt"

###############################################################################
# imports

library(ggplot2)
library(MASS)
library(plyr)
library(reshape2)

###############################################################################
# A couple of FP definitions

compose <- function (...)
{   # Function composition
    comp2 <- function (f, g) {
        force(f)
        force(g)
        function (x) f(g(x))
    }
    Reduce(comp2, list(...))
}

partial <- function (f, ...)
{   # Partial application
    capture <- list(...)
    function (...) do.call(f, c(list(...), capture))
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

dens <- with(dfs$s, kde2d(x, y, n=50, lims=c(xlim, ylim)))

# create a new data frame of that 2d density grid
# (needs checking that I haven't stuffed up the order here of z?)
gr <- with(dens,
           data.frame(expand.grid(x, y),
                      as.vector(z)))
names(gr) <- c("xgr", "ygr", "zgr")

mod <- loess(zgr ~ xgr*ygr, data=gr)

dfs$s$density <- predict(mod,
                         newdata=data.frame(xgr=dfs$s$x,
                                            ygr=dfs$s$y))

###############################################################################
# create a melted data.frame that represents a matrix of the densities

m <- dens$z
colnames(m) <- dens$y
rownames(m) <- dens$x

dens.df <- melt(m)
names(dens.df) <- c("x", "y", "z")

###############################################################################
# let's find where to draw labels with numbers on the plot

p <- ggplot(dens.df, aes(x, y, z=z)) + stat_contour()
val.df <- ggplot_build(p)$data[[1]]

# from `val.df`, we need an `x` present in all `level` and where each label
# appears only once

getLevels <- compose(sort,
                     unique,
                     partial(`[[`,
                             "level"))

# separate them by x value
by.x <- dlply(val.df, .(x), identity)
# get the ones that contain all levels
sel <- sapply(lapply(by.x,
                     getLevels),
              identical,
              getLevels(val.df))
# we keep the first one that contains all levels
txt.df <- by.x[sel][[1]]

###############################################################################
# Make the plot

myplot <- ggplot() +
    geom_point(data=dfs$s,
               mapping=aes(x=x, y=y, color=density)) +
    stat_contour(data=dens.df,
                 mapping=aes(x=x, y=y, z=z),
                 colour="#000000",
                 alpha=0.5) +
    geom_text(data=txt.df,
              aes(x=x, y=y, z=NULL, label=level),
              size=2.5,
              alpha=1) +
    geom_point(data=dfs$e, aes(x, y), size=2, alpha=0.5) +
    scale_colour_gradientn(colours=c("#0000FF",
                                     "#00FFFF",
                                     "#FFFF00",
                                     "#FF0000")) +
    xlab(xlab) +
    ylab(ylab) +
    coord_cartesian(xlim=xlim, ylim=ylim) +
    theme_bw()
myplot

###############################################################################
# And save it

ggsave(filename=outfile, plot=myplot)

###############################################################################
