#The following script outlines the code for main analysis in the manuscript
#"Traditional plant functional groups explain variation in economic 
# but not size-related traits across the tundra biome"

#Script written by Haydn J.D. Thomas
#August 2018

#Sleeping willow is feeling unique

# --- #

#1. Load packages and functions####
library(data.table)
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggbiplot) #Note: ggbiplot can be downloaded from https://github.com/vqv/ggbiplot
library(gridExtra)
library(mclust)
library(vegan)
library(ggdendro)
library(maps)
library(mapdata)
library(mapproj)
library(ggpubr)
se <- function(x) sd(x)/sqrt(length(x))
`%notin%` <- function(x,y) !(x %in% y)
normalise<- function(x) (x-min(x))/(max(x)-min(x))
mround <- function(x,base){ 
  base*round(x/base) 
} 

#Customised ggbiplot functions (different rotations)####
ggbiplot2<-function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                     obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                     ellipse = TRUE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                     alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                     varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                     color = muted("red"), # <- add new arguments to the function
                     linetype = "solid",
                     alpha_arrow = 1) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar*2), arrow = arrow(length = unit(1/2, 
                                                                                                    "picas")), color = color, linetype = linetype, alpha = alpha_arrow)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar*2, angle = angle, hjust = hjust), 
                       color = "black", size = varname.size)
  }
  return(g)
}

ggbiplot_reverse<-function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                            obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                            ellipse = TRUE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                            alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                            varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                            color = muted("red"), # <- add new arguments to the function
                            linetype = "solid",
                            alpha_arrow = 1) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/-xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(-xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar*2), arrow = arrow(length = unit(1/2, 
                                                                                                    "picas")), color = color, linetype = linetype, alpha = alpha_arrow)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar*2, angle = angle, hjust = hjust), 
                       color = "black", size = varname.size)
  }
  return(g)
}

ggbiplot_reverse_y<-function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                              obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                              ellipse = TRUE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                              alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                              varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                              color = muted("red"), # <- add new arguments to the function
                              linetype = "solid",
                              alpha_arrow = 1) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(-yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar*2), arrow = arrow(length = unit(1/2, 
                                                                                                    "picas")), color = color, linetype = linetype, alpha = alpha_arrow)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar*2, angle = angle, hjust = hjust), 
                       color = "black", size = varname.size)
  }
  return(g)
}

ggbiplot_reverse_xy<-function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                               obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                               ellipse = TRUE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                               alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                               varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                               color = muted("red"), # <- add new arguments to the function
                               linetype = "solid",
                               alpha_arrow = 1) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (225/pi) * atan(-yvar/-xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(-xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar*2), arrow = arrow(length = unit(1/2, 
                                                                                                    "picas")), color = color, linetype = linetype, alpha = alpha_arrow)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar*2, angle = angle, hjust = hjust), 
                       color = "black", size = varname.size)
  }
  return(g)
}





#2. Read in Data####

#Data is available in the following paper:

#Bjorkman, A.D. et al.,  2018. 
#Tundra Trait Team: A database of plant traits spanning the tundra biome. 
#Global Ecology and Biogeography.
#https://onlinelibrary.wiley.com/doi/10.1111/geb.12821

##try.ttt<-read.csv() #Read in data here
##try.ttt.original<-try.ttt #Keep original file for later

#3. Tranform data and calculate species means####
#Transform data
try.ttt$StdValue<-log10(try.ttt$StdValue)
try.ttt<-try.ttt[which((try.ttt)$StdValue!="-Inf"),]

#Standardise data
try.ttt<-try.ttt %>%
  group_by(TraitShort) %>%
  mutate(StdValueScaled = normalise(StdValue))

#Calculate species means
trait.meanstable<-ddply(try.ttt,.(TraitShort,AccSpeciesName),summarise,
                        mean = mean(StdValueScaled, na.rm=TRUE),
                        sd = sd(StdValueScaled, na.rm=TRUE),
                        se = se(StdValueScaled),
                        n= length(StdValueScaled))

trait.meanstable<- trait.meanstable[order(trait.meanstable$TraitShort),]

#Extract means
trait.meanstable<-trait.meanstable %>%
  select(TraitShort,AccSpeciesName,mean) %>%
  spread(TraitShort,mean)

#Add back in Functional Group
trait.meanstable$F_Group<-try.ttt$F_Group[match(trait.meanstable$AccSpeciesName,try.ttt$AccSpeciesName)]
rownames(trait.meanstable) <- trait.meanstable$AccSpeciesName

#Add missing groups
trait.meanstable[trait.meanstable$AccSpeciesName=="Calamagrostis purpurea",]$F_Group<-"GRAMINOID"
trait.meanstable[trait.meanstable$AccSpeciesName=="Rumex aquaticus",]$F_Group<-"FORBSV"

#4. Build permanova loops to calculate variance explained by groups####
#for all trait combinations

#How many species have data?
nspecies<-nrow(subset(trait.meanstable,complete.cases(trait.meanstable[,c(2:7)])))

rsq.table=NULL #Output table

traits<-colnames(trait.meanstable[,c(2:7)]) #Collect trait names

for(i in 2:6){ #Start at two as won't work for just one trait
  combinations<-combn(traits, i) #extract unique combination of traits for i samples
  combinations<-rbind(combinations,"F_Group") #Add F_Group column

  for(j in 1:ncol(combinations)){ #For number of unique combinations
    trait.data<-trait.meanstable[ , names(trait.meanstable) %in% as.vector(combinations[,j])] #Extract trait columns
    trait.data <- subset(trait.data,complete.cases(trait.data[,c(1:i)])) #Extract only complete cases

    #Draw random rows
    x=1 #Start repeats
    repeat {

      sample<-trait.data[sample(nrow(trait.data),nspecies), ] #Take maximum random samples

      mydata<-sample[,-(i+1)] #extract numerical data
      groups<-sample$F_Group #extract groups

      #Conduct permanova
      ad.run<-adonis(mydata ~ groups, permutations=999) #permanova - how much variation explained by groups

      #Extract variables
      rsq<-ad.run$aov.tab[1,5] #Extract R2
      nvar<-i #Note number of traits
      vars<-paste(unlist(colnames(mydata)), sep=",", collapse=", ") #Note trait identities

      #Combine in table
      rsq.table<-rbind(rsq.table,as.data.frame(cbind(nvar,rsq,vars))) #Combine

      x = x+1 #Break loop at required number of runs (1000)
      if (x == 11){ #Set number of runs (x+1)
        break
      }
    }
  }
}

#Clean output table

#Convert to numbers
rsq.table$rsq<-as.numeric(as.character(rsq.table$rsq))

rsq.table<-rsq.table %>%
  group_by(nvar,vars) %>%
  summarise (mean = mean(rsq),
             sd = sd(rsq),
             se = se(rsq))

rsq.table$vars<-as.character(rsq.table$vars)

#Process output table
rsq.table$w.LDMC<-apply(rsq.table,1,function(x) ifelse(grepl("LDMC",x[2]),1,0))
rsq.table$w.LA<-apply(rsq.table,1,function(x) ifelse(grepl("LeafArea",x[2]),1,0))
rsq.table$w.LN<-apply(rsq.table,1,function(x) ifelse(grepl("LeafN",x[2]),1,0))
rsq.table$w.Height<-apply(rsq.table,1,function(x) ifelse(grepl("PlantHeight",x[2]),1,0))
rsq.table$w.SLA<-apply(rsq.table,1,function(x) ifelse(grepl("SLA",x[2]),1,0))
rsq.table$w.SM<-apply(rsq.table,1,function(x) ifelse(grepl("SeedMass",x[2]),1,0))
rsq.table$w.SM.Height<-apply(rsq.table,1,function(x) ifelse(grepl("SeedMass",x[2])|grepl("PlantHeight",x[2]),1,0))
rsq.table$Type<-apply(rsq.table,1,function(x) ifelse(sum(c(as.numeric(x[6]),as.numeric(x[8]),as.numeric(x[10])))==0,"Structural",
                                                     ifelse(sum(c(as.numeric(x[7]),as.numeric(x[9]),as.numeric(x[11])))==0,"Economic","Both")))


#Plot change in explanatory power of functional group with trait information
(a<-ggplot(rsq.table,aes(nvar,mean*100,colour=factor(Type), shape=factor(Type)))+
   geom_point(cex=2.5)+
   theme_bw()+
   scale_shape_manual(values=c(1,19,19),
                      name="Trait Combination",
                      breaks = c("Economic","Structural","Both"),
                      labels=c("Only Economic Traits", "Only Morphological Traits","Economic & Morphological Traits"))+
   scale_colour_manual(values=c("black","deepskyblue2","red3"),
                       name="Trait Combination",
                       breaks = c("Economic","Structural","Both"),
                       labels=c("Only Economic Traits", "Only Morphological Traits","Economic & Morphological Traits"))+
   theme(legend.title=element_text(size=15) , legend.text=element_text(size=12), plot.title = element_text(size=18, vjust=1), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
   labs(x = "Number of traits", y = expression(paste("Variance Explained (%)")))+
   theme(axis.title.x = element_blank()))

(b<-ggplot(rsq.table,aes(nvar,mean*100,colour=factor(w.LDMC),shape=factor(w.LDMC)))+
   geom_point(cex=2.5)+
   theme_bw()+
   scale_shape_manual(values=c(1,19),
                      breaks = c(1,0),
                      name="Trait Combination",
                      labels=c("Includes LDMC", "Excludes LDMC"))+
   scale_colour_manual(values=c("black","deepskyblue2"),
                       breaks = c(1,0),
                       name="Trait Combination",
                       labels=c("Includes LDMC", "Excludes LDMC"))+
   theme(legend.title=element_text(size=15) , legend.text=element_text(size=12), plot.title = element_text(size=18, vjust=1), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
   labs(x = "Number of traits", y = expression(paste("Variance Explained (%)")))+
   theme(axis.title.x = element_blank()))

(c<-ggplot(rsq.table,aes(nvar,mean*100))+
    geom_point(cex=2.5,aes(colour=factor(w.SM.Height),shape=factor(w.SM.Height)))+
    theme_bw()+
    scale_shape_manual(values=c(1,19),
                       name="Trait Combination",
                       breaks = c(1,0),
                       labels=c("Includes Height / Seed Mass", "Excludes Height / Seed Mass"))+
    scale_colour_manual(values=c("black","red3"),
                        name="Trait Combination",
                        breaks = c(1,0),
                        labels=c("Includes Height / Seed Mass", "Excludes Height / Seed Mass"))+
    theme(legend.title=element_text(size=15) , legend.text=element_text(size=12), plot.title = element_text(size=18, vjust=1), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
    labs(x = "Number of traits", y = expression(paste("Variance Explained (%)"))))

#Reformat figure
gA <- ggplotGrob(a)
gB <- ggplotGrob(b)
gC <- ggplotGrob(c)

grid::grid.newpage()

grid::grid.draw(rbind(gA, gB, gC))

#5. Distribution of functional groups in trait-space####
#Part 1, functional groups only####
trait.meanstable_6<-trait.meanstable[,-8] #Remove wood density
trait.meanstable_6<-subset(trait.meanstable_6,complete.cases(trait.meanstable_6[,c(2:7)])) #Extract species with data for all traits

pca.data6<-trait.meanstable_6[,2:7] #trait data only
pca.species6<-trait.meanstable_6[,1] #species na,es
pca.groups6<-trait.meanstable_6[,8] #Functional groups
names(pca.data6)<-c("LDMC","LA","LN","PH","SM","SLA") #trait labels

#Conduct PCA
tundra.pca6 <- prcomp(pca.data6,
                      center = TRUE,
                      scale. = TRUE)

#Plot figure (no labels)
(groups.fig<-ggbiplot2(tundra.pca6, obs.scale = 1, var.scale = 1, ellipse = TRUE, 
                      circle = FALSE, alpha=0.9,group= pca.groups6, varname.size=0,varname.adjust = 2.4,color = "black")+
  scale_colour_manual(values=c("forestgreen","blue","darkorchid","orange","grey80"))+
  theme_bw()+
  theme(legend.position="none")+
  xlim(-6,6)+
  ylim(-5,5)+
  ggtitle("Traditional Functional Groups")+
  theme(plot.title = element_text(margin = margin(t = 10,r = 10,b = 10, l = 10)))+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")))

#Plot figure (axis labels)
(groups.fig.axes<-ggbiplot2(tundra.pca6, obs.scale = 1, var.scale = 1, ellipse = TRUE, 
                            circle = FALSE, alpha=0.8,group= pca.groups6, varname.size=6.5,varname.adjust = 3.5,color = "black")+
    scale_colour_manual(values=c("forestgreen","blue","darkorchid","orange","grey80"))+
    theme_bw()+
    theme(legend.position="none")+
    xlim(-6,6)+
    ylim(-6,6)+
    ggtitle("Whole Tundra (284 species)")+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(plot.title = element_text(size=15, vjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=15), axis.title=element_text(size=18,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")))

#Part 2. With k means####

rownames(trait.meanstable_6)<-trait.meanstable_6$AccSpeciesName #Set row names to species names
mydata <- trait.meanstable_6[c(-1,-8)] #trait data only

#Repeat groups to get best cluster 
#(k-means clustering method produces slightly different clusters each run due to clustering methodology:
#this approach uses the "best" cluster i.e. maximum variance explained)

x=1 #Start repeats
rsq<-NULL
rsq.table<-NULL
repeat {
  
  fit_kmeans <- kmeans(mydata, 4) # 4 cluster solution
  
  mydata_group <- data.frame(mydata, fit_kmeans$cluster)
  
  #Permanova
  ad.run<-adonis(mydata_group[,c(1:6)] ~ as.factor(mydata_group$fit_kmeans.cluster), permutations=999)
  
  rsq<-ad.run$aov.tab[1,5] #Extract R2
  
  #Combine in table
  rsq.table<-rbind(rsq.table,rsq)
  x = x+1 #Break loop at required number of runs (1000)
  if (x == 11){ #Set number of runs (x+1)
    break
  }
}

#Find max group
common.group<-max(rsq.table)

#Repeat clustering using "best" cluster to make sure working with this in
#later tests and figures

repeat {
  #k_means
  fit_kmeans<- kmeans(mydata, 4) # 4 cluster solution
  if (round(as.numeric(as.character(adonis(mydata_group[,c(1:6)] ~ as.factor(fit_kmeans$cluster), permutations=999)$aov.tab[1,5])),8) == round(common.group,8)){ #Break when right group
    break
  }
}

#append cluster assignment to table
mydata <- data.frame(mydata, fit_kmeans$cluster)
groups.kmean<-as.factor(fit_kmeans$cluster)

#Plot PCA
pca.mydata<-mydata[,c(1:6)]
pca.mydata.groups<-groups.kmean

mydata.pca <- prcomp(pca.mydata,
                     center = TRUE,
                     scale. = TRUE)

#NOTE: k-means clustering produces groups based on the data only, i.e. they do 
#not necessarily align with functional groups. The group numbering (1-4) is also 
#different for each run. These lines of code therefore match each group to the
#most closely aligned functional group. As numbering of k-means cluster
#differ each time, groups are hard coded based on species (have been tested for
#all group conbinations - see table S2)

#Name ESHRUB Group
es.group.k<-subset(pca.mydata.groups,names(pca.mydata.groups)=="Empetrum nigrum")
levels(pca.mydata.groups)[es.group.k]<-"ESHRUB"
#DSHRUB
es.group.k<-subset(pca.mydata.groups,names(pca.mydata.groups)=="Juncus trifidus")
levels(pca.mydata.groups)[es.group.k]<-"DSHRUB"
#GRAM
es.group.k<-subset(pca.mydata.groups,names(pca.mydata.groups)=="Betula glandulosa")
levels(pca.mydata.groups)[es.group.k]<-"GRAMINOID"
#FORB
es.group.k<-subset(pca.mydata.groups,names(pca.mydata.groups)=="Anthriscus sylvestris")
levels(pca.mydata.groups)[es.group.k]<-"FORBSV"

#Reorder factors
pca.mydata.groups<-factor(pca.mydata.groups, levels=c("ESHRUB","DSHRUB","GRAMINOID","FORBSV"))

#Plot figure
(kmeans.fig<-ggbiplot2(mydata.pca, obs.scale = 1, var.scale = 1, ellipse = TRUE, 
                       circle = FALSE, alpha=0.8,group= pca.mydata.groups,varname.size=0,varname.adjust = 2.4,color = "black")+
    scale_colour_manual(values=c("blue","forestgreen","orange","darkorchid","grey80"))+
    theme_bw()+
    theme(legend.position="none")+
    xlim(-6,6)+
    ylim(-5,5)+
    ggtitle("K-means")+
    theme(plot.title = element_text(margin = margin(t = 10,r = 10,b = 10, l = 10)))+
    theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")))

#Part 3. With Hierarchical Agglomerative####
mydata <- trait.meanstable_6[c(-1,-8)] #trait data only
rownames(mydata)<-trait.meanstable_6[,1] #Set rownames to species names

d <- dist(mydata, method = "euclidean") # distance matrix
fit_dendro <- hclust(d, method="ward.D") #Conduct clustering
groups <- cutree(fit_dendro, k=4) # cut tree into 4 clusters
mydata <- data.frame(mydata, groups) #Add groups to table

#Plot PCA
pca.mydata<-mydata[,c(1:6)]
groups.dendro<-as.factor(cutree(fit_dendro, k=4))
pca.mydata.groups<-groups.dendro

mydata.pca <- prcomp(pca.mydata,
                     center = TRUE,
                     scale. = TRUE)

#Name ESHRUB Group
es.group.k<-subset(pca.mydata.groups,names(pca.mydata.groups)=="Empetrum nigrum")
levels(pca.mydata.groups)[es.group.k]<-"ESHRUB"
#DSHRUB
es.group.k<-subset(pca.mydata.groups,names(pca.mydata.groups)=="Juncus trifidus")
levels(pca.mydata.groups)[es.group.k]<-"DSHRUB"
#GRAM
es.group.k<-subset(pca.mydata.groups,names(pca.mydata.groups)=="Betula glandulosa")
levels(pca.mydata.groups)[es.group.k]<-"GRAMINOID"
#FORB
es.group.k<-subset(pca.mydata.groups,names(pca.mydata.groups)=="Anthriscus sylvestris")
levels(pca.mydata.groups)[es.group.k]<-"FORBSV"

pca.mydata.groups<-factor(pca.mydata.groups, levels=c("ESHRUB","DSHRUB","GRAMINOID","FORBSV"))

#Plot figure
(dendro.fig<-ggbiplot2(mydata.pca, obs.scale = 1, var.scale = 1, ellipse = TRUE, 
                       circle = FALSE, alpha=0.8,group= pca.mydata.groups,varname.size=0,varname.adjust = 2.4,color = "black")+
    scale_colour_manual(values=c("blue","forestgreen","orange","darkorchid","grey80"))+
    theme_bw()+
    theme(legend.position="none")+
    xlim(-6,6)+
    ylim(-5,5)+
    ggtitle("Hierarchical Agglomerative")+
    theme(plot.title = element_text(margin = margin(t = 10,r = 10,b = 10, l = 10)))+
    theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")))

#Arrange left hand side of figure 5
grid.arrange(groups.fig, kmeans.fig, dendro.fig,nrow=3)

##6. Variance explained by groups####
#Part 1. For functional groups####
trait.meanstable.structure<-select(trait.meanstable_6,AccSpeciesName,LeafArea,PlantHeight,SeedMass,F_Group) #Select size traits
trait.meanstable.structure<-subset(trait.meanstable.structure,complete.cases(trait.meanstable.structure[,c(2:4)]))

trait.meanstable.economics<-select(trait.meanstable_6,AccSpeciesName,LDMC,SLA,LeafN,F_Group) #Select economic traits
trait.meanstable.economics<-subset(trait.meanstable.economics,complete.cases(trait.meanstable.economics[,c(2:4)]))

#First do with all traits
mydata<-trait.meanstable_6[,2:7]
pca.species6<-trait.meanstable_6[,1]
groups<-trait.meanstable_6[,8]

#Permanova
ad.run<-adonis(mydata ~ groups, permutations=999)
ad.run

rsq<-ad.run$aov.tab[1,5] #Extract R2
type<-"all"
cluster<-"fgroup"

#Combine in table
rsq.table.axes<-NULL
rsq.table.axes<-rbind(rsq.table.axes,as.data.frame(cbind(cluster,type,rsq))) #Combine

#Size traits
ad.run<-adonis(mydata ~ groups, permutations=999)
ad.run

rsq<-ad.run$aov.tab[1,5] #Extract R2
type<-"structure"
cluster<-"fgroup"

#Combine in table
rsq.table.axes<-rbind(rsq.table.axes,as.data.frame(cbind(cluster,type,rsq))) #Combine

#Economic traits
ad.run<-adonis(mydata ~ groups, permutations=999) #Run Permanova
ad.run
rsq<-ad.run$aov.tab[1,5] #Extract R2
type<-"economics"
cluster<-"fgroup"

#Combine in table
rsq.table.axes<-rbind(rsq.table.axes,as.data.frame(cbind(cluster,type,rsq))) #Combine

#Part 2. For k-means####
trait.meanstable.structure<-trait.meanstable_6
trait.meanstable.structure$F_Group<-as.factor(fit_kmeans$cluster)
trait.meanstable.structure<-select(trait.meanstable.structure,AccSpeciesName,LeafArea,PlantHeight,SeedMass,F_Group)
trait.meanstable.structure<-subset(trait.meanstable.structure,complete.cases(trait.meanstable.structure[,c(2:4)]))

trait.meanstable.economics<-trait.meanstable_6
trait.meanstable.economics$F_Group<-as.factor(fit_kmeans$cluster)
trait.meanstable.economics<-select(trait.meanstable.economics,AccSpeciesName,LDMC,SLA,LeafN,F_Group)
trait.meanstable.economics<-subset(trait.meanstable.economics,complete.cases(trait.meanstable.economics[,c(2:4)]))

#First do with all traits
mydata<-trait.meanstable_6[,2:7]
pca.species6<-trait.meanstable_6[,1]
groups.kmean<-as.factor(fit_kmeans$cluster)

#Permanova
ad.run<-adonis(mydata ~ groups.kmean, permutations=999)
ad.run

rsq<-ad.run$aov.tab[1,5] #Extract R2
type<-"all"
cluster<-"kmeans"

#Combine in table
rsq.table.axes<-rbind(rsq.table.axes,as.data.frame(cbind(cluster,type,rsq))) #Combine

#Size traits
ad.run<-adonis(mydata ~ groups.kmean, permutations=999)
ad.run

rsq<-ad.run$aov.tab[1,5] #Extract R2
type<-"structure"
cluster<-"kmeans"

#Combine in table
rsq.table.axes<-rbind(rsq.table.axes,as.data.frame(cbind(cluster,type,rsq))) #Combine

#Economic traits
ad.run<-adonis(mydata ~ groups.kmean, permutations=999) #Run Permanova
ad.run
rsq<-ad.run$aov.tab[1,5] #Extract R2
type<-"economics"
cluster<-"kmeans"

#Combine in table
rsq.table.axes<-rbind(rsq.table.axes,as.data.frame(cbind(cluster,type,rsq))) #Combine

#Part 3. For HA clusters####
trait.meanstable.structure<-trait.meanstable_6
trait.meanstable.structure$F_Group<-as.factor(cutree(fit_dendro, k=4))
trait.meanstable.structure<-select(trait.meanstable.structure,AccSpeciesName,LeafArea,PlantHeight,SeedMass,F_Group)
trait.meanstable.structure<-subset(trait.meanstable.structure,complete.cases(trait.meanstable.structure[,c(2:4)]))

trait.meanstable.economics<-trait.meanstable_6
trait.meanstable.economics$F_Group<-as.factor(cutree(fit_dendro, k=4))
trait.meanstable.economics<-select(trait.meanstable.economics,AccSpeciesName,LDMC,SLA,LeafN,F_Group)
trait.meanstable.economics<-subset(trait.meanstable.economics,complete.cases(trait.meanstable.economics[,c(2:4)]))

#First do with all traits
mydata<-trait.meanstable_6[,2:7]
pca.species6<-trait.meanstable_6[,1]
groups.dendro<-as.factor(cutree(fit_dendro, k=4))

#Permanova
ad.run<-adonis(mydata ~ groups.dendro, permutations=999)
ad.run

rsq<-ad.run$aov.tab[1,5] #Extract R2
type<-"all"
cluster<-"dendro"

#Combine in table
rsq.table.axes<-rbind(rsq.table.axes,as.data.frame(cbind(cluster,type,rsq))) #Combine

#Size traits
ad.run<-adonis(mydata ~ groups.dendro, permutations=999)
ad.run

rsq<-ad.run$aov.tab[1,5] #Extract R2
type<-"structure"
cluster<-"dendro"

#Combine in table
rsq.table.axes<-rbind(rsq.table.axes,as.data.frame(cbind(cluster,type,rsq))) #Combine

#Economic traits
ad.run<-adonis(mydata ~ groups.dendro, permutations=999) #Run Permanova
ad.run
rsq<-ad.run$aov.tab[1,5] #Extract R2
type<-"economics"
cluster<-"dendro"

#Combine in table
rsq.table.axes<-rbind(rsq.table.axes,as.data.frame(cbind(cluster,type,rsq))) #Combine

#Convert units
rsq.table.axes$rsq<-as.numeric(as.character(rsq.table.axes$rsq))
rsq.table.axes$type <- factor(rsq.table.axes$type, levels=c("economics","structure","all"))


#Plot middle column of Figure 5
groups.explained.fig<-ggplot(rsq.table.axes[rsq.table.axes$cluster=="fgroup",],aes(type,rsq*100))+
  geom_bar(stat = "identity",colour="black",fill=c("deepskyblue2","red3","white"),width=0.8)+
  coord_flip()+
  theme_bw()+
  ylim(0,100)+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(margin = margin(t = 10,r = 10,b = 10, l = 10)))+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
  labs(y="Variance Explained (%)",x="Trait Type")

kmeans.explained.fig<-ggplot(rsq.table.axes[rsq.table.axes$cluster=="kmeans",],aes(type,rsq*100))+
  geom_bar(stat = "identity",colour="black",fill=c("deepskyblue2","red3","white"),width=0.8)+
  coord_flip()+
  theme_bw()+
  ylim(0,100)+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(margin = margin(t = 10,r = 10,b = 10, l = 10)))+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
  labs(y="Variance Explained (%)",x="Trait Type")

dendro.explained.fig<-ggplot(rsq.table.axes[rsq.table.axes$cluster=="dendro",],aes(type,rsq*100))+
  geom_bar(stat = "identity",colour="black",fill=c("deepskyblue2","red3","white"),width=0.8)+
  coord_flip()+
  theme_bw()+
  ylim(0,100)+
  ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(margin = margin(t = 10,r = 10,b = 10, l = 10)))+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
  labs(y="Variance Explained (%)",x="Trait Type")

grid.arrange(groups.explained.fig,kmeans.explained.fig,dendro.explained.fig,nrow=3)

#7. Change in species composition among groups####

#Create table with all groups
trait.meanstable_6$kmeans<-as.numeric(as.character(groups.kmean))
trait.meanstable_6$dendro<-as.numeric(as.character(groups.dendro))

#Name ESHRUB Group
#Find group
es.group.k<-trait.meanstable_6[trait.meanstable_6$AccSpeciesName=="Empetrum nigrum",]$kmeans
es.group.d<-trait.meanstable_6[trait.meanstable_6$AccSpeciesName=="Empetrum nigrum",]$dendro

#Rename 'group'
trait.meanstable_6[trait.meanstable_6$kmeans==es.group.k,]$kmeans<-"ESHRUB"
trait.meanstable_6[trait.meanstable_6$dendro==es.group.d,]$dendro<-"ESHRUB"

#Name DSHRUB Group
#Find group
ds.group.k<-trait.meanstable_6[trait.meanstable_6$AccSpeciesName=="Juncus trifidus",]$kmeans
ds.group.d<-trait.meanstable_6[trait.meanstable_6$AccSpeciesName=="Juncus trifidus",]$dendro

#Rename 'group'
trait.meanstable_6[trait.meanstable_6$kmeans==ds.group.k,]$kmeans<-"DSHRUB"
trait.meanstable_6[trait.meanstable_6$dendro==ds.group.d,]$dendro<-"DSHRUB"

#Name GRAM Group
#Find group
gr.group.k<-trait.meanstable_6[trait.meanstable_6$AccSpeciesName=="Betula glandulosa",]$kmeans
gr.group.d<-trait.meanstable_6[trait.meanstable_6$AccSpeciesName=="Betula glandulosa",]$dendro

#Rename 'group'
trait.meanstable_6[trait.meanstable_6$kmeans==gr.group.k,]$kmeans<-"GRAMINOID"
trait.meanstable_6[trait.meanstable_6$dendro==gr.group.d,]$dendro<-"GRAMINOID"

#Name FORBSV Group
#Find group
fb.group.k<-trait.meanstable_6[trait.meanstable_6$AccSpeciesName=="Anthriscus sylvestris",]$kmeans
fb.group.d<-trait.meanstable_6[trait.meanstable_6$AccSpeciesName=="Anthriscus sylvestris",]$dendro

#Rename 'group'
trait.meanstable_6[trait.meanstable_6$kmeans==fb.group.k,]$kmeans<-"FORBSV"
trait.meanstable_6[trait.meanstable_6$dendro==fb.group.d,]$dendro<-"FORBSV"

trait.meanstable_6$F_Group<-as.character(trait.meanstable_6$F_Group)
trait.meanstable_6$kmeans<-as.character(trait.meanstable_6$kmeans)
trait.meanstable_6$dendro<-as.character(trait.meanstable_6$dendro)

#extract only cluster information, ready to colour rows
trait.meanstable_6_colour <- trait.meanstable_6 %>%
  select(F_Group,kmeans,dendro) %>%
  mutate(colour = rownames(.)) %>%
  melt(id.vars="colour")

#Set species order (same for each cluster)
trait.meanstable_6$F_Group <- factor(trait.meanstable_6$F_Group, levels=c("FORBSV", "GRAMINOID", "DSHRUB", "ESHRUB"))
trait.meanstable_6$kmeans <- factor(trait.meanstable_6$kmeans, levels=c("FORBSV", "GRAMINOID", "DSHRUB", "ESHRUB"))
trait.meanstable_6$dendro <- factor(trait.meanstable_6$dendro, levels=c("FORBSV", "GRAMINOID", "DSHRUB", "ESHRUB"))

#Rearrange table for plotting
species<-trait.meanstable_6 %>%
  arrange(dendro) %>%
  arrange(kmeans) %>%
  arrange(F_Group)

species<-species[, "AccSpeciesName"]

species<-factor(species, levels = species)

trait.meanstable_6_colour$colour<-factor(trait.meanstable_6_colour$colour, levels=species)


(group_plot<-ggplot(trait.meanstable_6_colour, aes(x=variable, y=colour, label=variable, fill=as.factor(value))) + 
    geom_tile()+
    theme_bw()+theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(axis.ticks = element_blank(), axis.text.y = element_blank())+
    labs(x = "\nGrouping Method", y = "Species")+
    geom_vline(xintercept = 1.5, linetype = "solid")+
    geom_vline(xintercept = 2.5, linetype = "solid")+  
    scale_x_discrete(labels=c("F_Group" = "Functional\n Groups", "kmeans" = "K-means",
                              "dendro" = "Hierarchical\n Agglomerative"))+
    scale_fill_manual(values=alpha(c("forestgreen","blue","darkorchid","orange"),0.7),
                      name="Functional Group\n",
                      breaks=c("ESHRUB", "DSHRUB", "GRAMINOID", "FORBSV"),
                      labels=c("Evergreen Shrub", "Deciduous Shrub", "Graminoid", "Forb"))+
    theme(legend.title=element_text(size=15) , legend.text=element_text(size=12))+
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14)))


#8. Plot all panels of Figure 5####
grid.arrange(groups.fig, groups.explained.fig, kmeans.fig, kmeans.explained.fig,dendro.fig,dendro.explained.fig,nrow=3)

gA <- ggplotGrob(groups.fig)
gB <- ggplotGrob(groups.explained.fig)

gC <- ggplotGrob(kmeans.fig)
gD <- ggplotGrob(kmeans.explained.fig)

gE <- ggplotGrob(dendro.fig)
gF <- ggplotGrob(dendro.explained.fig)

gH <- ggplotGrob(group_plot)

#Left two columns
grid::grid.newpage()
grid::grid.draw(cbind(rbind(gA, gC, gE),rbind(gB, gD,gF)))

#Right column
grid::grid.newpage()
grid::grid.draw(gH)

#9 Calculate similarity between groups####
trait.meanstable_6_compare<-trait.meanstable_6[,c(8:10)]

trait.meanstable_6_compare$Ftok<-ifelse(trait.meanstable_6_compare$F_Group == trait.meanstable_6_compare$kmeans,1,0) #Functional vs k-means
trait.meanstable_6_compare$Ftod<-ifelse(trait.meanstable_6_compare$F_Group == trait.meanstable_6_compare$dendro,1,0) #Functional vs HA
trait.meanstable_6_compare$dtok<-ifelse(trait.meanstable_6_compare$dendro == trait.meanstable_6_compare$kmeans,1,0) #k-means vs HA
trait.meanstable_6_compare$all<-ifelse(trait.meanstable_6_compare$F_Group == trait.meanstable_6_compare$kmeans & trait.meanstable_6_compare$F_Group == trait.meanstable_6_compare$dendro,1,0) #All clusters


#Calculate metrics
#Overall
ftok_tundra<-sum(trait.meanstable_6_compare$Ftok)/nrow(trait.meanstable_6_compare) #tundra, all
ftod_tundra<-sum(trait.meanstable_6_compare$Ftod)/nrow(trait.meanstable_6_compare) #tundra, all
dtok_tundra<-sum(trait.meanstable_6_compare$dtok)/nrow(trait.meanstable_6_compare) #tundra, all
all_tundra<-sum(trait.meanstable_6_compare$all)/nrow(trait.meanstable_6_compare) #tundra, all

#ESHRUB
ftok_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="ESHRUB",]$Ftok)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="ESHRUB",]) #tundra, E shrubs
ftod_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="ESHRUB",]$Ftod)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="ESHRUB",]) #tundra, E shrubs
dtok_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="ESHRUB",]$dtok)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="ESHRUB",]) #tundra, E shrubs
all_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="ESHRUB",]$all)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="ESHRUB",]) #tundra, E shrubs

#DSHRUB
ftok_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="DSHRUB",]$Ftok)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="DSHRUB",]) #tundra, D shrubs
ftod_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="DSHRUB",]$Ftod)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="DSHRUB",]) #tundra, D shrubs
dtok_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="DSHRUB",]$dtok)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="DSHRUB",]) #tundra, D shrubs
all_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="DSHRUB",]$all)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="DSHRUB",]) #tundra, D shrubs

#GRAMINOID
ftok_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="GRAMINOID",]$Ftok)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="GRAMINOID",]) #tundra, Graminoids
ftod_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="GRAMINOID",]$Ftod)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="GRAMINOID",]) #tundra, Graminoids
dtok_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="GRAMINOID",]$dtok)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="GRAMINOID",]) #tundra, Graminoids
all_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="GRAMINOID",]$all)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="GRAMINOID",]) #tundra, Graminoids

#FORBSV
ftok_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="FORBSV",]$Ftok)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="FORBSV",]) #tundra, Forbs
ftod_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="FORBSV",]$Ftod)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="FORBSV",]) #tundra, Forbs
dtok_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="FORBSV",]$dtok)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="FORBSV",]) #tundra, Forbs
all_tundra<-sum(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="FORBSV",]$all)/nrow(trait.meanstable_6_compare[trait.meanstable_6_compare$F_Group=="FORBSV",]) #tundra, Forbs

#Identify consistent species
consistent_species<-subset(trait.meanstable_6_compare,all==1) #Species names
groups.consistent<-as.factor(trait.meanstable_6_compare$all) #Grouping factor

#Plot PCA based on whether species are consistently categorised
(ggbiplot2(mydata.pca, obs.scale = 1, var.scale = 1, ellipse = TRUE, 
          circle = FALSE, alpha=0.8,group= groups.consistent,varname.size=4.5,varname.adjust = 1.5,color = "black")+
  scale_colour_manual(values=c("orange","darkorchid","forestgreen","blue","grey80"))+
  theme(legend.position="none")+
  ylim(-6,6)+
  xlim(-6,6)+
  theme_bw()+
  theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=10), axis.title=element_text(size=12,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm")))

#10. Single trait distribution figures####

#Draw for all traits
#Reload original (untransformed) trait data
try.ttt.original$F_Group<-factor(try.ttt.original$F_Group, levels = c("DSHRUB", "ESHRUB", "FORBSV","GRAMINOID"))

#Extract only species used in main analysis
species_used<-trait.meanstable_6$AccSpeciesName

try.ttt.6<-subset(try.ttt.original,AccSpeciesName %in% species_used)

#Create new meanstable without log-transforming data
orig.meanstable<-try.ttt.6 %>%
  group_by(AccSpeciesName,TraitShort) %>%
  summarise(StdValue = mean(StdValue)) %>%
  ungroup()

#Add kmeans & HA clusters
orig.meanstable$F_Group<-trait.meanstable_6$F_Group[match(orig.meanstable$AccSpeciesName,trait.meanstable_6$AccSpeciesName)]
orig.meanstable$kmeans<-trait.meanstable_6$kmeans[match(orig.meanstable$AccSpeciesName,trait.meanstable_6$AccSpeciesName)]
orig.meanstable$dendro<-trait.meanstable_6$dendro[match(orig.meanstable$AccSpeciesName,trait.meanstable_6$AccSpeciesName)]

orig.meanstable$F_Group<-factor(orig.meanstable$F_Group, levels = c("DSHRUB", "ESHRUB", "FORBSV","GRAMINOID"))

#Height
pairwise.wilcox.test(orig.meanstable[orig.meanstable$TraitShort=="PlantHeight",]$StdValue,orig.meanstable[orig.meanstable$TraitShort=="PlantHeight",]$F_Group)

(height<-ggplot(orig.meanstable[orig.meanstable$TraitShort=="PlantHeight",],aes(StdValue,fill=F_Group))+
    geom_density(adjust=2)+
    scale_x_log10(limits = c(0.005,1.5)) +
    scale_fill_manual(values=alpha(c("forestgreen","blue","darkorchid","orange"),0.6),
                      name="Functional Group\n",
                      breaks=c("DSHRUB", "ESHRUB", "GRAMINOID", "FORBSV"),
                      labels=c("DeciduousShrub", "Evergreen Shrub", "Graminoid", "Forb"))+
    theme_bw()+
    ylim(0,1.2)+
    theme(plot.title = element_text(size=18, vjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(0.4,0.4,0.4,0.4), units = , "cm"))+
    labs(x = "Plant Height (m)", y = "Density"))

#SLA
pairwise.wilcox.test(orig.meanstable[orig.meanstable$TraitShort=="SLA",]$StdValue,orig.meanstable[orig.meanstable$TraitShort=="SLA",]$F_Group)

(SLA<-ggplot(orig.meanstable[orig.meanstable$TraitShort=="SLA",],aes(StdValue,fill=F_Group))+
    geom_density(adjust=2)+
    scale_x_log10(limits=c(1,100)) +
    scale_fill_manual(values=alpha(c("forestgreen","blue","darkorchid","orange"),0.6),
                      name="Functional Group\n",
                      breaks=c("DSHRUB", "ESHRUB", "GRAMINOID", "FORBSV"),
                      labels=c("Deciduous Shrub", "Evergreen Shrub", "Graminoid", "Forb"))+
    theme_bw()+
    ylim(0,4.5)+
    theme(plot.title = element_text(size=18, vjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(0.4,0.4,0.4,0.4), units = , "cm"))+
    labs(x = expression(paste("Specific Leaf Area (",mm^{2}," ",mg^{-1},")")), y = "Density"))

#Leaf Area
pairwise.wilcox.test(orig.meanstable[orig.meanstable$TraitShort=="LeafArea",]$StdValue,orig.meanstable[orig.meanstable$TraitShort=="LeafArea",]$F_Group)

(LA<-ggplot(orig.meanstable[orig.meanstable$TraitShort=="LeafArea",],aes(StdValue,fill=F_Group))+
    geom_density(adjust=2)+
    scale_x_log10(limits = c(0.1,100000)) +
    scale_fill_manual(values=alpha(c("forestgreen","blue","darkorchid","orange"),0.6),
                      name="Functional Group\n",
                      breaks=c("DSHRUB", "ESHRUB", "GRAMINOID", "FORBSV"),
                      labels=c("Deciduous Shrub", "Evergreen Shrub", "Graminoid", "Forb"))+
    theme_bw()+
    ylim(0,0.8)+
    theme(plot.title = element_text(size=18, vjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(0.4,0.4,0.4,0.4), units = , "cm"))+
    labs(x = expression(paste("Leaf Area (",mm^{2},")")), y = "Density"))

#Leaf nitrogen
pairwise.wilcox.test(orig.meanstable[orig.meanstable$TraitShort=="LeafN",]$StdValue,orig.meanstable[orig.meanstable$TraitShort=="LeafN",]$F_Group)

(LN<-ggplot(orig.meanstable[orig.meanstable$TraitShort=="LeafN",],aes(StdValue,fill=F_Group))+
    geom_density(adjust=2)+
    scale_x_log10(limits = c(3,100)) +
    scale_fill_manual(values=alpha(c("forestgreen","blue","darkorchid","orange"),0.6),
                      name="Functional Group\n",
                      breaks=c("DSHRUB", "ESHRUB", "GRAMINOID", "FORBSV"),
                      labels=c("Deciduous Shrub", "Evergreen Shrub", "Graminoid", "Forb"))+
    theme_bw()+
    ylim(0,3.2)+
    theme(plot.title = element_text(size=18, vjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(0.4,0.4,0.4,0.4), units = , "cm"))+
    labs(x = expression(paste("Leaf N (mg ",g^{-1},")")), y = "Density"))

#LDMC
pairwise.wilcox.test(orig.meanstable[orig.meanstable$TraitShort=="LDMC",]$StdValue,orig.meanstable[orig.meanstable$TraitShort=="LDMC",]$F_Group)

(LDMC<-ggplot(orig.meanstable[orig.meanstable$TraitShort=="LDMC",],aes(StdValue,fill=F_Group))+
    geom_density(adjust=2)+
    scale_x_log10(limits = c(0.07,1)) +
    scale_fill_manual(values=alpha(c("forestgreen","blue","darkorchid","orange"),0.6),
                      name="Functional Group\n",
                      breaks=c("DSHRUB", "ESHRUB", "GRAMINOID", "FORBSV"),
                      labels=c("Deciduous Shrub", "Evergreen Shrub", "Graminoid", "Forb"))+
    theme_bw()+
    ylim(0,4.75)+
    theme(plot.title = element_text(size=18, vjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(0.4,0.4,0.4,0.4), units = , "cm"))+
    labs(x = expression(paste("Leaf Dry Matter Content (g ",g^{-1},")")), y = "Density"))

#Seed mass
pairwise.wilcox.test(orig.meanstable[orig.meanstable$TraitShort=="SeedMass",]$StdValue,orig.meanstable[orig.meanstable$TraitShort=="SeedMass",]$F_Group)

(SM<-ggplot(orig.meanstable[orig.meanstable$TraitShort=="SeedMass",],aes(StdValue,fill=F_Group))+
    geom_density(adjust=4)+
    scale_x_log10(limits = c(0.0001,1000)) +
    scale_fill_manual(values=alpha(c("forestgreen","blue","darkorchid","orange"),0.6),
                      name="Functional Group\n",
                      breaks=c("DSHRUB", "ESHRUB", "GRAMINOID", "FORBSV"),
                      labels=c("Deciduous Shrub", "Evergreen Shrub", "Graminoid", "Forb"))+
    theme_bw()+
    ylim(0,0.8)+
    theme(plot.title = element_text(size=18, vjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(0.4,0.4,0.4,0.4), units = , "cm"))+
    labs(x = "Seed Mass (mg)", y = "Density"))

#Combine figures
mylegend<-g_legend(height) #extract legend

grid.arrange(arrangeGrob(height + theme(legend.position="none"),
                         LA + theme(legend.position="none"),
                         SM + theme(legend.position="none"),
                         SLA + theme(legend.position="none"),
                         LDMC + theme(legend.position="none"),
                         LN + theme(legend.position="none"),
                         nrow=3))

grid.arrange(mylegend)


#Draw with all trait observations, not species means
#Add kmeans & HA
try.ttt.6$F_Group<-trait.meanstable_6$F_Group[match(try.ttt.6$AccSpeciesName,trait.meanstable_6$AccSpeciesName)]
try.ttt.6$kmeans<-trait.meanstable_6$kmeans[match(try.ttt.6$AccSpeciesName,trait.meanstable_6$AccSpeciesName)]
try.ttt.6$dendro<-trait.meanstable_6$dendro[match(try.ttt.6$AccSpeciesName,trait.meanstable_6$AccSpeciesName)]

#Height
pairwise.wilcox.test(orig.meanstable[orig.meanstable$TraitShort=="PlantHeight",]$StdValue,orig.meanstable[orig.meanstable$TraitShort=="PlantHeight",]$F_Group)

(height_all<-ggplot(try.ttt.6[try.ttt.6$TraitShort=="PlantHeight",],aes(StdValue,fill=F_Group))+
    geom_density(adjust=3)+
    scale_x_log10(limits = c(0.005,1.5)) +
    scale_fill_manual(values=alpha(c("forestgreen","blue","darkorchid","orange"),0.6),
                      name="Functional Group\n",
                      breaks=c("ESHRUB", "DSHRUB", "GRAMINOID", "FORBSV"),
                      labels=c("Evergreen Shrub", "Deciduous Shrub", "Graminoid", "Forb"))+
    theme_bw()+
    ylim(0,1.2)+
    theme(plot.title = element_text(size=18, vjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(0.4,0.4,0.4,0.4), units = , "cm"))+
    labs(x = "Plant Height (m)", y = "Density"))

#SLA
pairwise.wilcox.test(orig.meanstable[orig.meanstable$TraitShort=="SLA",]$StdValue,orig.meanstable[orig.meanstable$TraitShort=="SLA",]$F_Group)

(SLA_all<-ggplot(try.ttt.6[try.ttt.6$TraitShort=="SLA",],aes(StdValue,fill=F_Group))+
    geom_density(adjust=2)+
    scale_x_log10(limits=c(1,100)) +
    scale_fill_manual(values=alpha(c("forestgreen","blue","darkorchid","orange"),0.6),
                      name="Functional Group\n",
                      breaks=c("ESHRUB", "DSHRUB", "GRAMINOID", "FORBSV"),
                      labels=c("Evergreen Shrub", "Deciduous Shrub", "Graminoid", "Forb"))+
    theme_bw()+
    ylim(0,4.5)+
    theme(plot.title = element_text(size=18, vjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(0.4,0.4,0.4,0.4), units = , "cm"))+
    labs(x = expression(paste("Specific Leaf Area (",mm^{2}," ",mg^{-1},")")), y = "Density"))

#Leaf Area
pairwise.wilcox.test(orig.meanstable[orig.meanstable$TraitShort=="LeafArea",]$StdValue,orig.meanstable[orig.meanstable$TraitShort=="LeafArea",]$F_Group)

(LA_all<-ggplot(try.ttt.6[try.ttt.6$TraitShort=="LeafArea",],aes(StdValue,fill=F_Group))+
    geom_density(adjust=3)+
    scale_x_log10(limits = c(0.1,100000)) +
    scale_fill_manual(values=alpha(c("forestgreen","blue","darkorchid","orange"),0.6),
                      name="Functional Group\n",
                      breaks=c("ESHRUB", "DSHRUB", "GRAMINOID", "FORBSV"),
                      labels=c("Evergreen Shrub", "Deciduous Shrub", "Graminoid", "Forb"))+
    theme_bw()+
    ylim(0,0.8)+
    theme(plot.title = element_text(size=18, vjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(0.4,0.4,0.4,0.4), units = , "cm"))+
    labs(x = expression(paste("Leaf Area (",mm^{2},")")), y = "Density"))

#Leaf nitrogen
pairwise.wilcox.test(orig.meanstable[orig.meanstable$TraitShort=="LeafN",]$StdValue,orig.meanstable[orig.meanstable$TraitShort=="LeafN",]$F_Group)

(LN_all<-ggplot(try.ttt.6[try.ttt.6$TraitShort=="LeafN",],aes(StdValue,fill=F_Group))+
    geom_density(adjust=2)+
    scale_x_log10(limits = c(3,100)) +
    scale_fill_manual(values=alpha(c("forestgreen","blue","darkorchid","orange"),0.6),
                      name="Functional Group\n",
                      breaks=c("ESHRUB", "DSHRUB", "GRAMINOID", "FORBSV"),
                      labels=c("Evergreen Shrub", "Deciduous Shrub", "Graminoid", "Forb"))+
    theme_bw()+
    ylim(0,3.2)+
    theme(plot.title = element_text(size=18, vjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(0.4,0.4,0.4,0.4), units = , "cm"))+
    labs(x = expression(paste("Leaf N (mg ",g^{-1},")")), y = "Density"))

#LDMC
pairwise.wilcox.test(orig.meanstable[orig.meanstable$TraitShort=="LDMC",]$StdValue,orig.meanstable[orig.meanstable$TraitShort=="LDMC",]$F_Group)

(LDMC_all<-ggplot(try.ttt.6[try.ttt.6$TraitShort=="LDMC",],aes(StdValue,fill=F_Group))+
    geom_density(adjust=2)+
    scale_x_log10(limits = c(0.07,1)) +
    scale_fill_manual(values=alpha(c("forestgreen","blue","darkorchid","orange"),0.6),
                      name="Functional Group\n",
                      breaks=c("ESHRUB", "DSHRUB", "GRAMINOID", "FORBSV"),
                      labels=c("Evergreen Shrub", "Deciduous Shrub", "Graminoid", "Forb"))+
    theme_bw()+
    ylim(0,4.75)+
    theme(plot.title = element_text(size=18, vjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(0.4,0.4,0.4,0.4), units = , "cm"))+
    labs(x = expression(paste("Leaf Dry Matter Content (g ",g^{-1},")")), y = "Density"))

#Seed Mass
pairwise.wilcox.test(orig.meanstable[orig.meanstable$TraitShort=="SeedMass",]$StdValue,orig.meanstable[orig.meanstable$TraitShort=="SeedMass",]$F_Group)

(SM_all<-ggplot(try.ttt.6[try.ttt.6$TraitShort=="SeedMass",],aes(StdValue,fill=F_Group))+
    geom_density(adjust=4)+
    scale_x_log10(limits = c(0.0001,1000)) +
    scale_fill_manual(values=alpha(c("forestgreen","blue","darkorchid","orange"),0.6),
                      name="Functional Group\n",
                      breaks=c("ESHRUB", "DSHRUB", "GRAMINOID", "FORBSV"),
                      labels=c("Evergreen Shrub", "Deciduous Shrub", "Graminoid", "Forb"))+
    theme_bw()+
    ylim(0,0.8)+
    theme(plot.title = element_text(size=18, vjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(0.4,0.4,0.4,0.4), units = , "cm"))+
    labs(x = "Seed Mass (mg)", y = "Density"))

#Combine figures
grid.arrange(arrangeGrob(height_all + theme(legend.position="none"),
                         LA_all + theme(legend.position="none"),
                         SM_all + theme(legend.position="none"),
                         SLA_all + theme(legend.position="none"),
                         LDMC_all + theme(legend.position="none"),
                         LN_all + theme(legend.position="none"),
                         nrow=3))

# Code for supplementary analyses is not included here since this 
# makes use of additional tundra composition data. For detail on supplementary
# analyses please contact haydn.thomas@ed.ac.uk