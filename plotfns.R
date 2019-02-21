# green.shade
# red.shade
# rg.array.colors       Red/green array.
# by.array.colors       Blue/yellow array.
# rgb.colors            From red, green, blue colorwheel.
# ryb.colors            From red, yellow, blue colorwheel.
# matlab.colors         Like Matlab.
# yahoo.weather.colors  Like Yahoo Weather Map.
# genespring.colors
# 
# my.lineplot
# my.heatmap
# my.colorbar
# sample.colorbars      Print out all different color schemes.
#
# INTERNAL FUNCTIONS
# colorbar2rgb
# breaks2color
# wheel2bar


green.shade <- function(n) {
  sapply((0:(n-1))/(n-1), function(x) { rgb(0, x, 0) } )
}

red.shade <- function(n) {
  sapply((0:(n-1))/(n-1), function(x) { rgb(x, 0, 0) } )
}

rg.array.colors <- function(n) {
  colors <- c()
  shades <- ((0:(n-1)/(n-1) - 0.5)) * 2
  for(s in shades) {
    if(s >= 0)
      x <- rgb(s, 0, 0)
    else
      x <- rgb(0, -s, 0)
    colors <- c(colors, x)
  }
  colors
}

by.array.colors <- function(n) {
  colors <- c()
  shades <- ((0:(n-1)/(n-1) - 0.5)) * 2
  for(s in shades) {
    if(s >= 0)
      x <- rgb(s, s, 0)
    else
      x <- rgb(0, 0, -s)
    colors <- c(colors, x)
  }
  colors
}

rgb.colors <- function(n, shade=0, tint=0,
  min.r=0, max.r=1, min.g=0, max.g=1, min.b=0, max.b=1) {
  COLORWHEEL.RGB <- list(
    RED=function(pos) {
      breaks <- c(255, 255, 0, 0, 0, 255, 255)/255
      breaks2color(breaks, pos)
    },
    GREEN=function(pos) {
      breaks <- c(0, 255, 255, 255, 0, 0, 0)/255
      breaks2color(breaks, pos)
    },
    BLUE=function(pos) {
      breaks <- c(0, 0, 0, 255, 255, 255, 0)/255
      breaks2color(breaks, pos)
    }
  )
  COLORBAR.RGB <- wheel2bar(COLORWHEEL.RGB)
  I <- 0:(n-1) / (n-1)
  I <- rev(I)   # make red the brightest, blue darkest.
  sapply(I, function(x)
    colorbar2rgb(x, colorbar=COLORBAR.RGB,
      shade=shade, tint=tint, min.r=min.r, max.r=max.r, 
      min.g=min.g, max.g=max.g, min.b=min.b, max.b=max.b) 
  )
}

ryb.colors <- function(n, shade=0, tint=0,
  min.r=0, max.r=1, min.g=0, max.g=1, min.b=0, max.b=1) {
  COLORWHEEL.RYB <- list(
    RED=function(pos) {
      breaks <- c(255, 255, 255, 0, 0, 255, 255)/255
      breaks2color(breaks, pos)
    },
    GREEN=function(pos) {
      breaks <- c(0, 128, 255, 255, 0, 0, 0)/255
      breaks2color(breaks, pos)
    },
    BLUE=function(pos) {
      breaks <- c(0, 0, 0, 0, 255, 255, 0)/255
      breaks2color(breaks, pos)
    }
  )
  COLORBAR.RYB <- wheel2bar(COLORWHEEL.RYB)
  I <- 0:(n-1) / (n-1)
  I <- rev(I)   # make red the brightest, blue darkest.
  sapply(I, function(x)
    colorbar2rgb(x, colorbar=COLORBAR.RYB,
      shade=shade, tint=tint, min.r=min.r, max.r=max.r, 
      min.g=min.g, max.g=max.g, min.b=min.b, max.b=max.b)
  )
}

matlab.colors <- function(n, shade=0, tint=0,
  min.r=0, max.r=1, min.g=0, max.g=1, min.b=0, max.b=1) {
  COLORBAR.MAT <- list(
    RED=function(pos) {
      breaks <- c(127, 255, 255, 255, 127,   0,   0,   0,   0)/255
      breaks2color(breaks, pos)
    },
    GREEN=function(pos) {
      breaks <- c(  0,   0, 127, 255, 255, 255, 127,   0,   0)/255
      breaks2color(breaks, pos)
    },
    BLUE=function(pos) {
      breaks <- c(  0,   0,   0,   0, 127, 255, 255, 255, 143)/255
      breaks2color(breaks, pos)
    }
  )
  I <- 0:(n-1) / (n-1)
  I <- rev(I)   # make red the brightest, blue darkest.
  sapply(I, function(x)
    colorbar2rgb(x, colorbar=COLORBAR.MAT,
      shade=shade, tint=tint, min.r=min.r, max.r=max.r, 
      min.g=min.g, max.g=max.g, min.b=min.b, max.b=max.b)
  )
}

yahoo.weather.colors <- function(n, shade=0, tint=0,
  min.r=0, max.r=1, min.g=0, max.g=1, min.b=0, max.b=1) {
  COLORBAR <- list(
    RED=function(pos) {
      breaks <- c(209, 204, 255, 255, 255, 204, 84, 102, 153, 204, 255)/255
      breaks2color(breaks, pos)
    },
    GREEN=function(pos) {
      breaks <- c(73, 102, 153, 204, 255, 255, 169, 204, 255, 255, 255)/255
      breaks2color(breaks, pos)
    },
    BLUE=function(pos) {
      breaks <- c(73, 102, 102, 102, 103, 103, 255, 255, 255, 255, 255)/255
      breaks2color(breaks, pos)
    }
  )
  I <- 0:(n-1) / (n-1)
  I <- rev(I)   # make red the brightest, blue darkest.
  sapply(I, function(x)
    colorbar2rgb(x, colorbar=COLORBAR,
      shade=shade, tint=tint, min.r=min.r, max.r=max.r, 
      min.g=min.g, max.g=max.g, min.b=min.b, max.b=max.b)
  )
}

genespring.colors <- function(n, shade=0, tint=0,
  min.r=0, max.r=1, min.g=0, max.g=1, min.b=0, max.b=1) {
  COLORBAR.YWM <- list(
    RED=function(pos) {
      breaks <- c(255, 255, 0)/255
      breaks2color(breaks, pos)
    },
    GREEN=function(pos) {
      breaks <- c(0, 255, 0)/255
      breaks2color(breaks, pos)
    },
    BLUE=function(pos) {
      breaks <- c(0, 0, 255)/255
      breaks2color(breaks, pos)
    }
  )
  I <- 0:(n-1) / (n-1)
  I <- rev(I)   # make red the brightest, blue darkest.
  sapply(I, function(x)
    colorbar2rgb(x, colorbar=COLORBAR.YWM,
      shade=shade, tint=tint, min.r=min.r, max.r=max.r, 
      min.g=min.g, max.g=max.g, min.b=min.b, max.b=max.b)
  )
}

colorbar2rgb <- function(value, colorbar=COLORBAR.RGB, 
  shade=0, tint=0, min.r=0, max.r=1, min.g=0, max.g=1, min.b=0, max.b=1) {
  # value is 0 to 1, from RED around to RED again.
  # min.color, max.color, shade, tint from 0 to 1
  r <- colorbar$RED(value) + tint - shade
  g <- colorbar$GREEN(value) + tint - shade
  b <- colorbar$BLUE(value) + tint - shade
  r <- r*(max.r-min.r)+min.r
  g <- g*(max.g-min.g)+min.g
  b <- b*(max.b-min.b)+min.b
  r <- max(min(r, 1), 0)
  g <- max(min(g, 1), 0)
  b <- max(min(b, 1), 0)
  return(rgb(r, g, b, maxColorValue=1))
}

breaks2color <- function(breaks, pos) {
  # pos is [0, 1]
  # breaks[1] = pos 0.0
  # breaks[length(breaks)] = pos 1.0
  delta <- 1/(length(breaks)-1)
  i <- floor(pos/delta)+1
  if(i <= 0) return(breaks[1])
  if(i == length(breaks)) return(breaks[i])
  c1 <- breaks[i]
  c2 <- breaks[i+1]
  c1 + (pos-(i-1)*delta)/(delta/(c2-c1))
}
  
wheel2bar <- function(wheel) {
  # Only want the first 4/6 of the color wheel.  No numbers > 4/6.
  list(
    RED=function(pos) { wheel$RED(pos*4/6) },
    GREEN=function(pos) { wheel$GREEN(pos*4/6) },
    BLUE=function(pos) { wheel$BLUE(pos*4/6) }
  )
}

my.lineplot <- function(X, Y=NA, xlab=NA, ylab=NA, 
  colors=NA, SPACING=NA, lwd=1) {
  # X should be gene x samples matrix.
  # Y should be a vector of the 1-based line number for each row of X.
  if((length(Y) == 1) && is.na(Y))
    Y <- 1:nrow(X)
  num.lines <- max(Y)
  if((length(colors)==1) && is.na(colors))
    colors <- rep("#000000", nrow(X))
  if((length(ylab)==1) && is.na(ylab))
    ylab <- NA
  if((length(xlab)==1) && is.na(xlab))
    xlab <- NA
  if((length(SPACING) == 1) && is.na(SPACING))
    SPACING <- max(abs(X))

  xlim <- c(1, ncol(X))  
  ylim <- c(-SPACING, num.lines*SPACING)
  plot(NA, type="n", axes=FALSE, xlim=xlim, ylim=ylim, xlab="", ylab="")
  for(i in 1:nrow(X)) {
    offset <- (Y[i]-1) * SPACING
    x <- as.numeric(X[i,]) + offset
    lines(x, lwd=lwd, col=colors[i])
  }
  axis(1, at=1:ncol(X), labels=xlab)
  axis(2, at=seq(0, (num.lines-1)*SPACING, SPACING), labels=ylab)
  box()
}

normalize.one.mv <- function(x, M=0, V=1) {
  # Normalize a list of numbers so the mean is M and variance is V.
  V.0 <- var(x)
  M.0 <- mean(x)
  (x-M.0)*sqrt(V/V.0) + M
}

normalize.mv <- function(X, M=0, V=1) {
  t(apply(X, 1, function(x) normalize.one.mv(x, M, V)))
}

# matrix should contain values from 0 to 1.
my.heatmap <- function(matrix, col=rg.array.colors, xlab="", ylab="", 
  normalize=FALSE, scale=FALSE, cluster=FALSE) {
  # If normalize is TRUE, then will normalize to N(0, 1).  It can also
  # be a vector of (mean, variance) to specify the normalization
  # parameters.
  # If scale is TRUE, then will everything onto a 0 to 1 scale.  It
  # can also be a vector of (min, max, median) that indicates the
  # minimum, maximum, and median values that should correspond to 0,
  # 1, and 0.5.  If median is NA, then will not attempt to center the
  # median.
  # If cluster is TRUE, will cluster the rows and columns.
  if((length(normalize) == 1 && normalize == TRUE) || 
    (length(normalize) == 2)) {
    M <- 0; V <- 1
    if(length(normalize) == 2) {
      M <- normalize[1]; V <- normalize[2] 
    }
    m <- apply(matrix, 1, mean) - M
    matrix <- sweep(matrix, 1, m)
    matrix <- t(apply(matrix, 1, function(x) {
      V.0 <- var(x); M.0 <- mean(x); (x-M.0)*sqrt(V/V.0) + M.0 }))
  }

  if((length(scale) == 1 && scale == TRUE) || (length(scale) == 3)) {
    scale.min <- NA; scale.max <- NA; scale.med <- NA
    if(length(scale) != 1) {
      scale.min <- scale[1]; scale.max <- scale[2]; scale.med <- scale[3]
    }
    if(is.na(scale.min)) scale.min <- min(matrix)
    if(is.na(scale.max)) scale.max <- max(matrix)
    if(scale.max <= scale.min) stop("invalid scale parameters")
    matrix[matrix < scale.min] <- scale.min
    matrix[matrix > scale.max] <- scale.max
    if(!is.na(scale.med)) {
      # Center around 0, then scale so it's from -0.5 to +0.5.
      matrix <- matrix - scale.med
      x <- max(abs(c(min(matrix), max(matrix))))
      matrix <- matrix / x / 2 + 0.5
    } else {
      matrix <- matrix - min(matrix)
      matrix <- matrix / max(matrix)
    }
  }

  # col should be dark to bright.
  if(mode(col) == "character") {
    num.colors <- length(col)
  } else if(mode(col) == "list") {
    num.colors <- 256
    col <- rg.array.colors(num.colors)
  } else {
    num.colors <- 256
    col <- col(num.colors)
  }
  x.size <- 1
  y.size <- 1
  matrix <- as.matrix(matrix)

  # image treats the rows as x and columns as y.  So I need to
  # "rotate" the matrix 90 degrees clockwise.
  matrix <- matrix(t(matrix)[,nrow(matrix):1], ncol(matrix), nrow(matrix))

  # Weird.  If X11 is not running, then have to do it this way.
  # image puts row 1 on the bottom and column 1 on the right.  Flip
  # both the rows and columns.
  #matrix <- matrix[nrow(matrix):1,ncol(matrix):1]

  if(cluster) {
    h.r <- hclust(dist(matrix))
    h.c <- hclust(dist(t(matrix)))
    matrix <- matrix[h.r$order, h.c$order]
  }

  x <- x.size * 1:nrow(matrix)
  y <- y.size * 1:ncol(matrix)
  breaks <- (0:num.colors)/num.colors
  image(x, y, matrix, xlab=xlab, ylab=ylab, axes=FALSE, col=col, breaks=breaks)
  matrix
}

my.colorbar <- function(n=65, col=rg.array.colors) {
  # By default, Matlab plots n=65.
  I <- (n-1):0 / (n-1)
  M <- matrix(I, length(I), 1)
  my.heatmap(M, col=col)
}

sample.colorbars <- function() {
  M <- t(matrix(1:8, 4, 2))
  layout(M)
  col.fns <- c(green.shade, red.shade, rg.array.colors, by.array.colors,
    rgb.colors, ryb.colors, matlab.colors, genespring.colors)
  col.names <- c("Green", "Red", "Red/Green", "Yellow/Blue",
    "RGB", "RYB", "Matlab", "GeneSpring")
  # For some reason, I can't loop 1:length(col.fns).  col <-
  # col.fns[i] always gets me rg.array.colors.
  i <- 1
  for(col in col.fns) {
    my.colorbar(col=col)
    title(col.names[i])
    i <- i + 1
  }
}

