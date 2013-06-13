#!/usr/bin/env Rscript

## edens - electrophoretic densitometry
## Copyright (C) 2008 Paul Ryvkin
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## You should have received a copy of the GNU General Public License along
## with this program; if not, write to the Free Software Foundation, Inc.,
## 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.



require('EBImage')
require('methods')


nargs = length(commandArgs(T))
if (nargs < 1) {
  cat("USAGE: Rscript densitometry.R conf_file\n")
  q()
}

source(commandArgs(T)[1])
 
# median of each row in an image
profileimgy = function(img) {
  apply(img,2,median)
}

localmaxima = function(y, width=2, refinewidth=4) {
  # smooth it, find smoothed maxima, refine to get exact maxima
  s = smooth.spline(y, spar=0.5)
  dy = diff(s$y, width)/width
  dy = c(rep(dy[1],width/2), dy, rep(dy[length(dy)],width/2))
  dyt1 = c(dy[2:length(dy)], dy[length(dy)])
  smoothedmaxima = which( (dy > 0 & dyt1 < 0) | (dy == 0) )
  refinedmaxima = c()

  for ( i in smoothedmaxima ) {
    refinewindow = max(c(0,
      (i-refinewidth/2))):min(c((i+refinewidth/2),length(y)))
    possiblerefinedmaxima = which( y == max(y[refinewindow]) )
    refinedmaximum = floor(mean(
      possiblerefinedmaxima[possiblerefinedmaxima %in% refinewindow]))
    refinedmaxima = c(refinedmaxima, refinedmaximum)
  }
  refinedmaxima
}

img2matrix = function(img) {
    matrix(imageData(img), dim(img)[1], dim(img)[2])
}

###

pdf(paste(outprefix,'_','plots.pdf',sep=''))

orig = readImage(inputfn)

image(img2matrix(orig), main='original image', col=gray((0:256)/256))

lanesamplew = (lanesamplex[length(lanesamplex)] - lanesamplex[1])+1
lanesampleh = (lanesampley[length(lanesampley)] - lanesampley[1])+1

lanes = Image(0, c(lanesamplew, lanesampleh, numlanes))

laneprofiles = array(dim=c(numlanes, lanesampleh))

for(i in 1:numlanes) {
  lanes[,,i] = imageData(orig)[firstlanex + (i-1)*lanesep + lanesamplex,
         lanesampley, 1 ]

  laneprofiles[i,] = profileimgy(lanes[,,i])
}


opar=par()

layout(matrix(1:numlanes,numlanes,1,byrow=T))
par(mar=c(0,0,0,0))
for(i in 1:numlanes) {
  plot(laneprofiles[i,],type='l', xaxt='n',ylim=c(0,1))
}

par(opar)

### calibrate bp
# only want the top N corresponding to theoretical ladder peaks
peaks = localmaxima(laneprofiles[ladder,], 2, 12)
peaks = peaks[order(laneprofiles[ladder,peaks], decreasing=T)][1:length(ladderbp)]
if (exists("ladder.manual")) {
    peaks = ladder.manual
}

plot(laneprofiles[ladder,], t='l', main='ladder', xlab='pixel', ylab='ladder intensity')

s = smooth.spline(laneprofiles[ladder,],spar=0.5)
lines(s$y,col='red')

abline(v = peaks, col='blue')


rladderbp = rev(ladderbp)
ladderpixels = peaks[order(peaks)]

start = c(-200, -1/200, 10)
names(start) = c('p1','p2','p3')
weights = (1:length(rladderbp))**3  # weigh smaller bp more heavily
                                       # in regression
                                       # => more accuracy at smaller bp
fit = nls(rladderbp ~ p1 + exp(p2 * ladderpixels + p3), start=start, weights=weights )
bp = coef(fit)[1] + exp(coef(fit)[2] * (1:lanesampleh) + coef(fit)[3])

plot(ladderpixels, rladderbp, main='ladder calibration',xlab='pixel', ylab='bp')
lines(1:lanesampleh, bp)


# model background

bgroundx = c()
bgroundy = c()
bgroundz = c()

xsampleskip = 2
ysampleskip = 2

np = 1

for (laneno in bgroundlanes) {
  for(x in (firstlanex+(laneno-1)*lanesep + lanesamplex)[seq(1,length(lanesamplex),xsampleskip)] ) {
    for(y in lanesampley[seq(1,length(lanesampley), ysampleskip)] ) {
      bgroundx[np] = x
      bgroundy[np] = y
      np = np + 1
    }
  }
}

bgroundz = orig[ (bgroundy-1)*(dim(orig)[1]) + bgroundx]

bground = Image(0, dim=dim(orig))

bground[ (bgroundy-1)*(dim(orig)[1]) + bgroundx] = bgroundz


fit = nls(bgroundz ~ pAx * bgroundx**2 + pAy * bgroundy**2 + pBx * bgroundx + pBy * bgroundy + pC)

bgroundmodelled = Image(0, dim=dim(orig))

mdlxr = min(bgroundx):max(bgroundx)
mdlyr = min(bgroundy):max(bgroundy)
#mdlxr = 1:(dim(orig)[1])
#mdlyr = 1:(dim(orig)[2])

mdlx = rep(0, length(mdlxr)*length(mdlyr))
mdly = rep(0, length(mdlxr)*length(mdlyr))
mdlz = rep(0, length(mdlxr)*length(mdlyr))
np = 1

for(x in mdlxr) {
  for(y in mdlyr) {
    mdlx[np] = x #- min(bgroundx)
    mdly[np] = y #- min(bgroundy)
    np = np+1
  }
}

mdlz = coef(fit)[1] * mdlx**2 + coef(fit)[2] * mdly**2 +
  coef(fit)[3] * mdlx + coef(fit)[4] * mdly + coef(fit)[5]

#bgroundmodelled[(min(bgroundy)+mdly-1)*(dim(orig)[1]) +
#                min(bgroundx)+mdlx] = mdlz
bgroundmodelled[(mdly-1)*(dim(orig)[1]) + mdlx] = mdlz

origminusbground = orig - bgroundmodelled

image(img2matrix(bgroundmodelled), main='modelled bground',col=gray((0:256)/256))

image(img2matrix(origminusbground), main='modelled bground subtracted',
      col=gray((0:256)/256))


# correct specified lanes

for(i in correctlanes) {
  lanes[,,i] = imageData(origminusbground)[firstlanex + (i-1)*lanesep + lanesamplex,
         lanesampley, 1]

  laneprofiles[i,] = profileimgy(lanes[,,i])
}


# display them again

opar=par()

layout(matrix(1:numlanes,numlanes,1,byrow=T))
par(mar=c(0,0,1,0))
for(i in 1:numlanes) {

  col = 'black'

  if (i %in% correctlanes ) {
    col = 'red'
  }
  
  plot(laneprofiles[i,],type='l', xaxt='n',ylim=c(0,1), col=col)

  if (i==1) {
    title(paste('lanes ',min(correctlanes),'-',max(correctlanes), ' are bground corrected', sep=''))
  } else {
    par(mar=c(0,0,0,0))
  }

}

par(opar)


# now display the bp-calibrated distributions

minpix = 1

if (sum(bp >= maxbp) > 0) {
  minpix = max(which(bp >= maxbp))
}

maxpix = length(lanesampley)

# normalize to densities when x is scaled to [0,1]
lanedensities = array(0,dim=c(numlanes,length(minpix:maxpix)))

x = rev(bp[minpix:maxpix])
dx = diff(x)
dx = c(dx,dx[length(dx)])

for(i in datalanes) {
  y = rev(laneprofiles[i, minpix:maxpix])
  y[y < 0] = 0
  lanedensities[i,] = y / sum(y*abs(dx))
}

opar = par()

nplotrow = sqrt(length(datalanes))
nplotcol = nplotrow
if (nplotrow != round(nplotrow)) {
  nplotcol = floor(nplotrow)
  nplotrow = ceiling(nplotrow)
}

layout(matrix(1:numlanes,nplotrow,nplotcol,byrow=T))

for(i in 1:length(datalanes)) {
  maxdens = max(lanedensities[i,])

  laneno = datalanes[i]
  lanetitle = datalanetitles[i]
  y = lanedensities[laneno,]

  s = smooth.spline(x,y, spar=0.7)
  
#  plot(c(),xlim=c(0,maxbp), ylim=c(0,maxdens),
  plot(s,
       xlab='fragment size (bp)', ylab='density', xaxt='n',
       main=datalanetitles[i], t='l')
  axis(side=1, seq(0,max(x),5))
  par(tcl=-0.25)
  axis(side=1, seq(0,max(x),5),labels=F)

  
  cdf = cumsum(y*abs(dx))
  #quantiles = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  quantiles = c(0.1, 0.25, 0.5, 0.75, 0.9)
  for (q in quantiles) {
    show(c(q,x[which(cdf >= q)[1]]))
  }
  
  xatmax = x[y == max(y)]
  maxy = y[y == max(y)]
  
}

dev.off()
