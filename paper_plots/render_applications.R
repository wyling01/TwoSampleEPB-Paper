#!/usr/bin/env Rscript
# Render the proposed macrophage and ChIP-seq CDF/joint-prior pairs.
# Uses the saved fits, without refitting or renormalizing probability masses.

args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=', args)])
script <- gsub('~+~', ' ', script, fixed=TRUE)
script_dir <- dirname(normalizePath(script))
root <- dirname(script_dir)
data_root <- file.path(root, 'Real Dataset')
out <- file.path(script_dir, 'figures')
if(!nzchar(Sys.which('pdflatex'))) stop('Install pdflatex with TikZ, standalone, and Latin Modern.')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
out <- normalizePath(out)
required_data <- file.path(data_root, c(
  'GSE9101/Result/1D_GSE9101.csv', 'GSE9101/Result/2D_GSE9101.csv',
  'chip-seq/Result/1D_chipseq.csv', 'chip-seq/Result/2D_chipseq.csv'
))
if(any(!file.exists(required_data)))
  stop('Missing application data file(s):\n', paste(required_data[!file.exists(required_data)], collapse='\n'))
work <- tempfile('twosample-figures-')
dir.create(work)

compile <- function(stem) {
  log <- file.path(work, paste0(stem, '_compile.txt'))
  status <- system2(unname(Sys.which('pdflatex')),
    shQuote(c('-no-shell-escape', '-interaction=nonstopmode', '-halt-on-error',
              paste0('-output-directory=', work), file.path(work, paste0(stem, '.tex')))),
    stdout=log, stderr=log)
  if(status!=0L) stop('pdflatex failed:\n', paste(tail(readLines(log), 25), collapse='\n'))
}

# GM12891 / GM12892 variance ratio from fit.info$ratio.var in chip-seq/chipseq.R.
manorm_ratio <- 1.2756994798302
joint_specs <- list(
  macrophage=list(file='GSE9101/Result/2D_GSE9101.csv',
                  lim=c(-5.8,-2.8), ticks=c(-5,-4,-3),
                  zmax=.15, zticks=c(0,.05,.1,.15)),
  chipseq=list(file='chip-seq/Result/2D_chipseq.csv',
               lim=c(-6.2,.5), ticks=c(-6,-4,-2,0),
               zmax=.35, zticks=c(0,.1,.2,.3), ratio=manorm_ratio))
cdf_specs <- list(
  macrophage=list(file='GSE9101/Result/1D_GSE9101.csv',
                  lim=c(-2,2), ticks=-2:2),
  chipseq=list(file='chip-seq/Result/1D_chipseq.csv',
               lim=c(-.05,.35), ticks=c(0,.1,.2,.3), ratio=manorm_ratio))

# Orthographic camera, with equal scaling for the two log-variance axes.
project <- function(x,y,z=0,lim,zmax) {
  az <- -60*pi/180; el <- 30*pi/180
  xx <- (x-mean(lim))/diff(lim); yy <- (y-mean(lim))/diff(lim)
  zz <- .83*z/zmax
  cbind(x=135+130*(-sin(az)*xx+cos(az)*yy),
        y=95+130*(-sin(el)*cos(az)*xx-sin(el)*sin(az)*yy+cos(el)*zz),
        depth=cos(el)*cos(az)*xx+cos(el)*sin(az)*yy+sin(el)*zz)
}
fmt <- function(v) sprintf('%.5f',v)
pt <- function(x) sprintf('(%.5f,%.5f)',x[1],x[2])

render_joint <- function(name) {
  sp <- joint_specs[[name]]
  dat <- read.csv(file.path(data_root,sp$file))
  stopifnot(all(is.finite(as.matrix(dat[c('x','y','prob')]))),
            all(dat$x>0),all(dat$y>0),all(dat$prob>=0),abs(sum(dat$prob)-1)<1e-8)
  d <- data.frame(x=log(dat$x),y=log(dat$y),p=dat$prob)
  inside <- d$x>=sp$lim[1]&d$x<=sp$lim[2]&d$y>=sp$lim[1]&d$y<=sp$lim[2]
  retained <- sum(d$p[inside])
  stopifnot(retained>0,retained<=1+1e-8,max(d$p)<=sp$zmax)
  d <- d[inside & d$p>0,]
  pr <- function(x,y,z=0) project(x,y,z,sp$lim,sp$zmax)
  lo <- sp$lim[1]; hi <- sp$lim[2]
  a <- pr(lo,lo); b <- pr(hi,lo); c <- pr(hi,hi); e <- pr(lo,hi)
  tex <- c('\\documentclass[border=0pt]{standalone}',
           '\\usepackage[T1]{fontenc}', '\\usepackage{lmodern}', '\\usepackage{amsmath}',
           '\\usepackage{tikz}', '\\usetikzlibrary{shadings}',
           '\\definecolor{ink}{HTML}{27343C}', '\\definecolor{stem}{HTML}{315C73}',
           '\\definecolor{stemlight}{HTML}{7F9EAE}', '\\definecolor{stemdark}{HTML}{284E64}',
           '\\definecolor{stemtop}{HTML}{91ACBA}', '\\definecolor{axis}{HTML}{929DA5}',
           '\\definecolor{backedge}{HTML}{D2DBE1}', '\\definecolor{identity}{HTML}{7E8C96}',
           '\\definecolor{manorm}{HTML}{AC7451}',
           '\\definecolor{floor}{HTML}{F6F8FA}',
           '\\definecolor{grid}{HTML}{FFFFFF}',
           '\\begin{document}',
           '\\begin{tikzpicture}[x=1bp,y=1bp,every node/.style={inner sep=0pt,outer sep=0pt,text=ink,font=\\fontsize{8.5}{10}\\selectfont}]',
           '\\path[use as bounding box] (0,0) rectangle (241.2,205.2);',
           '\\begin{scope}', '\\clip (0,0) rectangle (241.2,205.2);')
  put <- function(s) tex <<- c(tex,s)
  line <- function(x,y,style) put(sprintf('\\draw[%s] %s -- %s;',style,pt(x),pt(y)))
  node <- function(pos,txt,opts='') put(sprintf('\\node[%s] at %s {%s};',opts,pt(pos),txt))
  put(sprintf('\\fill[floor] %s -- %s -- %s -- %s -- cycle;',pt(a),pt(b),pt(c),pt(e)))
  for(t in sp$ticks) {
    line(pr(lo,t),pr(hi,t),'grid,line width=.45bp')
    line(pr(t,lo),pr(t,hi),'grid,line width=.45bp')
  }
  line(a,e,'backedge,line width=.4bp');line(e,c,'backedge,line width=.4bp')
  line(a,b,'axis,line width=.5bp');line(b,c,'axis,line width=.5bp')
  line(a,pr(lo,lo,sp$zmax),'axis,line width=.5bp')
  line(a,c,'identity,line width=.6bp,dash pattern=on 3bp off 3bp')
  if(!is.null(sp$ratio)) {
    offset <- log(sp$ratio); xx <- c(max(lo,lo+offset),min(hi,hi+offset))
    stopifnot(all(xx>=lo & xx<=hi),all(xx-offset>=lo-1e-12 & xx-offset<=hi+1e-12))
    line(pr(xx[1],xx[1]-offset),pr(xx[2],xx[2]-offset),'manorm,line width=.75bp')
  }
  base <- pr(d$x,d$y); tip <- pr(d$x,d$y,d$p)
  # Draw back to front. Thin short stems proportionally; never inflate their height.
  for(j in order(base[,3])) {
    h <- tip[j,2]-base[j,2]
    diameter <- min(1.3,.3*h)
    if(h<=0 || diameter<=0) next
    x <- base[j,1]; y0 <- base[j,2]; y1 <- tip[j,2]
    r <- diameter/2; ry <- r*.5
    opacity <- min(1,d$p[j]/1e-4) # Fade numerical dust below 0.01% mass.
    if(h>=1) {
      put(sprintf('\\begin{scope}[transparency group,opacity=%.5f]',opacity))
      put(sprintf('\\shade[left color=stemdark,middle color=stemlight,right color=stem] %s ellipse[x radius=%sbp,y radius=%sbp];',pt(c(x,y0)),fmt(r),fmt(ry)))
      put(sprintf('\\shade[left color=stemdark,middle color=stemlight,right color=stem] %s rectangle %s;',pt(c(x-r,y0)),pt(c(x+r,y1))))
      put(sprintf('\\fill[stemtop] %s ellipse[x radius=%sbp,y radius=%sbp];',pt(c(x,y1)),fmt(r),fmt(ry)))
      put('\\end{scope}')
    } else {
      put(sprintf('\\fill[stem,opacity=%.5f] %s rectangle %s;',opacity,pt(c(x-r,y0)),pt(c(x+r,y1))))
    }
  }
  xdir <- b[1:2]-a[1:2]; ydir <- c[1:2]-b[1:2]
  xn <- c(xdir[2],-xdir[1])/sqrt(sum(xdir^2)); yn <- c(ydir[2],-ydir[1])/sqrt(sum(ydir^2))
  ticklabel <- function(t) sprintf('$%g$',t)
  for(t in sp$ticks) {
    p <- pr(t,lo)[1:2];line(p,p+2.5*xn,'axis,line width=.45bp');node(p+8*xn,ticklabel(t))
    p <- pr(hi,t)[1:2];line(p,p+2.5*yn,'axis,line width=.45bp')
    if(!(t==lo && hi %in% sp$ticks)) node(p+8*yn,ticklabel(t))
  }
  for(t in sp$zticks) {
    p <- pr(lo,lo,t)[1:2];line(p-c(2.5,0),p,'axis,line width=.45bp')
    node(p-c(5,0),sprintf('$%g\\%%$',100*t),'anchor=east')
  }
  node(pr(mean(sp$lim),lo)[1:2]+19*xn,'$\\log(\\sigma_{iA}^{2})$',
       sprintf('rotate=%.5f,font=\\fontsize{9.5}{11}\\selectfont',atan2(xdir[2],xdir[1])*180/pi))
  node(pr(hi,mean(sp$lim))[1:2]+19*yn,'$\\log(\\sigma_{iB}^{2})$',
       sprintf('rotate=%.5f,font=\\fontsize{9.5}{11}\\selectfont',atan2(ydir[2],ydir[1])*180/pi))
  node(c(11,mean(c(a[2],pr(lo,lo,sp$zmax)[2]))),'Probability mass','rotate=90')
  if(!is.null(sp$ratio)) {
    line(c(36,19),c(50,19),'identity,line width=.6bp,dash pattern=on 3bp off 3bp')
    node(c(55,19),'Equal variances','anchor=west,font=\\fontsize{8.5}{10}\\selectfont')
    line(c(126,19),c(140,19),'manorm,line width=.75bp')
    node(c(145,19),'MAnorm2 ratio','anchor=west,font=\\fontsize{8.5}{10}\\selectfont')
  } else {
    line(c(79,19),c(93,19),'identity,line width=.6bp,dash pattern=on 3bp off 3bp')
    node(c(98,19),'Equal variances','anchor=west,font=\\fontsize{8.5}{10}\\selectfont')
  }
  put('\\end{scope}'); put('\\end{tikzpicture}');put('\\end{document}')
  stem <- paste0(name,'_stems')
  writeLines(tex,file.path(work,paste0(stem,'.tex')))
  compile(stem)
  message(sprintf('%s: joint-prior window contains %.6f%% of fitted mass.',
                  name,100*retained))
}

xy <- function(x,y) sprintf('(%.5f,%.5f)', x,y)
render_cdf <- function(name) {
  sp <- cdf_specs[[name]]
  d <- read.csv(file.path(data_root,sp$file))
  d <- d[order(d$x),]
  stopifnot(all(is.finite(d$x)),all(d$x>0),all(is.finite(d$prob)),
            all(d$prob>=0),abs(sum(d$prob)-1)<1e-8)
  x <- log(d$x); cumulative <- cumsum(d$prob)
  lo <- sp$lim[1]; hi <- sp$lim[2]
  # Accumulate all atoms, including the tails outside the displayed window.
  # Each vertical segment is at the original support point, without merging.
  inside <- which(x>=lo & x<=hi)
  xx <- c(lo,rep(x[inside],each=2),hi)
  yy <- c(sum(d$prob[x<lo]),
          as.vector(rbind(cumulative[inside]-d$prob[inside],cumulative[inside])),
          sum(d$prob[x<=hi]))
  stopifnot(all(diff(xx)>=0),all(diff(yy)>=-1e-12),
            min(yy)>=0,max(yy)<=1+1e-8)
  left <- 46; right <- 225; bottom <- 64; top <- 180
  px <- function(x) left+(x-lo)/(hi-lo)*(right-left)
  py <- function(y) bottom+(y+.025)/1.05*(top-bottom)
  tex <- c('\\documentclass[border=0pt]{standalone}',
           '\\usepackage[T1]{fontenc}', '\\usepackage{lmodern}',
           '\\usepackage{amsmath}', '\\usepackage{tikz}',
           '\\definecolor{ink}{HTML}{27343C}',
           '\\definecolor{stem}{HTML}{315C73}',
           '\\definecolor{axis}{HTML}{929DA5}',
           '\\definecolor{identity}{HTML}{7E8C96}',
           '\\definecolor{manorm}{HTML}{AC7451}',
           '\\definecolor{grid}{HTML}{E8EDF0}',
           '\\begin{document}',
           '\\begin{tikzpicture}[x=1bp,y=1bp,every node/.style={inner sep=0pt,outer sep=0pt,text=ink,font=\\fontsize{8.5}{10}\\selectfont}]',
           '\\path[use as bounding box] (0,0) rectangle (241.2,205.2);')
  put <- function(s) tex <<- c(tex,s)
  line <- function(x1,y1,x2,y2,style) {
    put(sprintf('\\draw[%s] %s -- %s;',style,xy(x1,y1),xy(x2,y2)))
  }
  node <- function(x,y,txt,style='') {
    put(sprintf('\\node[%s] at %s {%s};',style,xy(x,y),txt))
  }
  for(t in seq(0,1,.25)) line(left,py(t),right,py(t),'grid,line width=.35bp')
  line(px(0),bottom,px(0),top,'identity,line width=.6bp,dash pattern=on 3bp off 3bp')
  put('\\begin{scope}')
  put(sprintf('\\clip %s rectangle %s;',xy(left,bottom),xy(right,top)))
  if(!is.null(sp$ratio)) {
    refx <- c(lo,log(sp$ratio),log(sp$ratio),hi)
    refy <- c(0,0,1,1)
    put(sprintf('\\draw[manorm,line width=.8bp,line join=round] %s;',
                paste(xy(px(refx),py(refy)),collapse=' -- ')))
  }
  put(sprintf('\\draw[stem,line width=1bp,line join=round] %s;',
              paste(xy(px(xx),py(yy)),collapse=' -- ')))
  put('\\end{scope}')
  line(left,bottom,right,bottom,'axis,line width=.5bp')
  line(left,bottom,left,top,'axis,line width=.5bp')
  for(t in sp$ticks) {
    line(px(t),bottom,px(t),bottom-2.5,'axis,line width=.45bp')
    node(px(t),bottom-9,sprintf('$%g$',t))
  }
  for(t in seq(0,1,.25)) {
    line(left-2.5,py(t),left,py(t),'axis,line width=.45bp')
    node(left-5,py(t),sprintf('$%g$',t),'anchor=east')
  }
  node(11,mean(c(bottom,top)),'CDF','rotate=90')
  node(mean(c(left,right)),36,'$\\log(\\lambda_i)$',
       'font=\\fontsize{9.5}{11}\\selectfont')
  if(is.null(sp$ratio)) {
    line(91,19,105,19,'stem,line width=1bp')
    node(110,19,'VREPB','anchor=west')
  } else {
    line(46,19,60,19,'stem,line width=1bp')
    node(65,19,'VREPB','anchor=west')
    line(127,19,141,19,'manorm,line width=.8bp')
    node(146,19,'MAnorm2 ratio','anchor=west')
  }
  put('\\end{tikzpicture}');put('\\end{document}')
  writeLines(tex,file.path(work,paste0(name,'_cdf.tex')))
  compile(paste0(name,'_cdf'))
}

render_pair <- function(name) {
  # Actual manuscript width: 5.7 inches, with two 2.8-inch panels.
  # Embedded PDFs remain entirely vector. No titles or subtitles are added.
  panel_width <- 201.6; gutter <- 7.2
  height <- 205.2*(2.8/3.35)
  tex <- c('\\documentclass[border=0pt]{standalone}',
           '\\usepackage[T1]{fontenc}', '\\usepackage{lmodern}',
           '\\usepackage{graphicx}', '\\usepackage{tikz}',
           '\\definecolor{ink}{HTML}{27343C}',
           '\\begin{document}',
           '\\begin{tikzpicture}[x=1bp,y=1bp,every node/.style={inner sep=0pt,outer sep=0pt}]',
           sprintf('\\path[use as bounding box] (0,0) rectangle (410.4,%.5f);',height),
           sprintf('\\node[anchor=south west] at (0,0) {\\includegraphics[width=2.8in]{\\detokenize{%s}}};',
                   file.path(work,paste0(name,'_cdf.pdf'))),
           sprintf('\\node[anchor=south west] at (208.8,0) {\\includegraphics[width=2.8in]{\\detokenize{%s}}};',
                   file.path(work,paste0(name,'_stems.pdf'))))
  for(j in 0:1) tex <- c(tex,
    sprintf('\\node[anchor=north west,text=ink,font=\\fontsize{8}{10}\\selectfont] at (%.5f,%.5f) {(%s)};',
            j*(panel_width+gutter)+7,height-5,c('a','b')[j+1]))
  tex <- c(tex,'\\end{tikzpicture}','\\end{document}')
  writeLines(tex,file.path(work,paste0(name,'_pair.tex')))
  compile(paste0(name,'_pair'))
  filename <- paste0(name,'_pair.pdf')
  if(!file.copy(file.path(work,filename),file.path(out,filename),overwrite=TRUE))
    stop('Could not write ',filename)
}

for(name in names(joint_specs)) {
  render_joint(name)
  render_cdf(name)
  render_pair(name)
}
message('Wrote macrophage_pair.pdf and chipseq_pair.pdf in ',out)

