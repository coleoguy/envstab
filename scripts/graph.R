#Andres Barboza

run <- readRDS("../data/res-mar15.rds")
library(ggplot2)
library(reshape)
library(plyr)


#### Dr. Blackmon ####

## Outputs matrix for heatmap, chrom# vs prob of change, for a selected stat

GetRes1 <- function(runs, model, stat, gens){
  getMeans <- function(x, stat, gens){
    vals <- c()
    for(i in 1:length(x)){
      vals[i] <- mean(x[[i]][[which(names(x[[i]]) == stat)]][gens[1]:gens[2]])
    }
    mean(vals)
  }
  res <- matrix(NA, nrow = length(runs[[1]]), ncol = length(runs))
  colnames(res) <- names(runs)
  rownames(res) <- names(runs[[1]])
  for(i in 1:ncol(res)){
    for(j in 1:nrow(res)){
      cur.res <- runs[[i]][[j]][[which(names(runs[[i]][[j]]) == model)]]
      res[j, i] <- getMeans(cur.res, stat, gens)
    }
  }
  return(res)
}


res.add <- GetRes(run,
                   model = "model additive",
                   stat = "mean.fitness",
                   gens = c(150,250))


res.epis <- GetRes(run,
                   model = "model epi.sign",
                   stat = "mean.fitness",
                   gens = c(150,250))

res.epii <- GetRes(run,
                   model = "model epi.inc",
                   stat = "mean.fitness",
                   gens = c(150,250))

res.epid <- GetRes(run,
                   model = "model epi.dec",
                   stat = "mean.fitness",
                   gens = c(150,250))



####Andres####

# Outputs table, with chrom# (col#1) the stat mean (col#2) for each iteration (each row)

GetRes2 <- function(runs, probc, model, stat, gens){
  value <- c()
  chromn <- c()
  tags <- c(names(runs))
  for(i in 1:length(runs)){
    for(j in 1:length(which(names(runs[[i]]) == probc))){
      cur.res <- runs[[i]][[which(names(runs[[i]]) == probc)]][[which(names(runs[[i]][[j]]) == model)]]
      vals <- c()
      for(k in 1:length(cur.res)){
        vals[k] <- mean(cur.res[[k]][[which(names(cur.res[[k]]) == stat)]][gens[1]:gens[2]])
        chromn <- append(chromn, tags[[i]])
      }
      value <- append(value, vals)
    }
  }
  res <- data.frame(chromn, value)
  return(res)
}

add0.001 <- GetRes2(run, probc = "prob.change 0.001",
                 model = "model additive",
                 stat = "mean.fitness",
                 gens = c(1,250))

add0.061 <- GetRes2(run, probc = "prob.change 0.0608",
                 model = "model additive",
                 stat = "mean.fitness",
                 gens = c(1,250))

add0.121 <- GetRes2(run, probc = "prob.change 0.1206",
                 model = "model additive",
                 stat = "mean.fitness",
                 gens = c(1,250))

add0.180 <- GetRes2(run, probc = "prob.change 0.1804",
                 model = "model additive",
                 stat = "mean.fitness",
                 gens = c(1,250))

add0.240 <- GetRes2(run, probc = "prob.change 0.2402",
                 model = "model additive",
                 stat = "mean.fitness",
                 gens = c(1,250))

add0.300 <- GetRes2(run, probc = "prob.change 0.3",
                 model = "model additive",
                 stat = "mean.fitness",
                 gens = c(1,250))



#### previous ####
data.add <- GetRes2 (run, model = "model additive",
                          stat = "mean.fitness",
                          gens = c(150,250))

data.epis <- GetRes2 (run, model = "model epi.sign",
                     stat = "mean.fitness",
                     gens = c(150,250))

data.epii <- GetRes2 (run, model = "model epi.inc",
                     stat = "mean.fitness",
                     gens = c(150,250))

data.epid <- GetRes2 (run, model = "model epi.dec",
                     stat = "mean.fitness",
                     gens = c(150,250))
#### prev####


#### ggplot2 ####

PlotLine <- function(chrom.index, change.index, model.index, iter.index) {
  par(mar=c(8, 5, 2, 2), xpd=TRUE)
  plot(run[[chrom.index]][[change.index]][[model.index]][[iter.index]][[2]], type='l', las=1, lwd=2 ,col = '#6464FF', ylim=c(0,20), xlim=c(1, 250), ylab='Phenotype', xaxt="n")
  lines(run[[chrom.index]][[change.index]][[model.index]][[iter.index]][[3]], type='l', lwd=1.5, col = 'black')
  axis(1, tcl=0.4)

  plot(run[[chrom.index]][[change.index]][[model.index]][[iter.index]][[1]], type='l', las=1, lwd=2, col = '#C03830', ylim=c(0,1), xlim=c(1, 250), ylab='Fitness', xlab='Generations')
  legend("bottomleft", inset=c(0, -0.6), legend=c("Optimal Phenotype", "Mean Phenotype", "Mean Fitness"), col=c('black', '#6464FF', '#C03830'), lwd=2, lty=1, cex=0.8)
  axis(1, tcl=0.4)
}


PlotHm <- function(run) {
  title <- deparse(substitute(run))
  data <- melt(run)
  colnames(data) <- c("x", "y", "value")
  mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
  ggplot(data, aes(x = x, y = y, fill = value)) + 
    geom_tile(color = "white",
              linewidth = 1.5,
              linetype = 1) +
    scale_fill_gradientn(colours = mycol) +
    geom_text(aes(label = value), color = "white", size = 4) +
    coord_fixed(TRUE) +
    ggtitle(title)
}


PlotDen <- function(run) {
  title <- deparse(substitute(run))
  means <- ddply(run, "chromn", summarise, grp.mean=mean(value))
  ggplot(run, aes(x=value, color=chromn)) +
    geom_density(alpha=0.4) +
    scale_color_brewer(palette="Dark2") +
    geom_vline(data=means, aes(xintercept=grp.mean, color=chromn),
               linetype="dashed") +
    xlim(0,1) +
    ylim(0,55) +
    labs(title = title)
}


#### Plotting the models

pheno <- abs(c(0:20))
add <- pheno
inc <- (20*(pheno/20)^2)
dec <- (13*log(pheno+1)/2)
par(mar=c(7, 5, 2, 2), xpd=TRUE)

plot(c(0:20), add,  type = 'l', lwd=1.5, col='black', las=1, xlim=c(0,20), ylab='Phenotype', xlab='Contributing Alleles')
lines(c(0:20), inc, type = 'l', lwd=1.5, col='#6464FF')
lines(c(0:20), dec, type = 'l', lwd=1.5, col='#C03830')





plot(c(0:20), (1-(abs((20-add)/20))),  type = 'l', lwd=1.5, col='black', las=1, xlim=c(0,20), ylab='Fitness', xlab='Contributing Alleles')
lines(c(0:20), (1-(abs((20-inc)/20))), type = 'l', lwd=1.5, col='#6464FF')
lines(c(0:20), (1-(abs((20-dec)/20))), type = 'l', lwd=1.5, col='#C03830')
legend("bottomleft", inset=c(0, -0.4), legend=c("Additive", "Increasing Epistasis", "Decreasing Epistasis"), col=c('black', '#6464FF', '#C03830'), lwd=2, lty=1, cex=0.8)


# if(model == "additive"){
#   pheno <- sum(genome[,imp.loci])
#   w <- 1 - (abs(fav.pheno - pheno)/20)
# }
# 
# if (model== "epi.inc"){
#   x <- sum(genome[,imp.loci])
#   pheno <- (20*(x/20)^2)
# }
# 
# if(model== "epi.dec"){
#   x <- sum(genome[,imp.loci])
#   pheno <- (13*log(x+1)/2)
# }
# w <- 1 - (abs(fav.pheno - pheno)/20)
# return(w)

foo <- pheno+1

























