#Andres Barboza

setwd("~/Documents/env-stab")
run <- readRDS("res-feb-16.rds")
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

data1 <- GetRes2(run, probc = "prob.change 0.001",
                 model = "model additive",
                 stat = "mean.fitness",
                 gens = c(150,250))

data2 <- GetRes2(run, probc = "prob.change 0.0608",
                 model = "model additive",
                 stat = "mean.fitness",
                 gens = c(150,250))

data3 <- GetRes2(run, probc = "prob.change 0.1206",
                 model = "model additive",
                 stat = "mean.fitness",
                 gens = c(150,250))

data4 <- GetRes2(run, probc = "prob.change 0.1804",
                 model = "model additive",
                 stat = "mean.fitness",
                 gens = c(150,250))

data5 <- GetRes2(run, probc = "prob.change 0.2402",
                 model = "model additive",
                 stat = "mean.fitness",
                 gens = c(150,250))

data6 <- GetRes2(run, probc = "prob.change 0.3",
                 model = "model additive",
                 stat = "mean.fitness",
                 gens = c(150,250))



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

PlotHm <- function(run, title) {
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
  means <- ddply(run, "chromn", summarise, grp.mean=mean(value))
  ggplot(run, aes(x=value, color=chromn)) +
    geom_density(alpha=0.4) +
    scale_color_brewer(palette="Dark2") +
    geom_vline(data=means, aes(xintercept=grp.mean, color=chromn),
               linetype="dashed") +
    xlim(0,1) +
    ylim(0,30)
}

