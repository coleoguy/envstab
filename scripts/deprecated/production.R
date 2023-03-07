#Andres Barboza

env_stab_Jan_11_results <- readRDS("env_stab_Jan_11_results.rds")
chrom.num <- c(5, 10, 15, 20, 25)
prob.change <- c(.01, .1, .3)

#### additive model ####

production.additive <- vector(mode="list", length=length(chrom.num))
names(production.additive) <- paste("chrom.num", chrom.num)
for(j in 1:length(chrom.num)){
  production.additive[[j]] <- vector(mode="list", length=length(prob.change))
  names(production.additive[[j]]) <- paste("prob.change",prob.change)
  for (k in 1:length(prob.change)) {
    fit <- c()
    for (m in 1:length(env_stab_Jan_11_results[[j]][[k]]$`model additive`)) {
      x <- env_stab_Jan_11_results[[j]][[k]]$`model additive`[[m]]$mean.fitness
      mean.last.100 <- mean(x[c(150:250)])
      fit[m] <- mean.last.100
    }
    production.additive[[j]][[k]] <- append(production.additive[[j]][[k]], mean(fit))
  }
}

production.additive <- as.data.frame(do.call(rbind, production.additive))
production.additive <- data.matrix(production.additive)
heatmap(production.additive)


#### epi.sign model ####

production.epi.sign <- vector(mode="list", length=length(chrom.num))
names(production.epi.sign) <- paste("chrom.num", chrom.num)
for(j in 1:length(chrom.num)){
  production.epi.sign[[j]] <- vector(mode="list", length=length(prob.change))
  names(production.epi.sign[[j]]) <- paste("prob.change",prob.change)
  for (k in 1:length(prob.change)) {
    fit <- list()
    for (m in 1:length(env_stab_Jan_11_results[[j]][[k]]$`model additive`)) {
      x <- env_stab_Jan_11_results[[j]][[k]]$`model epi.sign`[[k]]$mean.fitness
      mean.last.100 <- mean(x[c(150:250)])
      fit[m] <- mean.last.100
    }
    v <- unlist(fit)
    v <- as.vector(v,'numeric')
    production.epi.sign[[j]][[k]] <- append(production.epi.sign[[j]][[k]], mean(v))
  }
}

production.epi.sign <- as.data.frame(do.call(rbind, production.epi.sign))
production.epi.sign <- data.matrix(production.epi.sign)
heatmap(production.epi.sign)


#### epi.inc model ####

production.epi.inc <- vector(mode="list", length=length(chrom.num))
names(production.epi.inc) <- paste("chrom.num", chrom.num)
for(j in 1:length(chrom.num)){
  production.epi.inc[[j]] <- vector(mode="list", length=length(prob.change))
  names(production.epi.inc[[j]]) <- paste("prob.change",prob.change)
  for (k in 1:length(prob.change)) {
    fit <- list()
    for (m in 1:length(env_stab_Jan_11_results[[j]][[k]]$`model epi.inc`)) {
      x <- env_stab_Jan_11_results[[j]][[k]]$`model epi.inc`[[k]]$mean.fitness
      mean.last.100 <- mean(x[c(150:250)])
      fit[m] <- mean.last.100
    }
    v <- unlist(fit)
    v <- as.vector(v,'numeric')
    production.epi.inc[[j]][[k]] <- append(production.epi.inc[[j]][[k]], mean(v))
  }
}

production.epi.inc <- as.data.frame(do.call(rbind, production.epi.inc))
production.epi.inc <- data.matrix(production.epi.inc)
heatmap(production.epi.inc)


#### epi.dec model ####

production.epi.dec <- vector(mode="list", length=length(chrom.num))
names(production.epi.dec) <- paste("chrom.num", chrom.num)
for(j in 1:length(chrom.num)){
  production.epi.dec[[j]] <- vector(mode="list", length=length(prob.change))
  names(production.epi.dec[[j]]) <- paste("prob.change",prob.change)
  for (k in 1:length(prob.change)) {
    fit <- list()
    for (m in 1:length(env_stab_Jan_11_results[[j]][[k]]$`model epi.dec`)) {
      x <- env_stab_Jan_11_results[[j]][[k]]$`model epi.dec`[[k]]$mean.fitness
      mean.last.100 <- mean(x[c(150:250)])
      fit[m] <- mean.last.100
    }
    v <- unlist(fit)
    v <- as.vector(v,'numeric')
    production.epi.dec[[j]][[k]] <- append(production.epi.dec[[j]][[k]], mean(v))
  }
}

production.epi.dec <- as.data.frame(do.call(rbind, production.epi.dec))
production.epi.dec <- data.matrix(production.epi.dec)
heatmap(production.epi.dec)


