# Originally J. Strope now Andres and Heath

source("functions.R")
# parameters for testing
loci <- 100
popsize <- 100  #1000
chrom.num <- c(5, 10, 15, 20, 25)
generations <- 20   #250
#probbility that the favored phenotype for the environment will change 
prob.change <- c(.01, .1, .3)
#loci that will impact the pheonotype
imp.loci <- sample(1:loci, 10)
#what phenotype is favored in the environment
fav.pheno <- runif(1, min=0, max=20)
iterations <- 4 #100
mut.rate <- .01
model <- c("additive"  , "epi.sign", "epi.inc","epi.dec")

#creating a vector to store the results into and having vectors within it for 
#the prob change, model, and chrom num
results <- vector(mode="list", length=length(chrom.num))
names(results) <- paste("chrom.num", chrom.num)
for(j in 1:length(chrom.num)){
  results[[j]] <- vector(mode="list", length=length(prob.change))
  names(results[[j]]) <- paste("prob.change",prob.change)
  for(k in 1:length(prob.change)){
    results[[j]][[k]] <- vector(mode="list", length=length(model))
    names(results[[j]][[k]]) <- paste("model",model)
  }
}

#iterating through and storing the results for each chrom num, prob change, model 
#in "results" vector 
cur.result <- list()
for(i in 1:length(chrom.num)){
  for(j in 1:length(prob.change)){
    for(k in 1:length(model)){
      for(m in 1:iterations){
        print(paste("running generation", m))
        cur.result[[m]] <- simulate(loci, chrom.num[i], popsize, generations, prob.change[j], 
                                  mut.rate, model[k])
      }
      results[[i]][[j]][[k]] <- cur.result
    }
  }
}

foo <- simulate(loci=100, chrom.num =5, popsize=500, generations=200,
               prob.change=.05, mut.rate=.01, model="epi.sign")

# load("run2")

#additive########

#0.016
add01.mse <- as.data.frame(matrix(NA,250,5))
colnames(add01.mse)<-c(5,10,15,20,25)

#iterate 1:5 to fill in the specified results for each chromosome number
for(j in 1:5){
  cur.res <- results[[j]]$`prob.change 0.01`$`model additive`
  #setting up a matrix to hold the squared error for 500 generations and 1000 
  #iterations
  temptable <- matrix(NA,250,1000)
  for(k in 250:500){
    for(i in 1:1000){
      #filling in the squared error over 1000 iterations
      temptable[,i] <- cur.res[[i]]$sq.err[[k]]
    }
    #taking the mean of the squared error
    add01.mse[,j] <- rowMeans(temptable)
  }
}

#0.001
add001.mse <- as.data.frame(matrix(NA,250,5))
colnames(add001.mse)<-c(5,10,15,20,25)

for(j in 1:5){
  cur.res <- results[[j]]$`prob.change 0.001`$`model additive`
  temptable <- matrix(NA,250,1000)
  for(k in 250:500){
    for(i in 1:1000){
      temptable[,i] <- cur.res[[i]]$sq.err[[k]]
    }
    add001.mse[,j] <- rowMeans(temptable)
  } 
}


#0.1
add1.mse <- as.data.frame(matrix(NA,250,5))
colnames(add1.mse)<-c(5,10,15,20,25)

for(j in 1:5){
  cur.res <- results[[j]]$`prob.change 0.1`$`model additive`
  temptable <- matrix(NA,250,1000)
  for(k in 250:500){
    for(i in 1:1000){
      temptable[,i] <- cur.res[[i]]$sq.err[[k]]
    }
    add1.mse[,j] <- rowMeans(temptable)
  }
}

#make a table to compare#
add.001.mean <- c()
for(i in 1:5){
  add.001.mean[i] <- mean(add001.mse[,i])
}
add.01.mean <- c()
for(i in 1:5){
  add.01.mean[i] <- mean(add01.mse[,i])
}
add.1.mean <- c()
for(i in 1:5){
  add.1.mean[i] <- mean(add1.mse[,i])
}

comp.add.mse <- matrix(c(add.001.mean,add.01.mean,add.1.mean), byrow=FALSE,ncol=3)
colnames(comp.add.mse)<-c(0.001, 0.01, 0.1)
rownames(comp.add.mse)<-c(5,10,15,20,25)

###epi.inc#####

#0.01
inc01.mse <- as.data.frame(matrix(NA,250,5))
colnames(add01.mse)<-c(5,10,15,20,25)

for(j in 1:5){
  cur.res <- results[[j]]$`prob.change 0.01`$`model epi.inc`
  temptable <- matrix(NA,250,1000)
  for(k in 250:500){
    for(i in 1:1000){
      temptable[,i] <- cur.res[[i]]$sq.err[[k]]
    }
    inc01.mse[,j] <- rowMeans(temptable)
  }
}

#0.001
inc001.mse <- as.data.frame(matrix(NA,250,5))
colnames(inc001.mse)<-c(5,10,15,20,25)

for(j in 1:5){
  cur.res <- results[[j]]$`prob.change 0.001`$`model epi.inc`
  temptable <- matrix(NA,250,1000)
  for(k in 250:500){
    for(i in 1:1000){
      temptable[,i] <- cur.res[[i]]$sq.err[[k]]
    }
    inc001.mse[,j] <- rowMeans(temptable)
  }
}

#0.1
inc1.mse <- as.data.frame(matrix(NA,250,5))
colnames(inc1.mse)<-c(5,10,15,20,25)

for(j in 1:5){
  cur.res <- results[[j]]$`prob.change 0.1`$`model epi.inc`
  temptable <- matrix(NA,250,1000)
  for(k in 250:500){
    for(i in 1:1000){
      temptable[,i] <- cur.res[[i]]$sq.err[[k]]
    }
    inc1.mse[,j] <- rowMeans(temptable)
  } 
}


#make a table to compare#
inc.001.mean <- c()
for(i in 1:5){
  inc.001.mean[i] <- mean(inc001.mse[,i])
}
inc.01.mean <- c()
for(i in 1:5){
  inc.01.mean[i] <- mean(inc01.mse[,i])
}
inc.1.mean <- c()
for(i in 1:5){
  inc.1.mean[i] <- mean(inc1.mse[,i])
}

comp.inc.mse <- matrix(c(inc.001.mean,inc.01.mean,inc.1.mean), byrow=FALSE,ncol=3)
colnames(comp.inc.mse)<-c(0.001, 0.01, 0.1)
rownames(comp.inc.mse)<-c(5,10,15,20,25)



######epi.dec####

#0.01
dec01.mse <- as.data.frame(matrix(NA,250,5))
colnames(dec01.mse)<-c(5,10,15,20,25)

for(j in 1:5){
  cur.res <- results[[j]]$`prob.change 0.01`$`model epi.dec`
  temptable <- matrix(NA,250,1000)
  for(k in 250:500){
    for(i in 1:1000){
      temptable[,i] <- cur.res[[i]]$sq.err[[k]]
    }
    dec01.mse[,j] <- rowMeans(temptable)
  }
}

#0.001
dec001.mse <- as.data.frame(matrix(NA,250,5))
colnames(dec001.mse)<-c(5,10,15,20,25)

for(j in 1:5){
  cur.res <- results[[j]]$`prob.change 0.001`$`model epi.dec`
  temptable <- matrix(NA,250,1000)
  for(k in 250:500){
    for(i in 1:1000){
      temptable[,i] <- cur.res[[i]]$sq.err[[k]]
    }
    dec001.mse[,j] <- rowMeans(temptable)
  }
}

#0.1
dec1.mse <- as.data.frame(matrix(NA,250,5))
colnames(dec1.mse)<-c(5,10,15,20,25)

for(j in 1:5){
  cur.res <- results[[j]]$`prob.change 0.1`$`model epi.dec`
  temptable <- matrix(NA,250,1000)
  for(k in 250:500){
    for(i in 1:1000){
      temptable[,i] <- cur.res[[i]]$sq.err[[k]]
    }
    dec1.mse[,j] <- rowMeans(temptable)
  }
}


#make a table to compare#
dec.001.mean <- c()
for(i in 1:5){
  dec.001.mean[i] <- mean(dec001.mse[,i])
}
dec.01.mean <- c()
for(i in 1:5){
  dec.01.mean[i] <- mean(dec01.mse[,i])
}
dec.1.mean <- c()
for(i in 1:5){
  dec.1.mean[i] <- mean(dec1.mse[,i])
}

comp.dec.mse <- matrix(c(dec.001.mean,dec.01.mean,dec.1.mean), byrow=FALSE,ncol=3)
colnames(comp.dec.mse)<-c(0.001, 0.01, 0.1)
rownames(comp.dec.mse)<-c(5,10,15,20,25)


##########epi.sign####

#0.01
sign01.mse <- as.data.frame(matrix(NA,250,5))
colnames(sign01.mse)<-c(5,10,15,20,25)

for(j in 1:5){
  cur.res <- results[[j]]$`prob.change 0.01`$`model epi.sign`
  temptable <- matrix(NA,250,1000)
  for(k in 250:500){
    for(i in 1:1000){
      temptable[,i] <- cur.res[[i]]$sq.err[[k]]
    }
    sign01.mse[,j] <- rowMeans(temptable)
  }
}

#0.001
sign001.mse <- as.data.frame(matrix(NA,250,5))
colnames(sign001.mse)<-c(5,10,15,20,25)

for(j in 1:5){
  cur.res <- results[[j]]$`prob.change 0.001`$`model epi.sign`
  temptable <- matrix(NA,250,1000)
  for(k in 250:500){
    for(i in 1:1000){
      temptable[,i] <- cur.res[[i]]$sq.err[[k]]
    }
    sign001.mse[,j] <- rowMeans(temptable)
  }
}

#0.1
sign1.mse <- as.data.frame(matrix(NA,250,5))
colnames(sign1.mse)<-c(5,10,15,20,25)

for(j in 1:5){
  cur.res <- results[[j]]$`prob.change 0.1`$`model epi.sign`
  temptable <- matrix(NA,250,1000)
  for(k in 250:500){
    for(i in 1:1000){
      temptable[,i] <- cur.res[[i]]$sq.err[[k]]
    }
    sign1.mse[,j] <- rowMeans(temptable)
  } 
}


#make a table to compare#
sign.001.mean <- c()
for(i in 1:5){
  sign.001.mean[i] <- mean(sign001.mse[,i])
}
sign.01.mean <- c()
for(i in 1:5){
  sign.01.mean[i] <- mean(sign01.mse[,i])
}
sign.1.mean <- c()
for(i in 1:5){
  sign.1.mean[i] <- mean(sign1.mse[,i])
}

comp.sign.mse <- matrix(c(sign.001.mean,sign.01.mean,sign.1.mean), byrow=FALSE,ncol=3)
colnames(comp.sign.mse)<-c(0.001, 0.01, 0.1)
rownames(comp.sign.mse)<-c(5,10,15,20,25)


