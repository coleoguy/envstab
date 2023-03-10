# Originally J. Strope now Andres and Heath
#1/13/22

######### INTERNAL FUNCTIONS ########################

#function that gets fitness of each individual
GetFit <- function(genome,
                    imp.loci, # what genes we keep track of and
                    # what will impact phenotype
                    fav.pheno,
                    model){ #dominance factor
  # creating a vector in loci.w
  loci.w <- c()
  if(model == "additive"){
    pheno <- sum(genome[,imp.loci])
    w <- 1 - (abs(fav.pheno - pheno)/20)
  }
  if(model ==  "epi.sign"){
    pheno <- 0
    inv.eff <- colSums(genome[, imp.loci[1:5]])
    for(i in 1:length(inv.eff)){
      if(inv.eff[i] == 0){
        pheno <- pheno + sum(genome[ , imp.loci[(5 + i)]]) * 2
      }
      if(inv.eff[i] == 1){
        pheno <- pheno +  2
      }
      if(inv.eff[i] == 2){
        pheno <- pheno + sum(as.numeric(!genome[ , imp.loci[(5 + i)]])) * 2
      }
    }
  }
  if (model== "epi.inc"){
    x <- sum(genome[,imp.loci])
    pheno <- (20*(x/20)^2)
  }
  if(model== "epi.dec"){
    x <- sum(genome[,imp.loci])
    pheno <- (13*log(x+1)/2)
  }
  w <- 1 - (abs(fav.pheno - pheno)/20)
  return(w)
  #returns product(multiplication) of fitness locations
}

GetMeanPheno <- function(population, imp.loci, model){
  cur.phenos <- c()
  for(j in 1:length(population)){
    genome <- population[[j]]
    if(model == "additive"){
      pheno <- sum(genome[,imp.loci])
    }
    if(model ==  "epi.sign"){
      pheno <- 0
      inv.eff <- colSums(genome[, imp.loci[1:5]])
      for(i in 1:length(inv.eff)){
        if(inv.eff[i] == 0){
          pheno <- pheno + sum(genome[ , imp.loci[(5 + i)]]) * 2
        }
        if(inv.eff[i] == 1){
          pheno <-  pheno + 2
        }
        if(inv.eff[i] == 2){
          pheno <- pheno + sum(as.numeric(!genome[ , imp.loci[(5 + i)]])) * 2
        }
      }
    }
    if (model== "epi.inc"){
      x <- sum(genome[, imp.loci])
      pheno <- (20*(x/20)^2)
      }
    if(model== "epi.dec"){
      x <- sum(genome[, imp.loci])
      pheno <- (13*log(x+1)/2)
      }
    cur.phenos[j] <- pheno
  }
  return(mean(cur.phenos))
}

# get initial population structure
GetPopulation <- function(loci, popsize){
  #stores population in list
  population <- list()
  #generates random uniform distribution to be the favored allele
  p <- runif(100)
  #calculates unfavored allele based on the frequency of favored allele
  q <- 1-p
  #loop assesing 1 to the population size
  for(i in 1:popsize){
    genome <- matrix(data = NA,
                     nrow = 2, #2 rows
                     ncol = loci) #number of columns is number of loci
    for(j in 1:loci){
      genome[, j] <- sample(x = c(0, 1), #0 or 1 allele
                            size = 2, #diploid
                            prob = c(p[j], q[j]), #probablility that it will be 
                            #favored or unvavored
                            replace = T) #can resue 0,1, and 2
    }
    #genome stored for each individual in the population
    population[[i]] <- genome
  }
  return(population)
}

# function that gives recombination spots when given the chromosome 
# number and loci
GetRecomb <- function(chrom.num,
                      loci){
  #stores loci/ chromosome number in list depicting each chromosome size. floor 
  #rounds down so there is whole numbers
  chrom.size <- floor(loci/chrom.num)
  #loci%% is the leftover from the rounding of the fractions. it is added to the
  #end of the chromosome. because of this, the last chromosome will be longest
  extra.genes <- loci%%chrom.num
  #gen.struct stores chromosome length info. chrom.num is number of rows, 2 is 
  #number of columns
  
  # Run on a single chrom.num
  gen.struct <- matrix(NA, chrom.num, 2)
  for(i in 1:nrow(gen.struct)){
    # first column is length of first chromosome
    gen.struct[i, 1] <- chrom.size*(i-1)+1
    #second column is length of second chromosome
    gen.struct[i, 2] <- chrom.size*i
  }
  #storing loci spots in both chromosomes
  gen.struct[chrom.num, 2] <- loci
  #create a vector for rspots
  rspots <- c()
  #loop from 1 to the number of chromosomes
  for(i in 1:chrom.num){
    #this finds the recombination spots and stores it in rspots. finds it by 
    #sampling the first column [i,1] and the second column [i,2] and using +1 to 
    #indicate that recombination is happening at the loci after the division spot
    rspots[i] <- sample((gen.struct[i, 1]+1):(gen.struct[i, 2]-1), 1)
  }
  #get recombination spots and length of genome when use get recomb function
  results <- list(rspots, gen.struct)
  #return recombination spots and gene structure from this function
  return(results)
}
# function that returns the opposite value that is being used. 
# this is used to switch to the other chromosome
GetOp <- function(x){
  if(x == 1) return(2)
  if(x == 2) return(1)
}

GetGametes <- function(population, popsize, chrom.num, loci, fitnesses){
  #in parents list, sample 1- populationsize , size is twice the population size, 
  #replace= T becuase parents can have more than one offspring, probablility of 
  #parents reproducing is a function of fitness
  parents <- sort(sample(1:popsize, popsize*2, replace=T, prob=fitnesses))
  # get the chromosome structure for the genome
  
  #stores recomb spots for second chromosome in gen.info list ### INCORRECT
  gen.info <- GetRecomb(chrom.num, loci)[[2]]
  #each time we make a gamete it gets stored in a list
  gametes <- list()
  #loop
  for(j in 1:length(parents)){
    #this will work its way through each of the most fit parents we sampled
    foo <- population[[parents[j]]]
    # get the recombination points for a single meiosis
    #stores recombination spots for chromosome one in recomb.spots list ### INCORRECT
    recomb.spots <- GetRecomb(chrom.num, loci)[[1]]
    #start on row 1 or 2, sample all the chromosomes, can be used more than once
    s.strand <- sample(1:2, chrom.num, replace=T)
    #stores gamete in vector
    gamete <- c()
    for(i in 1:chrom.num){
      #this tells us to sample the strand up to the recombination spot then switch strands at these assigned spots
      #first starts at the beginning of the chromosome then goes to the recombination spot, 
      #switches chromosome then samples from the spot after the recombination spot up to the 
      #next stopping point the repeats up to the end of the chromosome
      gamete <- c(gamete, c(foo[s.strand[i], gen.info[i,1]:recomb.spots[i]],
                            foo[GetOp(s.strand[i]),(recomb.spots[i]+1):gen.info[i,2]]))
    }
    gametes[[j]] <- gamete
  }
  return(gametes)
}


Mutate <- function(population, mut.rate){
  mutants <- sample(1:popsize, popsize*mut.rate)
  #Uses the mutation rate to sample out mutants from total population by their index
  for(i in mutants){
    #### TODO make mutations less good (old note ~A)
    population[[i]][sample(1:2, 1), sample(1:loci, 1)] <- sample(0:1, 1)
    # population[[i]] picks the selected mutant
    # The first sample selects one of the two parental copies
    # The second sample selects which of the loci is going to be affected
    # The third sample that is assigned is the actual mutation being applied,
    ## changing the current allele (potentially to the original value) there
    ## are not deleterious mutations in this model
  }
  return(population)
}




######### INTERNAL FUNCTIONS ########################


simulate <- function(loci, chrom.num, popsize, generations, 
                    prob.change, mut.rate, model){
  population <- GetPopulation(loci, popsize)
  imp.loci <- sort(sample(1:100, 10))
  #assesses what the favored genotype is
  fav.pheno <- runif(1, min=0, max=20)
  #selection coeficcient- .01 is the minimum value for fitness and .4 is 
  # maxiumum value. runif() generates a random uniform distribution. 10 is 
  # the number of random samples

  mean.fitness <- c()
  mean.phenotypes <- c()
  cur.opt.pheno <- c()
  for(j in 1:generations){
    print(paste("running generation", j))
    #store fitnesses in a vector
    fitnesses <- c()
    if(runif(1) <= prob.change){
      fav.pheno <- runif(1, min=0, max=20)
    }
    #loop to get fitness for each individual in the population and storing it 
    # in matrix "fitnesses". the i makes it run through each specific 
    # indivial seperately
    for(i in 1:popsize){
      fitnesses[i] <- GetFit(genome = population[[i]],
                             imp.loci = imp.loci,
                             fav.pheno = fav.pheno,
                             model = model)
    }
    # return the mean of fitnesses from the loop and store it in mean.fitnesses 
    mean.fitness[j] <- mean(fitnesses)
    mean.phenotypes[j] <- GetMeanPheno(population, imp.loci, model)
    cur.opt.pheno[j] <- fav.pheno
    
    
    
    # mutations
    population <- Mutate(population, mut.rate)
    
    # gametogenesis
    #getting gametes when input information and then storing it in matrix "gametes"
    gametes <- GetGametes(population, popsize, chrom.num, loci, fitnesses)
    # fertilization
    #loop to bind the rows and make it into one row in a matrix- combines gametes 
    #to make a single individual. i used to go through each specific individual in loop
    for(i in 1:popsize){
      population[[i]] <- rbind(gametes[[i]],
                               gametes[[popsize+i]])#use 1000+i because there are 2000 gametes
    }
  }
  sq.err <- c()
  for(j in 1:generations){
    sq.err[j] <- (mean.phenotypes[j]-cur.opt.pheno[j])**2
  }
  x <- list(mean.fitness, mean.phenotypes, cur.opt.pheno, sq.err)
  names(x) <- c("mean.fitness","mean.phenotypes", "opt.pheno", "sq.err")
  return(x)
  
}




