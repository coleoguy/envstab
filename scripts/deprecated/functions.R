#juliette strope
#9/15/21

#function that gets fitness of each individual
GetFit <- function(genome,
                   imp.loci, # what genes we keep track of and
                             # what will impact phenotype
                   fav.geno, #favored genotype
                   s, #selection coefficient
                   h){ #dominance factor
  # creating a vector in loci.w
  loci.w <- c()
  for(i in 1:length(imp.loci)){
  	#the sum of the loci that will affect fitness geno is numeric 0:2
    geno <- sum(genome[,imp.loci[i]])
   #if favored genotype is 0, then the unfavored is assigned to 2
    if(fav.geno[i]==0){
      #as.character returns a string of 0,1,and 2 dependent on the statement provided
      #switch compares values in string and gives the correct fitness for that locus
      #assessing whether each loci is favored or not
      loci.w[i] <- switch(as.character(geno),
                          "0" = 1,
                          "1" = 1-h[i]*s[i],
                          "2" = 1-s[i])
    }
    #if favored genotype is 2, then the unfavored is assigned to 0
    if(fav.geno[i]==2){
      loci.w[i] <- switch(as.character(geno),
                          "0" = 1-s[i],
                          "1" = 1-h[i]*s[i],
                          "2" = 1)
    }
  }
  return(prod(loci.w))
  #returns product(multiplication) of fitness locations
}

#function that gets fitness of each individual
GetFit2 <- function(genome,
                    imp.loci, # what genes we keep track of and
                    # what will impact phenotype
                    fav.pheno, #favored genotype
                    s, #selection coefficient
                    h){ #dominance factor
  # creating a vector in loci.w
  loci.w <- c()
  pheno <- sum(genome[,imp.loci])
  w <- 1 - (abs(fav.pheno - pheno)/20)
  return(w)
  #returns product(multiplication) of fitness locations
}

GetMeanPheno <- function(population){
  phenos <- c()
  for(i in 1:length(population)){
    phenos[i] <- sum(population[[i]][,imp.loci])
  }
  pheno <- mean(phenos)
  return(pheno)
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
                            prob = c(p[j], q[j]), #probablility that it will be favored or unvavored
                            replace = T) #can resue 0,1, and 2
    }
    #genome stored for each individual in the population
    population[[i]] <- genome
  }
  return(population)
}

# function that gives recombination spots when given the 
# chromosome number and loci
GetRecomb <- function(chrom.num,
                      loci){
  #stores loci/ chromosome number in list depicting each chromosome size.
  #floor rounds down so there is whole numbers
  chrom.size <- floor(loci/chrom.num)
  #loci%% is the leftover from the rounding of the fractions. it is added to
  #the end of the chromosome. because of this, the last chromosome will be longest
  extra.genes <- loci%%chrom.num
  #gen.struct stores chromosome length info. chrom.num is number of rows,
    #2 is number of columns
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
#function that returns the opposite value that is being used. this is used to 
# switch to the other chromosome
GetOp <- function(x){
  if(x == 1) return(2)
  if(x == 2) return(1)
}


GetGametes <- function(population, popsize, chrom.num, loci, fitnesses){
	#in parents list, sample 1- population size , size is twice the population size,
  #replace= T because parents can have more than one offspring, probability of
  #parents reproducing is a function of fitness
  parents <- sample(1:popsize, popsize*2, replace=T, prob=fitnesses)
  # get the chromosome structure for the genome

 #stores recomb spots for second chromosome in gen.info list
  gen.info <- GetRecomb(chrom.num, loci)[[2]]
  #each time we make a gamete it gets stored in a list
  gametes <- list()
  #loop
  for(j in 1:length(parents)){
    #this will work its way through each of the most fit parents we sampled
    foo <- population[[parents[j]]]
    # get the recombination points for a single meiosis
    #stores recombination spots for chromosome one in recomb.spots list
    recomb.spots <- GetRecomb(chrom.num, loci)[[1]]
    #start on row 1 or 2, sample all the chromosomes, can be used more than once
        s.strand <- sample(1:2, chrom.num, replace=T)
  #stores gamete in vector
    gamete <- c()
    for(i in 1:chrom.num){
    #this tells us to sample the strand up to the recombination spot then
    #switch strands at these assigned spots
    #first starts at the beginning of the chromosome then goes to the
    #recombination spot, switches chromosome then samples from the spot after
    #the recombination spot up to the next stopping point the repeats up to
    #the end of the chromosome
      gamete <- c(gamete, c(foo[s.strand[i], gen.info[i,1]:recomb.spots[i]],
                            foo[GetOp(s.strand[i]),(recomb.spots[i]+1):gen.info[i,2]]))
    }
    gametes[[j]] <- gamete
  }
  return(gametes)
}


