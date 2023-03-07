# Heath Blackmon
# coleoguy@gmail.com
# example of popgen simulation

loci <- 100
chrom.num <- 3
popsize <- 1000
generations <- 300
#source("functions.R")
#function gets initial population structure when given loci and population size 
#and stores it in the list "population"
population <- GetPopulation(loci, popsize)


#set.seed makes it able to be reproduced
# set.seed(1)
#imp.loci is what genes we will keep track of and what will be impacted(no replace
#statement bc we are not reusing genes)
imp.loci <- sample(1:100, 10)
#assesses what the favored genotype is
##fav.geno is not used in the workflow, is only called by GetFit function which is not used
fav.geno <- sample(c(0,2), 10, replace = T)
##fav.pheno is used by GetFit 2 function, which is used to calculate fitness of individuals
## is changed to new value after 150 generations (half the run?) in line 43
fav.pheno <- runif(1, min=0, max=20)
#selection coefficient- .01 is the minimum value for fitness and .4 is maximum
#value. runif() generates a random uniform distribution. 10 is the number of random samples
s <- runif(10, min=.01, max=.4)
#dominance factor- 0.5 used because traits are additive
h <- rep(.5, 10)


# assess fitness
mean.fitness <- c()
mean.phenotypes <- c()

#starting a loop to get mean fitness over the span of multiple generations
for(j in 1:generations){
  print(paste("running generation", j))
  #store fitnesses in a vector
  fitnesses <- c()
  if(j == 150){
    fav.pheno <- runif(1, min=0, max=20)
  }
  #loop to get fitness for each individual in the population and storing it in
  #matrix "fitnesses". the i makes it run through each specific individual separately
  for(i in 1:popsize){
    fitnesses[i] <- GetFit2(genome = population[[i]],
                           imp.loci = imp.loci,
                           fav.pheno = fav.pheno,
                           s = s,
                           h = h)
  }
#return the mean of fitnesses from the loop and store it in mean.fitnesses
#(I don't know why j is used though?is it to differentiate it from the previous mean.fitnesses?)
  mean.fitness[j] <- mean(fitnesses)
  mean.phenotypes[j] <- GetMeanPheno(population)
  # gametogenesis
  #getting gametes when input information and then storing it in matrix "gametes"
  gametes <- GetGametes(population, popsize, chrom.num, loci, fitnesses)

  # fertilization
  #loop to bind the rows and make it into one row in a matrix- combines gametes
  #to make a single individual. i used to go through each specific individual in loop
  for(i in 1:popsize){
    population[[i]] <- rbind(gametes[[i]],
                             gametes[[1000+i]])
  #use 1000+i because there are 2000 gametes
  }
  # mutations
  #sampling a random 10 % of the population size to be considered mutations
  mutants <- sample(1:popsize, popsize*.1)
  #loop to pick random places in the genome for mutations to occur
  for(i in mutants){
  	#population[[i]] exists to give us genome of each individual 	that will be mutated
  	#sample one time to determine randomly if it will be row one 	or row two
  	#sample again one time to find a random column(from 1:loci)
    population[[i]][sample(1:2, 1), sample(1:loci, 1)] <- sample(0:1, 1)
  } # ^ randomly assigns 0 or 1 to a site

}
#plots the impact mutations have on fitness of the population
wmutes <- mean.fitness
plot(x=1:generations, y=mean.fitness,
     xlab="generation",
     ylab="fitness",
     type="l")
plot(x=1:generations, y=mean.phenotypes,
     xlab="generation",
     ylab="phenotype",
     type="l")

lines(x=1:generations, y=wmutes, col="blue")

# J
# 1) Find out how it crashes when we change
# chromosome number (2-10)
# 2) Look at how the number of loci impacts the
# the rate of increase in fitness

