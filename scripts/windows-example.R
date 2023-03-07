# How to run Rcode in parallel 
# if you are using this parallel script for the first time install the following
# packages. otherwise you can skip this step
# the installation step should take care of the dependencies which are
# foreach, snow and iterators
install.packages("doSNOW") 

# load libraries
library(doSNOW) # this is the library for parallel computing

# define number of clusters for parallel computing
# number of CPUs
nClust <- 4
# make connections for each CPU
# setting outfile = "" would print the output on the console
# however in Windows OS you will not see an output unless you
# run the R script through command prompt,powershell, or
# windows subsystem for linux
cl <- makeCluster(nClust, outfile = "")
registerDoSNOW(cl)

# now you are ready to run the code in parallel
# using the code example below you can run any script in parallel
# the below code is similar to a for loop but with some exceptions
# you would set the iterations as you would in a for loop

# setting .verbose = T will be useful for troubleshooting
# this parameter outputs if the run failed or succeed and
# the number of runs remaining

# you may need to call packages for each palatalization you do
# even if you called them before this loop
# for that use .packages parameter
# make sure include all the dependencies when you are loading a package
# within the loop.

# once the script is finished running the final result will be
# the last step of your script. Unfortunately you will not be able
# to get outputs from the intermediate steps.
x <- foreach(j = 1:4,                    
             .verbose = T) %dopar% {
              # you would see the output of print statement
              # if you run this code in a linux environment
              # or via windows powershell
             }
# once you finish your script you need to close connections to the 
# CPUs
# stop the cluster
stopCluster(cl)
# you could see that in the above example your output is the square root
# of y only. you do not get the value of y itself
