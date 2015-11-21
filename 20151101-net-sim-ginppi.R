
##########################################
#20151101 modified for bio386 XSEDE exercise
#2015Oct13, use numeric lookup table for essential genes.
#2014March3, redo ginppi simulation witht same parameter for ms02networks. 
# 2014 Feb 17, change name "20131221.DIPandGIN.sim.aging_v2.R" to "net-aging-sim-2014Feb17.R"
# 2013 Dec 20, merge DIP PPI and Genetic Inxt Net -> Multi-net approach
rm(list=ls())

myhost = 'gordon'  # 'byte' 'blacklight' 'mactower' 'greenfield'
mydir = '/oasis/scratch/hongqin/temp_project/Slim-network-aging-bio386-fall-2015'
#mydir = "/crucible/mc48o9p/hqin2/mactower-network-failure-simulation-master/ms02GINPPI"
if (myhost == 'byte') {  mydir = "/Users/hqin/github/Slim-network-aging-bio386-fall-2015"
} else if (myhost == 'helen') { mydir = "";  
} 
setwd(mydir)

single_network_failure_v2b = function(lambda1, lambda2=lambda1/10, threshold=4, p, pairs, essenLookupTb ) {
  #updated 20151101
  # single network failure simulation, 20151013Tue
  # lambda1: First exponential constant failure rate for edges with degree > threshold
  # lambda2: Second exponential constant failure rate for edges with degree <= threshold
  # threshold: degree threshold for lambda1 and lambda2
  # pairs: network in pairwide format, using numeric NOs 20151013
  # essenLookupTb: lookup table for essential and nonessential genes, numeric values 
  # Limitatation: For essential-essential gene interaction, 
  #    its ages were estimated were estiamted twice for each essential gene.
  ## for debug:   lambda1 = 1/50; lambda2= lambda1/10; threshold=4; p=0.99
  
  inpairs = pairs[,3:4] #Column 3,4 are numeric ids,   
  names(inpairs) = c('No1','No2')
  
  #get connectivities per node
  degreeTb = data.frame( table(c(inpairs$No1, inpairs$No2)))
  names(degreeTb) = c("No", "degree")
  
  # 20151101. increasing order of degrees
  # How does this influence essen-essen intxn age calculation? 
  degreeTb = degreeTb[ order(degreeTb$degree), ] 
  
  #initiate a empty storage for Modular-Ages
  degreeTb$moduleAge = NA;
  
  #loop over every node to calculate Modular-Ages
  for( i in 1:length(degreeTb[,1])){
    if ( essenLookupTb[ degreeTb$No[i] ]) { #essential node, do calculations
      lambda = ifelse( degreeTb$degree[i] >= threshold, lambda1, lambda2)
      age = rexp( degreeTb$degree[i], rate=lambda ) #assuming exponential decaying (non-aging intxn)
      if(degreeTb$degree[i] >= threshold){ # apply stochasticity
        active = runif(degreeTb$degree[i])  #uniform interaction stochasticity
        active = ifelse( active<=p, 1, NA  ) #pick active interactions
        if( sum(active, na.rm=T) > 0 ){ #there should be at least 1 active intxn
          age = age * active # only active interactions for modular age estimation
          degreeTb$moduleAge[i] = max(age, na.rm=T) #maximum intxn age is the module age
        } else {# when no active intxn is available 
          degreeTb$moduleAge[i] = 0; #this module is born dead.
        }
      } else { # for degree < threshold, no stochasticity is applied. 
        degreeTb$moduleAge[i] = max(age, na.rm=T) #maximum intxn age is the module age
      }
    } else {# non-essential node, skip-over with NA
      degreeTb$moduleAge[i] = NA 
    }
  }
  
  summary(degreeTb)
  currentNetworkAge = min(degreeTb$moduleAge, na.rm=T)
}

# R -f file --args lambda1 lambda2 degreeThreshold p popSize
# R -f 20151101-net-sim-ginppi.R --args 0.002 0.0002 4 0.99 5
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
lambda1 = as.numeric(args[1]); lambda1;
lambda2 = as.numeric(args[2]); lambda2;
degreeThreshold = as.integer((args[3])); degreeThreshold
p = as.numeric(args[4]); p;
popSize = as.numeric(args[5]); popSize;

list.files(pattern = 'csv' )
debug = 0; 
#for debug
#  lambda1=0.002; lambda2= lambda1/10; degreeThreshold=5; p=0.99; popSize=10; threshold=5

#essential gene info
 #essenTb = read.csv("SummaryRegressionHetHomFactorized2015Oct13.csv", colClasses=rep('character', 9))
 essenLookupTb = read.csv("essntialGeneLookupTable_20151013.csv", row.names=1)
 essenLookupTb = essenLookupTb[,1]

#pair-wise gene networks
 infile = "merged_PPIGIN_Factorized2015Oct13.csv"
 pairs = read.csv(infile)
 names(pairs) = c("id1",'id2', "No1", "No2")
 print(head(pairs))
 if(debug==9) {     pairs = pairs[1:1000,]  }
 pairs = pairs[ pairs$No1 != pairs$No2, ]  

# label essential nodes, remove nonesse-nonessen pairs
 pairs$essen1 = essenLookupTb[pairs$No1]
 pairs$essen2 = essenLookupTb[pairs$No2]
#remove nonessen <-> nonessen intxn because they do not affect aging. 
 pairs$remove = ifelse( pairs$essen1==0 & pairs$essen2==0, T, F  )
 pairs= pairs[! pairs$remove, ]  

# TODO: how many essen-essen intx? 
 
  #get connectivities per node
  degreeTb = data.frame( table(c(pairs$No1, pairs$No2)))
  summary(degreeTb); 
  degreeTb[1:10,]
  #median degree =5, mean=12
  #for one ms02, media =6, mean=12.02, so orginal network is power-law like, skew at two ends. 

  full_age_dir = "ori.ginppi.bio386"  
  
  popAges = numeric(popSize)
  time1 = date()
  j=1; count = 0; 
  while ((j <= popSize) && ( count < popSize*30)) {
    count = count + 1;      
    print(paste("count=",count))
    currentNetworkAge = single_network_failure_v2b(lambda1, lambda2, degreeThreshold, p, pairs, essenLookupTb)
    if (currentNetworkAge > 0) {
      popAges[j] = currentNetworkAge      
      j = j+1
    } #else do nothing for 'j'. Go to next 'count'.  
  }# end of j while-loop, population loop
    
  #output popAges      
  current_user = Sys.info()["user"]; host= Sys.info()["nodename"];
  timestamp = format(Sys.time(), "%Y%b%d_%H%M%S")
  age.file.name=paste("threshold", degreeThreshold, "p", p, "lambda1", lambda1, 
                         "lambda2", lambda2,'popsize',popSize,current_user, timestamp,host,"csv", sep="." )
  full_age_file = paste( full_age_dir,'/', age.file.name, sep='')
  write.csv( popAges, full_age_file, row.names=F)
      
