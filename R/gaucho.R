#' Genetic Algorithm for Understanding Clonal Heterogeneity and Ordering (GAUCHO)
#' 
#' Use a genetic algorithm to find the relationships between the values in an input file - the package was written to deal with single
#' nucleotide variants (SNVs) in mixtures of cancer cells, but it will work with any mixture.  It will calculate appropriate phylogenetic
#' relationships between clones them and the proportion of each clone that each sample is composed of.  For detailed usage, please
#' read the accompanying vignette.
#' 
#' @param observations              Observation data frame where each row represents an SNV and each column represents
#' a discrete sample separated by time or space.  Note that the data frame must have column names and row names.  
#' Every value must be a proportion between 0 and 1.  See details
#' @param number_of_clones          An integer number of clones to be considered
#' @param pop_size                  The number of individuals in each generation
#' @param mutation_rate             The likelihood of each individual undergoing mutation per generation
#' @param iterations                The maximum number of generations to run
#' @param stoppingCriteria          The number of consecutive generations without improvement that will stop the algorithm.
#' Default value is 20\% of iterations.
#' @param parthenogenesis           The number of best-fitness individuals allowed to survive each generation
#' @param nroot                     Number of roots the phylogeny is expected to have.When nroot=0, a random integer between 1 and the number of clones is generated for each phylogeny
#' @param contamination             Is the input contaminated?  If set to 1, an extra clone is created in which to place inferred contaminants
#' @param check_validity            Unless set to false, eliminate any clones with no new mutations, disallow those clones. Increases computational overheads.
#' @details The input data should be a data.frame containing proportions of cells that contain a feature.  There are a number of ways
#' to create these data, including merging the balance of alleles and copy number of an SNV using the equation min(1,r*CN/(r+R)), where
#' CN is the copy number, r is the number of non-reference reads and R is the number of reference reads.  For example, if a site
#' were sequenced to a depth of 100x, with 25 non-reference reads and 75 reference reads and diploid copy number, the result would be
#' min(1,25*2/(25+75)) = 0.5.  Therefore, 50\% of the cells in the sample contain the SNV.
#' Further details are available in the accompanying vignette.
#' @return Returns an object of class ga \code{\link{ga-class}}.  Note that the number of clones and number of cases are 
#' stored in the unused min and max slots of the output object.
#' @export
#' @import compiler GA graph heatmap.plus png Rgraphviz
# @references gaucho paper here!
#' @author Alex Murison \email{Alexander.Murison@@icr.ac.uk} and Christopher Wardell \email{Christopher.Wardell@@icr.ac.uk}
#' @seealso \code{\link{ga-class}}, \code{\link{ga}}, \code{\link{gauchoReport}}, \code{\link{gaucho_simple_data}},
#' \code{\link{gaucho_hidden_data}}, \code{\link{gaucho_synth_data}},
#' \code{\link{gaucho_synth_data_jittered}}, \code{\link{BYB1_G07_pruned}}
#' @examples
#' ## The vignette provides far more in-depth explanation and examples ##
#'   
#' ## Load the included simple example data
#' gaucho_simple_data = read.table(file.path(system.file("extdata",package="gaucho"),"gaucho_simple_data.txt"),header=TRUE,row.names=1)
#' 
#' ## Run gaucho using 3 clones and a phylogeny with a single root
#' solution=gaucho(gaucho_simple_data, number_of_clones=3,nroot=1,iterations=1000)
#' 
#' ## Create the four output plots
#' gauchoReport(gaucho_simple_data,solution,outType="fitness")
#' gauchoReport(gaucho_simple_data,solution,outType="heatmap")
#' gauchoReport(gaucho_simple_data,solution,outType="phylogeny")
#' gauchoReport(gaucho_simple_data,solution,outType="proportion")
#' 
#' ## Output the solution and plots in the current working directory
#' # gauchoReport(gaucho_simple_data,solution)



gaucho<-function(observations, number_of_clones, pop_size=100, mutation_rate=0.8, iterations=1000,
               stoppingCriteria=round(iterations/5), parthenogenesis=2,nroot=0, contamination=0, check_validity=TRUE) {

  ## Load all libraries
  
  ##############################
  ## Start internal functions ##
  ##############################
  
  ## Crossover as gabin_SpCrossover but do not select a crossover point within a phylogeny
  phylo_cross<-function(object, parents) {
    fitness <- object@fitness[parents]
    parents <- object@population[parents, , drop = FALSE]
    n <- ncol(parents)
    children <- matrix(NA, nrow = 2, ncol = n)
    fitnessChildren <- rep(NA, 2)
    # Here we make sure that the crossover cannot happen in the phylogeny
    crossOverPoint <- sample((number_of_clones+1):n, size = 1)
    if (crossOverPoint == 0) {
      children[1:2, ] <- parents[2:1, ]
      fitnessChildren[1:2] <- fitness[2:1]
    }
    else if (crossOverPoint == n) {
      children <- parents
      fitnessChildren <- fitness
    }
    else {
      children[1, ] <- c(parents[1, 1:crossOverPoint], parents[2, 
                                                               (crossOverPoint + 1):n])
      children[2, ] <- c(parents[2, 1:crossOverPoint], parents[1, 
                                                               (crossOverPoint + 1):n])
      fitnessChildren <- NA
      # Check this hasn't created empty clones
      
    }
    if (check_validity) {
      while (is_invalid(children[1,]) || is_invalid(children[2,])) {
        crossOverPoint <- sample((number_of_clones+1):n, size = 1)
        if (crossOverPoint == 0) {
          children[1:2, ] <- parents[2:1, ]
          fitnessChildren[1:2] <- fitness[2:1]
        }
        else if (crossOverPoint == n) {
          children <- parents
          fitnessChildren <- fitness
        }
        else {
          children[1, ] <- c(parents[1, 1:crossOverPoint], parents[2, 
                                                                   (crossOverPoint + 1):n])
          children[2, ] <- c(parents[2, 1:crossOverPoint], parents[1, 
                                                                   (crossOverPoint + 1):n])
          fitnessChildren <- NA
          # Check this hasn't created empty clones
          
        }
      }
    }
    out <- list(children = children, fitness = fitnessChildren)
    return(out)
  }
  ## /end of phylo_cross

  ## Check whether any of the clones have no new mutations
  is_invalid<-function(input) {
    
    #Get the list of mutations
    mutation_list<-input[((number_of_clones)+(pseudo_number_of_clones*number_of_cases)+1):length(input)]
    
    #multiply 1 by the number of mutations new to each clone
    is_it_true<-1
    for (current_clone in 1:number_of_clones) {    
      is_it_true<-length(which(mutation_list==current_clone))*is_it_true
    }
    
    #if the result is 0 then return TRUE (invalid) or return FALSE (valid)
    if (is_it_true==0) { 
      return(TRUE)
    }    else {return(FALSE)}
    
  }
  
  ## Calculate fitness of phylogeny
  fit_phylogeny<-function(input) {
    
    ### Work out the phylogeny matrix
    dep_mat<-generate_phylo_matrix(input[1:number_of_clones],number_of_clones)
    
    ### Calculate the proportions of each clone in each case
    proportion_matrix<-matrix(input[(number_of_clones+1):((number_of_cases*pseudo_number_of_clones)+number_of_clones)], 
                              ncol=pseudo_number_of_clones, nrow=number_of_cases, byrow=TRUE)
    for (rows in 1:nrow(proportion_matrix)) {
      scale<-sum(proportion_matrix[rows,])
      proportion_matrix[rows,]<-proportion_matrix[rows,]/scale
    }
    
    tot_score<-0
    
    # For each mutation, obtain the clones it appears in from the phylogeny matrix
    # Calculate the expected frequency in each case by summing clone frequencies for each mutation in each case
    prediction_matrix<-matrix(ncol=number_of_cases, nrow=number_of_mutations)
    rownum<-1
    for (mutation in input[((number_of_clones)+(pseudo_number_of_clones*number_of_cases)+1):length(input)]) {
      affected_clones<-dep_mat[mutation,]
      for (case in 1:number_of_cases) {
        score<-0
        for (clone in 1:number_of_clones) {         
          score<-score+(affected_clones[clone]*proportion_matrix[case,clone]) ## 16% of total load
        }
        prediction_matrix[rownum,case]<-score
        
      }
      rownum<-rownum+1
    }
    
    # calculate the deviance for each case/mutation from the observations
    # return -1 times the total deviance for all mutations
    vals<-prediction_matrix-observation_matrix ## 11% of total load
    for (plink in 1:nrow(vals)) {
      tot_score<-tot_score+sum(abs(vals[plink,])) ## 62% of total load
    }
    
    return(-1*tot_score)
    
  }
  # /end of fit_phylogeny
  
  ## Generates a population of individuals
  generate_phylogeny_aware_population<-function(object) {
    #print("popgen")
    population <- matrix(NA, nrow = object@popSize, ncol = object@nBits)
    for (j in 1:object@popSize) {
      population[j, ] <- generate_phylogeny_aware_individual()
      if (check_validity) {
        while (is_invalid(population[j,])) {
        population[j, ] <- generate_phylogeny_aware_individual()
        }
      }
    }
    return(population)
  }
  # /end of generate_phylogeny_aware_population
  
  ## Generates an individual
  generate_phylogeny_aware_individual<-function() {
    #  print("indivgen")
    
    # work out how many proportions
    clo_cas<-(pseudo_number_of_clones*number_of_cases)
    #work out how many mutations
    how_big<-number_of_clones+clo_cas+number_of_mutations
    # create our new individual
    new_individual<-vector(length=how_big)
    
    #Assign parents from 0 to number of clones
    new_individual[1:number_of_clones]<-new_phylo(nroot,number_of_clones)
    # assign proportions between 0 and 5
    new_individual[(number_of_clones+1):(number_of_clones+clo_cas)]<-round(runif(clo_cas,0,5))
    # assign mutation first appears as 1:number of clones
    new_individual[(number_of_clones+clo_cas+1):how_big]<-round(runif(number_of_mutations,0.5,(number_of_clones+0.49)))
    # return new individual
    return(new_individual)
  }
  # /end of generate_phylogeny_aware_individual
  
  ## Mutation function
  phylo_mutate<-function(object, parent) {
    
    #  print("mutator")
    # Turn the string into a vector
    victim_ready <- parent <- as.vector(object@population[parent, ])
    #<-as.numeric(unlist(strsplit(unsuspecting_victim, "")))
    
    # work out where we are going to mutate - 
    # 0 a proportion of clonal population
    # 1 mutation presence/absence in a clone
    # Randomly select with equal probability of both
    what<-round(runif(1,0,2))
    if (what == 0) {
      # Mutate Proportions
      
      ## We allow multiple mutations of the proportions string
      proportion_mutations=round(runif(1,1,number_of_clones))
      
      for(proportion_mutation in 1:proportion_mutations){
        # Randomly select a proportion
        where<-round(runif(1,(number_of_clones+0.5),((number_of_cases*pseudo_number_of_clones)+number_of_clones+0.49)))
        new<-victim_ready[where]
        
        ## Define what the increment or decrement will be:
        decide_delta=runif(1,0,1)
        if(decide_delta<=0.5){delta=1}
        if(decide_delta>0.5 & decide_delta<0.8){delta=2}
        if(decide_delta>=0.8){delta=3}
        
        # if proportion is already 0 you may ONLY increase it
        if (new-delta<0) {
          new<-new+delta
        }else{ 
          # Otherwise randomly increment or decrement by 1
          do_what<-round(runif(1,0,1))
          if (do_what==1) {new<-new+delta}
          if (do_what==0) {new<-new-delta}
        }
        victim_ready[where]<-new
        
      } # end of proportion_mutations loop
      
    }
    else if (what == 1 ) {
      # Mutate parents
      
      victim_ready[1:number_of_clones]<-new_phylo(nroot,number_of_clones)
    }
    else if (what == 2 ) {
      # Mutate a first_appears call
      where<-round(runif(1,((number_of_cases*pseudo_number_of_clones)+number_of_clones+0.5),length(victim_ready)))
      old<-victim_ready[where]
      victim_ready[where]<-round(runif(1,0.5,(number_of_clones+0.49)))
      if (check_validity) {
        while (is_invalid(victim_ready)) {
          victim_ready[where]<-old
          where<-round(runif(1,((number_of_cases*pseudo_number_of_clones)+number_of_clones+0.5),length(victim_ready)))
          old<-victim_ready[where]
          victim_ready[where]<-round(runif(1,0.5,(number_of_clones+0.49)))
        }
      }
    }
    
    # Return it
    return(victim_ready)
  }
  # /end of phylo_mutate
  
  
  ############################
  ## End internal functions ##
  ############################
  
  
  #########################
  ## Start of main logic ##
  #########################
  
  ## Get Observation File and associated values
  ## We ignore the first column, as it is the object names
  #observation_matrix<-observations[,2:ncol(observations)] 
  observation_matrix<-observations
  number_of_mutations<-as.numeric(nrow(observation_matrix))
  number_of_cases<-as.numeric(ncol(observation_matrix))

  if (contamination == 1) {
    pseudo_number_of_clones<-number_of_clones+1
    numBits<-(pseudo_number_of_clones*number_of_cases)+(number_of_mutations+number_of_clones)    
  } else {
    numBits<-(number_of_clones*number_of_cases)+(number_of_mutations+number_of_clones)
    pseudo_number_of_clones<-number_of_clones
  }
  
  ### Compile the genetic algorithm - this increases speed significantly
  setCompilerOptions(suppressAll=TRUE)
  ga=cmpfun(ga)
    
  ### RUN THE GENETIC ALGORITHM
  goo<-ga(type="binary", 
          fitness = fit_phylogeny, 
          population = generate_phylogeny_aware_population,
          selection = gabin_tourSelection,
          crossover = phylo_cross,
          mutation = phylo_mutate,
          popSize = pop_size,
          nBits = numBits,
          pmutation=mutation_rate,
          maxiter=iterations,
          elitism = parthenogenesis,
          run=stoppingCriteria
  )
  ## Add annotation to the object
  #goo@names=as.character(observations[,1])
  goo@names=rownames(observations)
    
  ## Add number of clones to min slot.  Note that this is an incorrect use of the slot
  goo@min = number_of_clones
  ## Add number of cases (e.g. timepoints) to max slot.  Note that this is an incorrect use of the slot
  goo@max= number_of_cases

  return(goo)

  #######################
  ## End of main logic ##
  #######################
}
