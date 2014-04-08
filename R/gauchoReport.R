#' View solutions contained within gaucho output
#' 
#' After running gaucho() on data, this function provides a convenient way to view the solutions and also export them as separate
#' text files and images.  For detailed usage, please read the accompanying vignette.
#' 
#' @param gauchoInput             Raw data analysed by gaucho() 
#' @param gauchoOutput            Object of class ga produced by gaucho()
#' @param outType               Type of output desired - must be one of the following: "complete","fitness","heatmap","phylogeny","proportion"
#' @param yRange                Y-axis range when plotting fitness of individuals.  Default is c(-250,0)
#' @param output_file_prefix    Optional prefix for all output files
#' @details This method reports data for the fittest individual; in the event of there being multiple individuals with identical
#' fitness, up to five individuals will be reported.  This function's output is governed by the outType argument.
#' All options except for the default "complete" value result in plotting the desired output to the current R session.
#' When outType=="complete", the following output is created for each individual:  the full length string, the phylogeny
#' matrix, the proportion matrix, the presence matrix, a heatmap of the raw data with the assigned clones as coloured bars
#' at the side, a stacked barplot showing the proportion of each clone at each timepoint and a plot showing the phylogenetic
#' relationship between the clones.  Note that the colours of the clones are consistent across all plots and that the contamination
#' clone (if present) is always the last clone.  Also produced is a plot illustrating the change in fitness as the generations evolved.
#' @return Nothing is returned.
#' @export
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


gauchoReport<-function(gauchoInput,gauchoOutput,outType="complete",yRange=c(-250,0),output_file_prefix="") {
  
  #################
  ## Input check ##
  #################
  
  ## Check that the outType argument is correctly specified
  if(!(outType %in% c("complete","fitness","heatmap","phylogeny","proportion"))){
    stop("outType incorrect: must be one of the following: \"complete\",\"fitness\",\"heatmap\",\"phylogeny\",\"proportion\"")
  }
  
  #####################
  ## End input check ##
  #####################
  
  ##################################################
  ## Fetch, calculate and otherwise set variables ##
  ##################################################
  
  ## Set observation matrix for heatmap
  observation_matrix = gauchoInput
  #observation_matrix = gauchoInput[,2:ncol(gauchoInput)]
  
  ## Get number of mutations from gauchoOutput
  number_of_mutations=length(gauchoOutput@names)
  
  ## Get number of clones from gauchoOutput
  number_of_clones=gauchoOutput@min
  
  ## Get number of cases (e.g. timepoints) from gauchoOutput
  number_of_cases=gauchoOutput@max
  
  ## We can now work out contamination status.  This is the reverse of calculating numBits in gaucho()
  if(((gauchoOutput@nBits-number_of_clones-number_of_mutations)/number_of_cases)==number_of_clones){
    contamination=0
  }else{
    contamination=1
  }
    
  ## Now we now the contamination status, we can calculate the pseudo_number_of_clones
  if (contamination == 1) {
    pseudo_number_of_clones<-number_of_clones+1
  } else {
    pseudo_number_of_clones<-number_of_clones
  }
      
  ## Set colour schemes
  heatmapColours=c("#313695","#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090","#fdae61","#f46d43","#d73027","#a50026") # blue yellow red 
  primary=c("#C9DD03","#FFD602","#F9A100","#EE7EA6","#A71930","#616365") # primary ICR colours
  secondary=c("#726E20","#6E273D","#F00034","#ADAFAF","#003D4C") # secondary ICR colours
  colours=c(primary,secondary[c(1,2,3,5)]) # exclude ICR light grey
  
  ######################
  ## End of variables ##
  ######################

  
  ########################################
  ## Begin looping over TOP 5 solutions ##
  ########################################
   
  ## Look at top solution(s) and parse into usable objects
  ## NOTE: if there are several equally good solutions, it will print (up to) the top five solutions
  top_number<-min(5,nrow(gauchoOutput@solution))
  for (number_of_answers in 1:top_number) {
    
    #######################################
    ## Begin plotting / producing output ##
    #######################################
    
    answers<-gauchoOutput@solution[number_of_answers,]
    phylogeny_matrix<-generate_phylo_matrix(answers[1:number_of_clones],number_of_clones)
    phylogeny_matrix<-rbind(phylogeny_matrix, answers[1:number_of_clones])
    if(contamination==1) {phylogeny_matrix=phylogeny_matrix[1:(nrow(phylogeny_matrix)-1),]}
    #colnames(phylogeny_matrix)=LETTERS[1:pseudo_number_of_clones]
    colnames(phylogeny_matrix)=LETTERS[1:number_of_clones]
    
    ## Only the proportion_matrix needs to include the pseudo_proportion_of_clones
    proportion_matrix<-matrix(answers[(number_of_clones+1):((number_of_cases*pseudo_number_of_clones)+number_of_clones)], 
                              ncol=pseudo_number_of_clones, nrow=number_of_cases, byrow=TRUE)
    colnames(proportion_matrix)=LETTERS[1:pseudo_number_of_clones] 
    for (rows in 1:nrow(proportion_matrix)) {
      scale<-sum(proportion_matrix[rows,])
      proportion_matrix[rows,]<-proportion_matrix[rows,]/scale
    }
    presence_matrix<-matrix(nrow=number_of_mutations, ncol=number_of_clones)
    for (rows in 1:number_of_mutations) {
      presence_matrix[rows,]<-phylogeny_matrix[answers[((number_of_cases*pseudo_number_of_clones)+number_of_clones+rows)],]
    }
    colnames(presence_matrix)=LETTERS[1:number_of_clones]
    
    ## Produce output text files only if outType is complete
    if(outType %in% c("complete")){
      message("Outputting the highest scoring solutions as text files")
      write.table(answers,paste0(output_file_prefix,number_of_clones,"clones.complete.solution",number_of_answers,"of",top_number,".txt"),sep="\t", row.names=FALSE,col.names=FALSE)
      write.table(phylogeny_matrix,paste0(output_file_prefix,number_of_clones,"clones.phylogeny.solution.",number_of_answers,"of",top_number,".txt"),sep="\t", row.names=FALSE)
      write.table(proportion_matrix,paste0(output_file_prefix,number_of_clones,"clones.proportions.solution.",number_of_answers,"of",top_number,".txt"),sep="\t", row.names=FALSE)
      write.table(presence_matrix,paste0(output_file_prefix,number_of_clones,"clones.mutation.matrix.",number_of_answers,"of",top_number,".txt"),sep="\t", row.names=FALSE)
      #write.table(gauchoOutput@best,paste0(output_file_prefix,number_of_clones,"clones.highest.score.solution.",number_of_answers,"of",top_number,".txt"),sep="\t", row.names=FALSE,col.names=FALSE)
    }
  
    ## Output fitness convergence plot
    if(outType %in% c("complete","fitness")){
      message("Outputting fitness convergence plot")
      if(outType=="complete"){ png(paste0(output_file_prefix,number_of_clones,"clones.fitnessconvergence.png")) }
      plot(gauchoOutput,ylim=yRange)
      if(outType=="complete"){ dev.off() }
    }
      
    ## Produce heatmap
    if(outType %in% c("complete","heatmap")){
      message("Outputting heatmap")
      ## Iterate through presence_matrix data and replace values with hex colours
      clones=presence_matrix
      colnames(clones)=1:ncol(clones) 
      for(i in 1:ncol(clones)){
        ## Use the modulo remainder NOT i, so we can recycle the same colour object
        ## regardless of the dimensions of presence_matrix
        recycler=ifelse(i%%length(colours)==0,length(colours),i%%length(colours))
        clones[clones[,i]=="1",i]=colours[recycler]
        clones[clones[,i]=="0",i]=secondary[4] # background is ICR light grey
      }
      if(outType=="complete"){ png(paste0(output_file_prefix,number_of_clones,"clones.solution",number_of_answers,"of",top_number,".heatmap.png")) }
      heatmap.plus(as.matrix(observation_matrix),Colv=NA,col=heatmapColours,scale="none",Rowv=NULL,labRow=gauchoOutput@names,RowSideColors=clones,cexRow=0.9)
      if(outType=="complete"){ dev.off() }
    }  
    
    ## Produce proportion stacked barplot
    if(outType %in% c("complete","proportion")){
      message("Outputting proportion barplot")
      rownames(proportion_matrix)=colnames(observation_matrix)
      if(outType=="complete"){ png(paste0(output_file_prefix,number_of_clones,"clones.solution",number_of_answers,"of",top_number,".proportions.png")) }
      ## Produce a plot with TWO areas. The top (1/10th) contains the legend, the bottom (9/10ths) contains the barplot
      layout(rbind(1,2), heights=c(1,9))
      old.par <- par(mar = c(0, 0, 0, 0)) # Save margins before we change them
      plot.new() # A blank canvas to hold the legend
      legend("center",LETTERS[1:pseudo_number_of_clones], pch = 15,
             col=colours,ncol=pseudo_number_of_clones,bty ="n",cex=1.5)
      par(old.par) # Put the old margins back
      barplot(t(proportion_matrix),col=colours,border=NA,xlab="Sample name",ylab="Proportion of cells")
      layout(1) # Restore layout
      if(outType=="complete"){ dev.off() }
    }
    
    ## Produce phylogeny plot
    if(outType %in% c("complete","phylogeny")){
      message("Outputting phylogeny plot")
          
      ## The phylogeny matrices are not compliant with "graph", so we edit them:
      ## 1.) We remove the bottom row (which is the original encoded string from the solution)
      pm=phylogeny_matrix[1:ncol(phylogeny_matrix),]
      ## 2.) Matrices must have 0 diag, sum(column)<=1, keeping the leftmost value
      diag(pm)=0
      for(i in 1:nrow(pm)){
        for(j in ncol(pm):1){
          if(sum(pm[,j])>1){
            pm[i,j]=0
          }
        }
      }
      
      ## Set lineage names and cast to graphNEL object
      rownames(pm)=LETTERS[1:nrow(pm)]
      colnames(pm)=rownames(pm)
      pg=as(pm,"graphNEL")
      
      ## Set up plotting parameters
      nAttrs = list() ## per node attributes
      nAttrs$fillcolor = rep(NA,nrow(pm))
      for(i in 1:nrow(pm)){
        ## Use the modulo remainder NOT i, so we can recycle the same colour object
        recycler=ifelse(i%%length(colours)==0,length(colours),i%%length(colours))
        nAttrs$fillcolor[i]=colours[recycler]
      }
      names(nAttrs$fillcolor) = LETTERS[1:nrow(pm)]
      nAttrs$color = rep(NA,nrow(pm))
      names(nAttrs$color) = LETTERS[1:nrow(pm)] ## Assign NA color to shape to ensure clean fill
      
      ## We store the original margins before plotting as stupid RGraphViz screws them up
      old.par <- par(mar = c(0, 0, 0, 0))
      if(outType=="complete"){ png(paste0(output_file_prefix,number_of_clones,"clones.solution",number_of_answers,"of",top_number,".phylogeny.png")) }
      plot(pg, nodeAttrs = nAttrs, attrs=list(edge = list(color=secondary[4],lwd=4))) # "attrs" sets global attributes
      if(outType=="complete"){ dev.off() }
      par(old.par)
    }
    
    #####################################
    ## End plotting / producing output ##
    #####################################
    
  }  
  
  ######################################
  ## End looping over TOP 5 solutions ##
  ######################################
  
}
# /end of gauchoReport



