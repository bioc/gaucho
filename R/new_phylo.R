## Generates phylogenies with one or more root nodes
new_phylo<-function(nroot,number_of_clones) {
  if (nroot > number_of_clones) {stop("Terminating Script: Number of roots cannot be greater than number of clones")}
  inp<-vector(length=number_of_clones)
  
  ## Maximum number of roots is the number of clones
  nroot=min(nroot,number_of_clones) 
  
  ## When nroot is 0, we randomise the number of roots 
  if (nroot==0) { 
    i<-round(runif(1,0.5, number_of_clones+0.49)) 
    inp[1:i]<-0
    last_gen_end<-i
    spaces_left<-number_of_clones-i
  } else { 
    inp[1:nroot]<-0
    last_gen_end<-nroot
    spaces_left<-number_of_clones-nroot
  }
  
  last_gen_start<-1
   
  while (spaces_left > 0) {
    length_next_gen<-round(runif(1,0.5, spaces_left+0.49))
    current_sample<-sample(rep(last_gen_start:last_gen_end, number_of_clones))
    inp[(last_gen_end+1):(last_gen_end+length_next_gen)]<-current_sample[1:length_next_gen]
    last_gen_start<-last_gen_end+1
    last_gen_end<-last_gen_end+length_next_gen
    spaces_left<-spaces_left-length_next_gen
  }

  return(inp)
}
## /end of new_phylo