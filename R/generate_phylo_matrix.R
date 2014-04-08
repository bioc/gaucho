## Creates a matrix where rows indicate the clones in which a mutation must occur.
## Identify the root nodes, set these to have 1 in the diagonal corresponding to their column
## For each root node, find their children, copy the parents column and add 1 to their diagonal
## find the children of those children and repeat until matrix constructed
generate_phylo_matrix<-function(input,number_of_clones) {
  pmat<-matrix(rep(0, (number_of_clones^2)), ncol=number_of_clones, byrow=TRUE)
  next_generation<-which(input==0)
  new_next_generation<-c()  
  for (node in next_generation) {
    pmat[node,node]<-1
    new_next_generation<-c(node, new_next_generation) 
  }
  next_generation<-new_next_generation
  while (length(next_generation>0)) {
    new_next_generation<-c() 
    for (i in next_generation) {
      getme<-which(input==i)
      for (j in getme) {
        pmat[,j]<-pmat[,i]
        pmat[j,j]<-1
        new_next_generation<-c(j, new_next_generation) 
      }
    }
    next_generation<-new_next_generation
  }
  return(pmat)
}
## /end of generate_phylo_matrix