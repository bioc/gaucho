\name{gaucho}
\alias{gaucho}
\title{Genetic Algorithm for Understanding Clonal Heterogeneity and Ordering (GAUCHO)}
\usage{
gaucho(observations, number_of_clones, pop_size = 100, mutation_rate = 0.8,
  iterations = 1000, stoppingCriteria = round(iterations/5),
  parthenogenesis = 2, nroot = 0, contamination = 0,
  check_validity = TRUE)
}
\arguments{
  \item{observations}{Observation data frame where each row
  represents an SNV and each column represents a discrete
  sample separated by time or space.  Note that the data
  frame must have column names and row names. Every value
  must be a proportion between 0 and 1.  See details}

  \item{number_of_clones}{An integer number of clones to be
  considered}

  \item{pop_size}{The number of individuals in each
  generation}

  \item{mutation_rate}{The likelihood of each individual
  undergoing mutation per generation}

  \item{iterations}{The maximum number of generations to
  run}

  \item{stoppingCriteria}{The number of consecutive
  generations without improvement that will stop the
  algorithm. Default value is 20\% of iterations.}

  \item{parthenogenesis}{The number of best-fitness
  individuals allowed to survive each generation}

  \item{nroot}{Number of roots the phylogeny is expected to
  have.When nroot=0, a random integer between 1 and the
  number of clones is generated for each phylogeny}

  \item{contamination}{Is the input contaminated?  If set
  to 1, an extra clone is created in which to place
  inferred contaminants}

  \item{check_validity}{Unless set to false, eliminate any
  clones with no new mutations, disallow those clones.
  Increases computational overheads.}
}
\value{
Returns an object of class ga \code{\link{ga-class}}.  Note
that the number of clones and number of cases are stored in
the unused min and max slots of the output object.
}
\description{
Use a genetic algorithm to find the relationships between
the values in an input file - the package was written to
deal with single nucleotide variants (SNVs) in mixtures of
cancer cells, but it will work with any mixture.  It will
calculate appropriate phylogenetic relationships between
clones them and the proportion of each clone that each
sample is composed of.  For detailed usage, please read the
accompanying vignette.
}
\details{
The input data should be a data.frame containing
proportions of cells that contain a feature.  There are a
number of ways to create these data, including merging the
balance of alleles and copy number of an SNV using the
equation min(1,r*CN/(r+R)), where CN is the copy number, r
is the number of non-reference reads and R is the number of
reference reads.  For example, if a site were sequenced to
a depth of 100x, with 25 non-reference reads and 75
reference reads and diploid copy number, the result would
be min(1,25*2/(25+75)) = 0.5.  Therefore, 50\% of the cells
in the sample contain the SNV. Further details are
available in the accompanying vignette.
}
\examples{
## The vignette provides far more in-depth explanation and examples ##

## Load the included simple example data
gaucho_simple_data = read.table(file.path(system.file("extdata",package="gaucho"),"gaucho_simple_data.txt"),header=TRUE,row.names=1)

## Run gaucho using 3 clones and a phylogeny with a single root
solution=gaucho(gaucho_simple_data, number_of_clones=3,nroot=1,iterations=1000)

## Create the four output plots
gauchoReport(gaucho_simple_data,solution,outType="fitness")
gauchoReport(gaucho_simple_data,solution,outType="heatmap")
gauchoReport(gaucho_simple_data,solution,outType="phylogeny")
gauchoReport(gaucho_simple_data,solution,outType="proportion")

## Output the solution and plots in the current working directory
# gauchoReport(gaucho_simple_data,solution)
}
\author{
Alex Murison \email{Alexander.Murison@icr.ac.uk} and
Christopher Wardell \email{Christopher.Wardell@icr.ac.uk}
}
\seealso{
\code{\link{ga-class}}, \code{\link{ga}},
\code{\link{gauchoReport}},
\code{\link{gaucho_simple_data}},
\code{\link{gaucho_hidden_data}},
\code{\link{gaucho_synth_data}},
\code{\link{gaucho_synth_data_jittered}},
\code{\link{BYB1_G07_pruned}}
}

