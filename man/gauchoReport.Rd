\name{gauchoReport}
\alias{gauchoReport}
\title{View solutions contained within gaucho output}
\usage{
gauchoReport(gauchoInput, gauchoOutput, outType = "complete",
  yRange = c(-250, 0), output_file_prefix = "")
}
\arguments{
  \item{gauchoInput}{Raw data analysed by gaucho()}

  \item{gauchoOutput}{Object of class ga produced by
  gaucho()}

  \item{outType}{Type of output desired - must be one of
  the following:
  "complete","fitness","heatmap","phylogeny","proportion"}

  \item{yRange}{Y-axis range when plotting fitness of
  individuals.  Default is c(-250,0)}

  \item{output_file_prefix}{Optional prefix for all output
  files}
}
\value{
Nothing is returned.
}
\description{
After running gaucho() on data, this function provides a
convenient way to view the solutions and also export them
as separate text files and images.  For detailed usage,
please read the accompanying vignette.
}
\details{
This method reports data for the fittest individual; in the
event of there being multiple individuals with identical
fitness, up to five individuals will be reported.  This
function's output is governed by the outType argument. All
options except for the default "complete" value result in
plotting the desired output to the current R session. When
outType=="complete", the following output is created for
each individual: the full length string, the phylogeny
matrix, the proportion matrix, the presence matrix, a
heatmap of the raw data with the assigned clones as
coloured bars at the side, a stacked barplot showing the
proportion of each clone at each timepoint and a plot
showing the phylogenetic relationship between the clones.
Note that the colours of the clones are consistent across
all plots and that the contamination clone (if present) is
always the last clone.  Also produced is a plot
illustrating the change in fitness as the generations
evolved.
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

