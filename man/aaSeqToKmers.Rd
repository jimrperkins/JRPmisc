\name{aaSeqToKmers}
\alias{aaSeqToKmers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to convert character strings into a kmer instance matrix
}
\description{
Takes as input a vector of character strings (such as amino acid seqs) and convert it to a matrix of kmer instances.
First, the function finds all kmers in any of the the strings, (bag of words) and uses this for the column of the matrix, with each row representing an element of the original vector (i.e. an a.a. sequence)
}
\usage{
aaSeqToKmers(seqs, k = 5, allPoss = FALSE, rev = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqs}{
  Vector of character strings.
}
  \item{k}{
  Size of the kmers found.
}
  \item{allPoss}{
  Do we use all potential kmers as columns, even those that don't exist in the input sequences/strings? Warning - these will be number of letters^k in size.
}
  \item{rev}{
  Look at the reverse sequences as well.
}
}
\details{
We are essentially wrapping \code{\link[tcR]{get.kmers}} to produce a big matrix of kmer instances.
}
\value{
A matrix of kmer instances. All kmers in any of the the strings as columns, with each row representing an element of the original vector (i.e. an input a.a. sequence).
}
\references{
 To update
}
\author{
James R Perkins \email{jimrperkins@gmail.com}
}
\note{
  If we use the allPoss, we will potentially generate VERY large matrices and it will take a long time to run.
}

\seealso{
    \code{\link[tcR]{get.kmers}}
}
\examples{
# Get all 5-mers in TNF
aaTNF <- c("MSTESMIRDVELAEEALPKKTGGPQGSRRCLFLSLFSFLIVAGATTLFCLLHFGVIGPQREEFPRDLSLISPLAQAVRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL")
rep(aaTNF, 2)
ee <- aaSeqToKmers(rep(aaTNF, 2), k=5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }% use one of  RShowDoc("KEYWORDS")
