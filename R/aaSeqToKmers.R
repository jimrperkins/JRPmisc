aaSeqToKmers <- function(seqs, k=5, allPoss=FALSE, rev=FALSE) {
    # First step - make "bag of words/dictionary of all kmers
    if(allPoss == FALSE) {
        kmerAllSeqs <- get.kmers(seqs, .k=k)
        bag <- sort(kmerAllSeqs$Kmers)
    } else if (allPoss == TRUE) {
        allLetters <- sort(unique(strsplit(paste(seqs, collapse = ""), "")[[1]]))
        gridBag <- expand.grid(rep(list(allLetters), k), stringsAsFactors = FALSE)
        # Use C implementation of paste(collapse="") as orders of magnitude quicker
        bag <- as.character(apply(gridBag, 1, str_c, collapse=""))
    }
    # Second step - use lapply to get all kmers for each sequence as a list
    kFreqL <- lapply(as.character(seqs), get.kmers, .k=k)

    # Third step - match kmers to the bag and make a vector
    kmerV <- lapply(kFreqL, function(x){
        countV <- vector(length=length(bag))
        countV[match(x$Kmers, bag)] <- x$Count
        return(countV)
        })
    kmerM <- do.call("rbind", kmerV)
    return(kmerM)
}
