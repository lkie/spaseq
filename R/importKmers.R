#' Importing counts of kmers from structure of directories
#'
#' Function imports data generated from FASTQ files by bash scripts
#' and creates a standardized R list.
#'
#' @param dataDir Path to the main directory with the processed data (do not include "/" at the end)
#' @param indices Indexes to be imported. If not provided, all dirs from the dataDir will be used.
#' @param constructs Names of the constructs to be imported. 
#' If not provided, all files with ".txt.gz" at the proper location extension will be imported.
#' @param maxKmerLen The longest kmer length that will be imported (default = 8)
#' @return all_kmers Multi-level list structure. The first level are different constructs,
#' the second level is sequencing index of the sample, 
#' the third level is kmer length, 
#' which contains a matrix of counts of different kmers at different positions of the constructs
#' @author Lukasz Jan Kielpinski
#' @export importKmers

importKmers <- function(dataDir = ".", indices, constructs, maxKmerLen = 8){
  if(missing(indices)){indices <- list.dirs(dataDir, recursive = FALSE, full.names = FALSE)}
  if(missing(constructs)){
    constructs <- c()
    for(index in indices){
      local_files <- list.files(paste(dataDir, "/", index, sep = ""))
      matching_files <- local_files[grep(local_files, pattern = ".txt.gz")]
      constructs <- c(constructs, substr(matching_files, start = 1, stop = nchar(matching_files) - 7)) # "- 7" to strip the construct name from ".txt.gz"
    }
    constructs <- unique(constructs)
  }
  
  all_kmers <- list()
  for(construct in constructs){
    kmers <- list()
    
    #Find out how long is the construct's randomized index:
    for(index in indices){
      file_name <- paste(dataDir, "/", index, "/", construct, ".txt.gz", sep = "")
      if(file.exists(file_name)){
        insertLength <- nchar(scan(file = file_name, what = "character", n = 1, quiet = TRUE))
        break
      }
    }
    
    for(index in indices){
      file_name <- paste(dataDir, "/", index, "/", construct, ".txt.gz", sep = "")
      if(file.exists(file_name)){
        kmers[[index]] <- list()
        for(kmerLen in 1:maxKmerLen){
          kmers[[index]][[kmerLen]] <- matrix(ncol = insertLength - kmerLen + 1, nrow = 4**kmerLen)
          rownames(kmers[[index]][[kmerLen]]) <- makeAllKmers(kmerLen)
          for(pos in 1:(insertLength - kmerLen + 1)){
            kmers_summary <- read.table(paste(dataDir, "/", index, "/", kmerLen,"/",pos, ".", construct, ".txt.gz", sep = ""))
            where_to_add <- match(x = kmers_summary[,2], table = rownames(kmers[[index]][[kmerLen]]))
            kmers[[index]][[kmerLen]][where_to_add, pos] <- kmers_summary[,1]
          }
          kmers[[index]][[kmerLen]][is.na(kmers[[index]][[kmerLen]])] <- 0
        }
      }
    }
    all_kmers[[construct]] <- kmers
  }
  all_kmers
}

