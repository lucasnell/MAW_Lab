#!/usr/bin/env Rscript


# Getting directory this script is in
allArgs <- commandArgs(trailingOnly = FALSE)
scriptPath <- gsub('--file=', '', allArgs[grep('--file', allArgs)])
separatedPath <- unlist(strsplit(scriptPath, '/'))
scriptDir <- paste0(head(separatedPath, -1), collapse = '/')
# Loading table
tableName <- paste(scriptDir, 'convert_to_Glazer_table.RData', sep = '/')
load(tableName)


# Getting filename passed to script
filename <- commandArgs(trailingOnly = TRUE)

# colnames should be c('Chromosome', 'Position') [order is important, not actual names]

inTable <- read.table(filename, sep = '\t', colClasses = c('character', 'integer'), 
                      header = TRUE)


convertCoordinate <- function(chr, pos, direction = 'old2new'){
  
  translate <- function(pos, startA, startB, endB, orientation){
    # Translates from one coordinate system to a second
    if (!orientation == 'reverse'){
      pos2 <- startB + (pos - startA)
    } else {
      pos2  <- endB - (pos - startA)
    }
    return(pos2)
  }
  
  if (direction == 'old2new'){
    # Pull out right scaffold
    x <- scafTable[scafTable$OldChr==chr & scafTable$OldStart<=pos &
                     scafTable$OldEnd>=pos,]
    # Make sure there's exactly 1 scaffold meeting criteria
    if (!nrow(x) == 1){
      return (c(NA, NA))
    } else {
      # Calculate new coordinate
      newChr <- x[1,'NewChr']
      newPos <- translate(pos, x[1,'OldStart'], x[1,'NewStart'], x[1,'NewEnd'],
                          x[1,'NewOrientation'])
      return(c(newChr, newPos))
    }
  } else if (direction == 'new2old'){
    # Pull out right scaffold
    x <- scafTable[scafTable$NewChr == chr & scafTable$NewStart <= pos & 
                     scafTable$NewEnd >= pos,]
    # Make sure there's exactly 1 scaffold meeting criteria
    if (!nrow(x) == 1){
      return (c(NA, NA))
    } else {
      # Calculate new coordinate
      oldChr <- x[1,'OldChr']
      oldPos <- translate(pos,x[1,'NewStart'],x[1,'OldStart'],x[1,'OldEnd'],
                          x[1,'NewOrientation'])
      return(c(oldChr, oldPos))
    }
  }
  else{
    print('Direction must be old2new or new2old')
    return(c(NA, NA))
  }
}

outMatrix <- t(mapply(convertCoordinate, chr = inTable[,1], 
                      pos = inTable[,2], USE.NAMES = FALSE))

outDF <- as.data.frame(outMatrix)
colnames(outDF) <- colnames(inTable)
outDF[,2] <- as.integer(outDF[,2])

outFileName <- gsub('.txt', '_NEW.txt', filename)

if (outFileName == filename){
  stop('File must end in ".txt".')
}


write.table(outDF, outFileName, sep = '\t', 
            row.names = FALSE, quote = FALSE, na = '')




# # For testing...
# inTable <- data.frame(Chromosome = c('3', 'Un', 'Un', '1', '1'),
#                       Position = c(1538202, 37499024, 23343225, 541084, 540083))
# # Should output, in order
# #   same position
# #   now on chr 1
# #   now on chr 2
# #   different location on chr 1
# #   in between contigs in original assembly, so NA




# ============================================================================
# ============================================================================

#     Original commenting inside function

# ============================================================================
# ============================================================================

# This R function converts between the 'old' and 'new' stickleback assembly 
#   coordinate systems. The 'old' coordinate system is the assembly described in the 
#   Jones et al 2012 stickleback genome paper. It requires access to the FileS4 
#   NewScaffoldOrder.csv file.
# It has 4 inputs: chr, pos, direction, and scafFile. It returns a list of 
#   [chromosome, position].

# Inputs:
# chr is a number or string (e.g. 1, '1', 'Un') of the starting chromosome.
# pos is a number of the starting position.
# direction is either 'old2new' or 'new2old'.
# scafFile gives the path and file name to the file 'FileS4 NewScaffoldOrder.csv'

# Output:
# List of [chromosome, position] of the converted coordinate.

# Examples:
# convertCoordinate(3,1538202,'old2new','Path/to/FileS4 NewScaffoldOrder.csv')
# same position
# convertCoordinate('Un',37499024,'old2new','Path/to/FileS4 NewScaffoldOrder.csv') 
# now on chr 1
# convertCoordinate('Un',23343225,'old2new','Path/to/FileS4 NewScaffoldOrder.csv') 
# now on chr 2
# convertCoordinate("1",541084,'old2new','Path/to/FileS4 NewScaffoldOrder.csv') 
# different location on chr 1
# convertCoordinate('1',680442,'new2old','Path/to/FileS4 NewScaffoldOrder.csv') 
# reverse of previous line
# convertCoordinate(12,594205,'new2old','Path/to/FileS4 NewScaffoldOrder.csv') 
# used to be on Un
# convertCoordinate(1,540083,'old2new','Path/to/FileS4 NewScaffoldOrder.csv') 
# in between contigs in original assembly, so NA
