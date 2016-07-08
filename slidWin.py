#!/usr/local/apps/anaconda/3-2.2.0/bin/python

"""Sliding Window function on raw read depth files from alignments to mouse chr Y or X.

Note: No normalization based on chr19 reads is done here. That is saved for R scripts.
"""

import numpy as np
import pandas as pd
import os
import argparse
from time import strftime
from multiprocessing import Pool

__author__ = 'Lucas Nell'


###############################
#  Arguments needed...
###############################
# Directory to find the raw coverage file
# Chromosome aligned to
# Window length (optional, default = 1000)
# Increment by which to slide window (optional, default = 500)
# Filename(s) of raw depths (all columns should have header_i == sample name_i and
#                            values should be simply read depth at that position)
# Should it be run with each file on a separate thread? (optional)




# ========================================================================================
# ========================================================================================
# ========================================================================================
# ========================================================================================

#               READING INPUT ARGUMENTS, CREATING OBJECTS

# ========================================================================================
# ========================================================================================
# ========================================================================================
# ========================================================================================


# -----------------------------
# Setting up parser
# -----------------------------
ScriptDescript = '''Sliding Window function on raw read depth files'''

Parser = argparse.ArgumentParser(description = ScriptDescript)
Parser.add_argument('-d', '--directory',  type = str, metavar = 'D',
                    help = 'Full path to raw-depth file', required = True)
Parser.add_argument('-c', '--chromosome', type = str, metavar = 'C',
                    choices = ['X', 'Y', 'Xb9'], required = True,
                    help = 'The chromosome this raw-depth file was aligned to.' + \
                            'Choices are X, Y, or Xb9. The latter is chrX from ' + \
                            'mouse genome build mm9.')
Parser.add_argument('-w', '--window', type = int, metavar = 'W',
                    required = False, default = 1000,
                    help = 'Length of the sliding window')
Parser.add_argument('-i', '--increment', type = int, metavar = 'I',
                    required = False, default = 500,
                    help = 'Increment by which to slide window')
Parser.add_argument('files', type = str, metavar = 'F', nargs = '+',
                    help = 'Raw depth (ending in ".txt.gz") file(s)')
Parser.add_argument('-t', '--threads', type = int, metavar = 'T',
                    required = False, default = 1,
                    help = 'How many threads can this use? The default is 1 ' +\
                           '(not useful to have more threads than files).')



# -----------------------------
# Now reading the arguments
# -----------------------------
args = vars(Parser.parse_args())


# From this command line input...
# $py -i 500 bam.bam foo.bam ffoo.bam -d /lustre1/lan/MouseGenome -c Xb9 -w 1000
# print(args) gives you...
# {'directory': '/lustre1/lan/MouseGenome', 'chromosome': 'Xb9', 'parallel': False,
#  'increment': 500, 'files': ['bam.bam', 'foo.bam', 'ffoo.bam'], 'window': 1000}


# -----------------------------
#   Making objects
# -----------------------------

Dir = args['directory']
if Dir[(len(Dir)-1)] != '/':
    Dir += '/'
Chr = args['chromosome']
Threads = args['threads']
Window = args['window']
Incr = args['increment']
Files = args['files']
del args  # No longer needed

# Total length of chromosome
Length = {'X': 171031299, 'Y': 89532639, 'Xb9': 166650296}[ Chr ]

# Length of each chunk DataFrame
ChunkLen = np.ceil(Length / Incr)

# Number of chunks per Window
chuWin = int(Window / Incr)



# ========================================================================================
# ========================================================================================
# ========================================================================================
# ========================================================================================

#               FUNCTIONS

# ========================================================================================
# ========================================================================================
# ========================================================================================
# ========================================================================================

# -----------------------------
#       Raw depths --> Chunks
# -----------------------------

# Function to get chunk averages from a pandas 'chunked' array of raw depths
def GetChunks(RawDepths, Names):

    """Get chunk averages from a pandas 'chunked' array of raw depths."""

    ChunkAvgs = np.empty([ChunkLen, len(Names)])
    # For loop down the array
    for i, chunk in enumerate(RawDepths):
        ChunkAvgs[i] = [ chunk[ name ].mean() for name in Names ]
    return ChunkAvgs


# Function to read raw depths and create a gzipped txt file of chunk means
# Outputs filename for next step
def AllTheChunks(File):

    """Function to read raw depths and create a gzipped txt file of chunk means.

    Outputs filename for next step."""

    RawDepthTab = pd.read_table(File, compression = 'gzip',
                                delim_whitespace = True, chunksize = Incr)
    ColNames = list(RawDepthTab.get_chunk(1))
    ChunkAvgs = GetChunks(RawDepthTab, ColNames)
    TxtFile = str(Incr) + 'chunk' + '_' + File.replace('.gz', '')
    # 'Code' to make a float of necessary # decimals
    fCode = str('%-0.' + str(np.ceil(np.log10(Window))) + 'f')
    np.savetxt(TxtFile, ChunkAvgs, fmt = fCode, delimiter = "\t",
           header = "\t".join(ColNames), comments = "")
    os.system('gzip -f ' + TxtFile)
    return str(str(Incr) + 'chunk' + '_' + File)



# -----------------------------
#       Chunks  -->  Windows
# -----------------------------


def CleanEnd(RollingChunks, ChunkAvgs, ind):
    """Adjusts last rolling-chunk average for uneven last-chunk length.

    Args:
        RollingChunks (pd.DataFrame): Table of rolling chunk averages.
        ChunkAvgs (pd.DataFrame): Table of chunk-mean averages.
        ind (int): Index of column you're fixing.

    Returns:
        (float): New rolling-chunk average.
    """
    Partial = Length % Incr  # Size of last chunk (!=Incr)
    Avg0 = RollingChunks.iloc[-1, ind]
    End = ChunkAvgs.iloc[-1, ind]
    # Sum of all chunk avgs other than end one
    NonEnd = (Avg0 * chuWin) - End
    # New sum of chunk averages, where End is underweighted
    Total = ( End * (Partial / Incr) ) + NonEnd
    # New chuWin, where End is underweighted
    EndN = (chuWin - 1) + (Partial / Incr)
    return Total / EndN



def SlidingWin(InFile):
    """Create, save table of read-depth averages across a sliding window.

    Args:
        InFile (str): Location of tab-delimited table of chunk averages.

    Returns:
        Nothing.
        Writes output to file.
    """
    # Input dataframe of chunk averages
    ChunkAvgs = pd.read_table(InFile)
    ColNames = list(ChunkAvgs)
    # Calculates a rolling mean of length 'chuWin', then shifts it so it starts with
    # numbers, then it removes the NaN at the end
    RollingChunks = pd.rolling_mean(ChunkAvgs, chuWin).dropna().reset_index(drop = True)
    for i in range(len(ColNames)):
        RollingChunks.iloc[-1, i] = CleanEnd(RollingChunks, ChunkAvgs, i)
    OutFile = '../Win/' + InFile.replace(str(Incr) + 'chunk', 'Win')
    RollingChunks.to_csv(OutFile, index = False, sep = '\t', compression = 'gzip',
                         float_format = '%-0.5f')
    print('Written to', OutFile)
    return



# -----------------------------
#       Putting steps together
# -----------------------------

# Messages to estimate time usage
def StartMessage(ID):
    print('~'*25)
    print(' * SlidWin.py on', ID)
    print(' Start Time =', strftime("%H:%M:%S"))
    print('\n', ' .' * 3)
def EndMessage(ID):
    print('\n', ' .' * 3, '\n')
    print(' * SlidWin.py on', ID)
    print(' End Time =', strftime("%H:%M:%S"))
    print('~'*25)


# Give it a file name, it'll run AllTheChunks and SlidingWin
def ForTheWin(File):
    StartMessage(File)
    NewFile = AllTheChunks(File)
    print('Now sliding a window across %s chunks . . .' % File)
    SlidingWin(NewFile)
    EndMessage(File)




# ========================================================================================
# ========================================================================================
# ========================================================================================
# ========================================================================================

#               RUNNING PROGRAM ON INPUT DATA

# ========================================================================================
# ========================================================================================
# ========================================================================================
# ========================================================================================



# Now doing it for all filenames given, parallel or serial
if __name__ == '__main__':
    os.chdir(Dir)
    if Threads > 1:
        with Pool(processes = Threads) as pool:
            pool.map(ForTheWin, Files)
    else:
        for f in Files:
            ForTheWin(f)


# os.system('''scp /Users/lucasnell/uga/Python/MouseCoverage/slidWin.py \
# lan@xfer2.gacrc.uga.edu:~/CovTools''')





# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # For testing...
# args = {'directory': '/lustre1/lan/MouseGenome', 'chromosome': 'X9', 'parallel': False,
#         'increment': 500, 'files': ['bam.txt.gz', 'foo.txt.gz', 'ffoo.txt.gz'],
#         'window': 1000}
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
