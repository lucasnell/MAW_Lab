#!/usr/bin/env python3

"""Script to convert Princeton site's probability weight matrices (PWMs) to MEME files."""


import pandas as pd
import argparse as ap
import os

__author__ = 'Lucas Nell'


memeHead = '''MEME version 4

ALPHABET= ACGT

strands: + -
'''
"""Header at the top of each MEME file."""

motBlock = '''
MOTIF %(mot)s
letter-probability matrix:
%(mat)s
'''
"""Block of string for each motif."""


pyFolder = os.path.dirname(os.path.abspath(__file__))
"""Folder this file is in. This is the default directory to work from."""



# -----------------------------
# Setting up parser (what allows this script to take arguments)
# -----------------------------

ScriptDescript = '''Create MEME file from probability weight matrix.'''

Parser = ap.ArgumentParser(description = ScriptDescript)
Parser.add_argument('files', type = str, metavar = 'F', nargs = '+',
                    help = 'Probability weight matrix file name(s).')
Parser.add_argument('-d', '--directory',  type = str, metavar = 'D', required = False,
                    help = 'Full path to PWM files. Default is folder this file is in.', 
                    default = pyFolder)

# Dictionary of arguments supplied to script.
args = vars(Parser.parse_args())

# Extracting info from args dictionary.
pwmFiles = args['files']
fDir = args['directory']

# Change working directory to default or the supplied directory.
os.chdir(os.path.expanduser(fDir))


# Function to make a MEME file from a PWM
def makeMEME(pwmFile, motifName = None):
    if motifName is None:
        motName = pwmFile.replace('_pwm', '').replace('.txt', '')
    else:
        motName = motifName
    # Only keeping the probabilities, in order A, C, G, T
    pwm = pd.read_table(pwmFile, delim_whitespace = True, nrows = 4,
                        skiprows = 1, header = None, 
                        encoding = 'utf-8').drop(0, axis=1).T
    # Turning to string, removing leading spaces on all lines
    pwdTabStr = pwm.to_string(header = False, index = False)
    # Final string with everything together
    memeStr = str(memeHead + motBlock % dict(mat = pwdTabStr, mot = motName))
    with open(str('%s_MEME.txt' % motName), mode = 'wt') as f:
        f.write(memeStr)
    print('Wrote to %s_MEME.txt' % motName)
    return


# Doing function on all supplies files
for p in pwmFiles:
    makeMEME(str(p))






