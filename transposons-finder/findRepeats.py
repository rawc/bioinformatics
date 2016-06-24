#!/usr/bin/env python3
# Name: Louis Chatta (lchatta)
# Group Members: (lchatta)

"""
This program finds all Inverted Repeats and Direct Repeats 
in a sequnce. Allowing the user to find all possible DNA transposons
and Retrotransposons. 

Input: 
FASTA File

Output:
All found IRs and DRs in a sequnce to STDOUT

"""

class Util():

    def getInvertedRepeat(self,string):
        return self.getComplementStrand(string[::-1])

    def getComplementStrand(self,bases):
        '''
        A nice function that returns the reverse strand 
        '''
        reverseBase = []
        for base in bases:
            if base == 'A':
                reverseBase.append('T')
            elif base == 'T':
                reverseBase.append('A')
            elif base == 'C':
                reverseBase.append('G')
            elif base == 'G':
                reverseBase.append('C')
            else:
                reverseBase.append(base.upper())
        return reverseBase
      

    def getReverseStrand(bases):
        '''
        A nice function that returns the reverse strand 
        '''
        reverseBase = []
        for base in bases:
            if base == 'A':
                reverseBase.append('T')
            if base == 'T':
                reverseBase.append('A')
            if base == 'C':
                reverseBase.append('G')
            if base == 'G':
                reverseBase.append('C')
        return reverseBase[::-1]
#######

########################################################################
# CommandLine
########################################################################
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - Find all Inverted and Direct repeats', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >'
                                             )
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

# Used to find all repeating elements
# Will find all of the Direct Repeats and Inverted Repeats
class RepeatFinder():

    def __init__(self, seq, header,revSeq=False):
        self.seq = seq
        self.revSeq = revSeq
        self.header = header

        self.subStrings = []
        self.util = Util()

# generates all substrings
    def getAllSubstrings(self,string):
      length = len(string)
      alist = []
      for i in range(length):
        for j in range(i,length):
          if len(string[i:j]) <1 : continue
          newString = string[i:j + 1]
          if newString not in alist:
             alist.append((newString, i, j+1) )
      self.subStrings = alist
      return alist

# generates all repeating elements
    def findAllRepeats(self):
        seenRepeats = []
        for subString, startPos, endPos in self.subStrings:
            foundRepeat = ([item for item in seenRepeats if item[0] == subString])
            if not foundRepeat:
                seenRepeats.append((subString,startPos,endPos))
            else:
                print('    ', self.buildSequence((foundRepeat[0])[0],subString,(foundRepeat[0])[1],endPos))
          
# generates all repeating inversion elements
    def findAllInversionRepeats(self):
        seenRepeats = []
        for subString, startPos, endPos in self.subStrings:
            invertedRepeat = self.util.getInvertedRepeat(subString)
            if subString == invertedRepeat: continue
            foundRepeat = ([item for item in seenRepeats if item[0] == invertedRepeat])

            if foundRepeat:
                print('    ', self.buildSequence((foundRepeat[0])[0],subString,(foundRepeat[0])[1],endPos))
            else:
                seenRepeats.append((subString,startPos,endPos))

# creates a sequncee based on a provided start and stop position 
    def buildSequence(self,startRepeat,endRepeat, startPos, endPos):

        startStrting = ''.join(startRepeat)
        endString = ''.join(endRepeat)

        innerSequence = ''.join((self.seq[startPos+len(startRepeat):endPos-len(endRepeat)]))

        return '{} {} {}'.format(startStrting, innerSequence ,endString)


# prints all of the transposons
    def printTransposons(self):
        self.getAllSubstrings(self.seq)
        print ('Direct Repeats')
        self.findAllRepeats()
        print('Inverted Repeats')
        self.findAllInversionRepeats()
        print(' ')

########################################################################
# Main
# Here is the main program
# 
#
########################################################################
   
from sequenceAnalysis import NucParams,FastAreader,ProteinParam
def main(myCommandLine=None):
    '''
    Implements the Usage exception handler that can be raised from anywhere in process.  

    # '''

    myCommandLine = CommandLine(myCommandLine)
    reader = FastAreader(myCommandLine.args.inFile)

    # loop through all of the sequences in a file
    # and calculate the respective ORF and write them
    for head, seq in reader.readFasta():
        nucParams = NucParams('')
        nucParams.addSequence(seq)
        nucParams.buildNucleotide()
        bases = list((''.join(nucParams.codons)))
        print(' '.join(bases))
        reverseBases = Util.getReverseStrand(bases)
        
        finder = RepeatFinder(bases,head,revSeq=reverseBases)
        finder.printTransposons()

import sys

if __name__ == "__main__":
    main(sys.argv[1:]) # only send the command line args that are not the program name


