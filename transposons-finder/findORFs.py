#!/usr/bin/env python3
# Name: Louis Chatta (lchatta)
# Group Members: (lchatta, cjavery)

"""
This program finds the open reading frame of sequences:
The Extra credit has been done through the command line 

Here is some psuedocode for finding the ORF:
    loop through the sequence
	    extract codon
	    if the codon is a start codon:
	        add it to the list
	  	
	  	if the codon is a stop codon and we have seen a start previously:
	        save orf

	after loop, sort all the orfs by length and save to file 

"""
outFile = None


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
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        self.parser.add_argument('outFile', action = 'store', help='output file name') 
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= range(0, 2000), action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', nargs='?', help='start Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

class RepeatFinder():
    frameList = [] # class var that holds all of the frames to print

    def __init__(self, seq, header,minGene=0,start=["ATG"],longestGene=False,revSeq=False):
        self.seq = seq
        self.revSeq = revSeq
        self.header = header
        self.minGene = minGene
        if start:
            self.start = start
        else: 
            self.start = ["ATG"]
        if minGene:
            self.minGene = minGene
        else:
            self.minGene = 0
        self.longestGene = longestGene
        self.frameList = []
        self.util = Util()


    def get_all_substrings(self,string):
      length = len(string)
      alist = []
      for i in range(length):
        for j in range(i,length):
          if len(string[i:j]) >7  or len(string[i:j]) <1 : continue
          # self.orderedSubSequences[string[i:j +1]] = i
          newString = string[i:j + 1]
          if newString not in alist:
             alist.append((newString, i, j+1) )
      self.subStrings = alist
      return alist

    def findAllRepeats(self):
        seenRepeats = []
        for subString, startPos, endPos in self.subStrings:
            foundRepeat = ([item for item in seenRepeats if item[0] == subString])
            if not foundRepeat:
                seenRepeats.append((subString,startPos,endPos))
            else:
                print('    ', self.buildSequence((foundRepeat[0])[0],subString,(foundRepeat[0])[1],endPos))
          

    def findAllInversionRepeats(self):
        seenRepeats = []
        for subString, startPos, endPos in self.subStrings:
            invertedRepeat = self.util.getInvertedRepeat(subString)
            if subString == invertedRepeat: continue
            foundRepeat = ([item for item in seenRepeats if item[0] == invertedRepeat])

            if foundRepeat:
                print('    ', self.buildSequence((foundRepeat[0])[0],subString,(foundRepeat[0])[1],endPos))
                # print('    ',foundRepeat, subString,endPos)
            else:
                seenRepeats.append((subString,startPos,endPos))


    def buildSequence(self,startRepeat,endRepeat, startPos, endPos):

        startStrting = ''.join(startRepeat)
        endString = ''.join(endRepeat)

        innerSequence = ''.join((self.seq[startPos+len(startRepeat):endPos-len(endRepeat)]))

        return '{} {} {}'.format(startStrting, innerSequence ,endString)


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

    '''

    if myCommandLine is None:

        myCommandLine = CommandLine([ 'tass2.fa',
                                      'tass2ORFdata-ATG-100.txt',
                                      '--longestGene',
                                      '--start=ATG',
                                      '--minGene=100'])
    else :

        myCommandLine = CommandLine(myCommandLine)

    myCommandLine.args.inFile #has #the input file name
    outFile = myCommandLine.args.outFile  #the output file name
    myCommandLine.args.longestGene #is True if only the longest Gene is desired
    myCommandLine.args.start #is a list of start codons
    myCommandLine.args.minGene #is the minimum Gene length to include
    

    # Clear the file if created previously, 
    # and open it 
    orfReader = FastAreader(myCommandLine.args.inFile)
    open(myCommandLine.args.outFile, 'w').close()
    f = open(outFile, 'a')

    # loop through all of the sequences in a file
    # and calculate the respective ORF and write them
    for head, seq in orfReader.readFasta():
        nucParams = NucParams('')
        nucParams.addSequence(seq)
        nucParams.buildNucleotide()
        bases = list((''.join(nucParams.codons)))
        print(' '.join(bases))
        reverseBases = getReverseStrand(bases)
        
        finder = RepeatFinder(bases,head,minGene=myCommandLine.args.minGene, longestGene=myCommandLine.args.longestGene, start=myCommandLine.args.start,revSeq=reverseBases)
        finder.get_all_substrings(finder.seq)
        print ('Direct Repeats')
        finder.findAllRepeats()
        print('Inverted Repeats')
        finder.findAllInversionRepeats()

        print(' ')


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

import sys

if __name__ == "__main__":
    main(sys.argv[1:]) # only send the command line args that are not the program name


