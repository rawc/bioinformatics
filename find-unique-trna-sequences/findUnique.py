#!/usr/bin/env python3
# Name: Louis Chatta (lchatta)
# Group Members: (lchatta, cjavery)

"""
This program finds all the unique trna subsequences
The Trna Class does all the heavy lifting
The main just loads the file data

INPUT:
FILENAME or STDIN

OUTPUT
FILENAME or STDOUT

"""
outFile = None


# Useful class that helps get substringsg
class util:
    def get_all_substrings(string):
        length = len(string)
        for i in range(length):
            for j in range(i + 1, length + 1):
                yield(string[i:j]) 

# class that holds all the trna data
class trnaFinder:
    def __init__(self,bases, head, allSequences):
        self.orderedSubSequences = {} # used for iterating through in an ordered fashion
        self.originalSequence = bases # the original sequence
        self.subSequences = (set(self.get_all_substrings(bases))) # the set off all subSequences
        self.header = head # the file header
        self.allSequences = (self.subSequences-allSequences) # all unique subSequences

        self.uniqueNonSuperSubstrings = self.removeAllSuperStrings(list(self.allSequences)) # the subSequences after removing super strings
        self.sortableHeader = head.split('|')[1] # the sortable component of the header

    # this loops through all the substrings and removes any substring that is a superString of a substring
    def removeAllSuperStrings(self,strings):
        all_strings = strings[:]
        index = 0
        for string in strings:
            for i in range(1,len(strings)):
                string2 = strings[i]
                if string in string2 and string != string2:
                    if string2 in all_strings:
                        all_strings.remove(string2)
                        index += 1
        return all_strings


    # This function writes the subsequnces to a file
    def writeSubstringsToFile(self,fil):
        '''
        here we loop through all the frames sorted by length in reverse order
        and print it in the appropriate format
        '''
        fil.write('{}\n'.format(self.header))
        fil.write('{}\n'.format(self.originalSequence))
        sorted_x = sorted(self.orderedSubSequences.items(), key=lambda x: x[1])
        for string in sorted_x:
            if string[0] in self.uniqueNonSuperSubstrings:
                fmtString = ''
                for i in range(0, string[1]):
                    fmtString += ('.')
                fil.write('{} {} \n'.format(fmtString,string[0]))

    # This function prints the subsequences to the console
    def printSubstrings(self):
        print (self.header)
        print(self.originalSequence)
        sorted_x = sorted(self.orderedSubSequences.items(), key=lambda x: x[1])
        for string in sorted_x:
            if string[0] in self.uniqueNonSuperSubstrings:
                fmtString = ''
                for i in range(0, string[1]):
                    fmtString += ('.')
                print('{}{}'.format(fmtString,string[0]))

    # This function generates all the substrings
    # And discards all subsequnces over 7
    def get_all_substrings(self,string):
      length = len(string)
      alist = []
      for i in range(length):
        for j in range(i,length):
          if len(string[i:j]) >7 : continue
          self.orderedSubSequences[string[i:j +1]] = i
          alist.append(string[i:j + 1]) 
      return alist


########################################################################
# Main
# Here is the main program
# 
#
########################################################################
   
from sequenceAnalysis import NucParams,FastAreader,ProteinParam
def main():

    # Clear the file if created previously, 
    # and open it 
    # for line in sys.stdin: 
    #     print (line)
    # print(sys.argv[0])
    orfReader = FastAreader('foo')


    # First we must extract all the fasta data
    allSequencesHead = []
    allSequences = []
    for head, seq in orfReader.readFastaStdIn(sys.stdin):
        nucParams = NucParams('')
        nucParams.addSequence(seq)
        nucParams.buildNucleotide()
        bases = (''.join(nucParams.codons))
        # convert these to sets to remove duplicate substrings
        seq = set(util.get_all_substrings(bases)) 
        allSequencesHead.append((bases,head))
        allSequences.append(seq)

    # Then we calculate the union of everything except the sequence we are on
    # And create a trnaFinder with the respective data

    trnaList = []
    baseIndex = 0
    for base,head in allSequencesHead:
        tmpSeq = allSequences.copy()
        del tmpSeq[baseIndex] # make sure to delete the current sequnce

        tmpUnique = set().union(*tmpSeq)
        trnaList.append(trnaFinder(base,head,tmpUnique))
        baseIndex += 1

    # write all the data to a file
    for trna in sorted(trnaList, key=lambda x: x.sortableHeader):
        trna.printSubstrings()


#######

import sys

if __name__ == "__main__":
    main() 


