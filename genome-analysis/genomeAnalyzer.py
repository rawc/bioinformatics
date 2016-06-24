#!/usr/bin/env python3
# Name: Louis Chatta (lchatta)
# Group Members: (lchatta, cjavery)
from sequenceAnalysis import NucParams,FastAreader,ProteinParam

"""
This program outputs the sequence length, GC content, and relative codon frequency for a given file containing genome information
""" 

HUNDRED = 100
NUC_DIVISOR = 1000000

class genomeCompare:
	def __init__ (self, filenames=['testGenome.fa'], compare=False):
		if compare:
			self.compareAnalysis(filenames)
		else:
			for fa in filenames:
				self.individualAnalysis(fa)

	def compareAnalysis(self,filenames):
		assert(len(filenames) == 2)
		print('{} {}'.format(filenames[0],filenames[1]))
		firstReader = FastAreader (filenames[0])
		secondReader = FastAreader(filenames[1])

		firstLength, firstNucParams = self.sequenceLength(firstReader)
		secondLength, secondNucParams = self.sequenceLength(secondReader)

		print('1st sequence length = {:.2f}Mb ----- 2nd sequence length = {:.2f}Mb'.format(firstLength,secondLength))  
		print ('sequence difference = {:.2f}Mb'.format(firstLength-secondLength))

		firstGC = self.gcContent(firstNucParams)
		secondGC = self.gcContent(secondNucParams)

		print('1st GC content = {:.1f} ----- 2nd GC content = {:.1f}%'.format(firstGC,secondGC))
		print ('GC difference = {:.1f}'.format(firstGC-secondGC))

		firstCodonCount = firstNucParams.codonComposition()
		firstAaComp = firstNucParams.aaComposition()
		firstNucComp = firstNucParams.nucComposition()

		secondCodonCount = secondNucParams.codonComposition()
		secondAaComp = secondNucParams.aaComposition()
		secondNucComp = secondNucParams.nucComposition()	

		for firstCodon,secondCodon in zip(sorted(firstCodonCount),sorted(secondCodonCount)):
			firstTotal = firstCodonCount[firstCodon]
			firstAmino = firstNucParams.rnaCodonTable[firstCodon]
			firstAminoCount = firstNucParams.aminoAcid.count(firstAmino)

			secondTotal = secondCodonCount[secondCodon]
			secondAmino = secondNucParams.rnaCodonTable[secondCodon]
			secondAminoCount = secondNucParams.aminoAcid.count(secondAmino)

			firstFinalTotal = firstTotal/firstAminoCount * HUNDRED
			secondFinalTotal = secondTotal/secondAminoCount * HUNDRED

			print('1st: {} : {} {:5.1f}% ({:6d}) ----- 2nd: {} : {} {:5.1f}% ({:6d}) ----- difference: {:5.1f}%'.format(firstCodon,firstAmino,firstFinalTotal,firstTotal,secondCodon,secondAmino,secondFinalTotal,secondTotal,firstFinalTotal-secondFinalTotal))	


	def individualAnalysis(self, fa):
		print ('Reading from file {}'.format(fa))
		myReader = FastAreader (fa)
		length,nucParams = self.sequenceLength(myReader)
		print('sequence length = {:.2f}Mb'.format(length))  

		#calculate and print GC content
		gc = self.gcContent(nucParams)
		print('GC content = {:.1f}%'.format(gc)) 

		#calculate and print relative codon usage
		codonCount = nucParams.codonComposition()
		aaComp = nucParams.aaComposition()
		nucComp = nucParams.nucComposition()
		for codon in sorted(codonCount):
			total = codonCount[codon]
			amino = nucParams.rnaCodonTable[codon]
			aminoCount = nucParams.aminoAcid.count(amino)

			finalTotal = total/aminoCount * HUNDRED
			print('{} : {} {:5.1f}% ({:6d})'.format(codon,amino,finalTotal,total))		


	def gcContent(self,nucParams):
		nucComp = nucParams.nucComposition()
		gc = nucComp['G'] + nucComp['C']
		gc = gc/nucParams.nucCount() * HUNDRED
		return gc

	def sequenceLength(self,myReader):

		nucParams = NucParams('')
		for head, seq in myReader.readFasta() :
			nucParams.addSequence(seq)
		nucParams.buildNucleotide()

		length = nucParams.nucCount()/NUC_DIVISOR
		return length, nucParams

genomeCompare()
