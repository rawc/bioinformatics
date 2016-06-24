#!/usr/bin/env python3
# Name: Louis Chatta (lchatta)
# Group Members: (lchatta, cjavery)
 
"""
This program calculates the physical-chemical properties of a protein sequence

input: a protein sequence 
    - as a one - letter encoding of the amino acid
    - a valid amino is a capital letter and is in the table
    - anything that violates this will not be stored in the protein string
output: number of amino acids and total molecular weight
       Molar extinction coefficient and Mass extinction coefficient
       theoretical isoelectric point(pI)
       amino acid composition
"""

class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O
 
    aa2mw = {
    'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
    'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
    'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
    'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }
    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34
    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}

    # the __init__ method requires a protein string to be provided, either as a
    # string or list of strings that will be concatenated

    """
    Here I loop through all the characters in the protein and and check if 
    they are allowed input, if so add to the splitAminos array. Then create
    the protein string
    """

    def __init__ (self, protein):
        splitAminos = []
        allowed_aminos = self.aa2mw.keys()
        for char in protein:
            if char in allowed_aminos:
                splitAminos.append(char)
    
        l = ''.join(splitAminos).split()
        self.protString = ''.join(l).upper()
        
        print(self.protString)


#Here the Count of the protein is returned
    def aaCount (self):
        return len(self.protString)

#Here the PI is calculated by finding the ph at which the net charge is closest to 0 and non-negative
    def pI (self):
        ph = 0.00
        lowestCharge = (99999,0.0)
        while(ph <= 14.0):
            currentCharge = self._charge_(ph)
            if currentCharge < lowestCharge[0] and self._charge_(ph) >=0: 
                lowestCharge = (currentCharge,ph)
            ph += 0.01

        return lowestCharge[1]

    """
    Here the dictionary of amino acids to their count is created and returned
    """
    def aaComposition (self) :
        compDict = {}
        for amino in self.aa2mw.keys():
            compDict[amino] = self.protString.count(amino)
        return compDict

    """
    Here the charge is calculated via the formula provided. Moreover, the acutual work is 
    done in the dotPh
    """
    def _charge_ (self, pH):
        firstCharge = self.dotPh(['R','K','H'],pH) + (10**self.aaNterm)/(10**self.aaNterm+10**pH)
        secondCharge = self.dotPh(['D','E','C','Y'],pH) + (10**pH)/(10**self.aaCterm+10**pH)
        return firstCharge - secondCharge
    """
    Here the inner sum of each of the aminos is calculated
    """
    def dotPh(self,allowed_aminos,pH):
        totalFirstCharge = 0
        for amino in self.protString:
            if amino in allowed_aminos:
                if amino in self.aa2chargePos.keys():
                    t = 10**self.aa2chargePos[amino]
                    b = 10**self.aa2chargePos[amino] + 10**pH
                    totalFirstCharge += t/b
                elif amino in self.aa2chargeNeg.keys():

                    t = 10**pH
                    b = 10**self.aa2chargeNeg[amino] + 10**pH
                    totalFirstCharge += t/b
                else:
                    t = 10**0
                    b = 10**0 + 10**pH
                    totalFirstCharge += t/b
        return totalFirstCharge

    """
    Here the associated count of the aminos is multiplied by the associated constant
    according to the formula provided.
    """
    def molarExtinction (self):
       return float(self.protString.count('Y')*self.aa2abs280['Y'] 
        + self.protString.count('W')*self.aa2abs280['W']
        + self.protString.count('C')*self.aa2abs280['C'])

    def massExtinction (self):
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    """
    Here the total molecularWeight is calculated by looping 
    through the sum of the molecularWeights associated to the amino acid
    and then subtracting the waters for the peptide bonds. 
    As can be seen in the fomula.
    """
    def molecularWeight (self):
        firstMolecWeight = 0.0
        length = sum(len(s) for s in self.protString)

        for amino in self.protString:
            firstMolecWeight += self.aa2mw[amino]
        return (firstMolecWeight - self.mwH2O * (length-1))
        

# Please do not modify any of the following.  This will produce a standard output that can be parsed

import sys
for inString in sys.stdin :
    myParamMaker = ProteinParam(inString)
    myAAnumber = myParamMaker.aaCount()
    print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
    print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
    print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
    print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
    print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
    print ("Amino acid composition:")
    myAAcomposition = myParamMaker.aaComposition()
    keys = list(myAAcomposition.keys())
    keys.sort()
    if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
    for key in keys :
        print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))