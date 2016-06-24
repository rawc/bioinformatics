#!/usr/bin/env python3
# Name: Louis Chatta(lchatta)
# Group Members: none
"""
This program takes a string input and loops through
all of the characters in the string as a capitalizied 
string. It then checks every chararacter to see if it is 
allowed. If it is not, then it starts counting the number
of bases to remove. If it is and the number of bases to remove
is greator than Zero, then insert the number in a {%d} format, and
insert the new character, and reset the counter number. If the 
number of bases to remove is still 0, then we are encountering
another base and we should simply insert the base. 
Notice that the each allowed base and counter string is added
to fixedDna array, in order to manipulate the new string. When 
we are done adding to the fixedDna, then we create a string 
from the array.
"""
dna = input("Enter a DNA sequence: ")
# the allowed bases saved in an array
allowedBases = ['A','C','G','T']
removedBaseCount = 0
# the data structure used to hold the new string
fixedDna = []
# loop through input as uppercase
for base in dna.upper():
	# increment counter if not in allowed bases
   if base not in allowedBases:
   		removedBaseCount += 1
   else:
   		# add the the number and {} if the next number
   		# is a base and we have already seen at least 
   		# one non-base
   		if removedBaseCount > 0:
   			fixedDna.append('{{{}}}'.format(removedBaseCount))
   			removedBaseCount = 0
   			fixedDna.append(base)
   		# if we have not seen at least one non base
   		# then just add the base
   		else:
   			fixedDna.append(base)
# join all the elements in the array and print
print (''.join(fixedDna))



