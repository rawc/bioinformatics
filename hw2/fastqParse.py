#!/usr/bin/env python3
# Name: Louis Chatta(lchatta)
# Group Members: none
"""
A fastQ sequence is taken as an input, with the 
assumed formatting of : split elements. With 
this assumption, we split the string on ':' and
then loop through all of the split items. 
Since the first item has a @ prepended to the string,
we check if the loop is at the first item, and if so we 
go ahead and print the element starting at index 1, in 
order to skip the '@'. In order to save lines of code, 
I saved all the print strings in an array that coordinates
with the index number. 
"""

fastQ = input("Enter a FASTQ sequence: ")
#split on ':'
fastQList = fastQ.split(':')
index = 0
#pre formatted print array
printList = ['Instrument = {}','Run ID = {}','Flow Cell ID = {}','Flow Cell Lane = {}','Tile Number = {}','X-coord = {}','Y-coord = {}']
for identifier in fastQList:
    if index == 0:
        print (printList[index].format(identifier[1:]))
    else:
        print (printList[index].format(identifier))
    index += 1


