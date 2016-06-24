fastQ = input("Enter a FASTQ sequence: ")
fastQList = fastQ.split(':')
int index == 0
for identifier in fastQList:
    if index == 0:
        print 'Instrument = {}'.format(identifier[1:])
    if index == 1:
        print 'Run ID = {}'.format(identifier)
    if index == 2:
        print 'Flow Cell ID = {}'.format(identifier)
    if index == 3:
        print 'Flow Cell Lane = {}'.format(identifier)
    if index == 4:
        print 'Tile Number = {}'.format(identifier)
    if index == 5:
        print 'X-coord = {}'.format(identifier)                         
    if index == 6:
        print 'Y-coord = {}'.format(identifier)
