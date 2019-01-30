import sys
# We have a file with a bunch of stuff per feature, and we want the best per feature.
# Takes four arguments: the file, if there's a header, the column of features to unique, and the column of values to get
# the *lowest* of.
bests = {}

isHeader = bool(sys.argv[2])
keyIndex = int(sys.argv[3])
compareIndex = int(sys.argv[4])


for line in open(sys.argv[1], "r"):
	if isHeader:
		isHeader = False
		sys.stdout.write(line)
		continue
	fields = line.rstrip().split()
	key = fields[keyIndex]
	thisCompare = float(fields[compareIndex])
	if key in bests:
		oldCompare, oldLine = bests[key]
		if thisCompare < oldCompare:
			bests[key] = (thisCompare, line)
	else:
		bests[key] = (thisCompare, line)

for compare, line in bests.values():
	sys.stdout.write(line)

