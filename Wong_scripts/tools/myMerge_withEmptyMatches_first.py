import sys, gzip

def myMerge(
        f1,
        f2,
        fout,
        keyFunction1=lambda x: x[0],
        keyFunction2=lambda y: y[0],
        printFunction=lambda x, y: "\t".join(x + y)+"\n",
        elsePrintFunction=None
        ):
    '''Some function of f1 and f2 are keys; print if matched
    
    keyFunctions takes as input the list of fields and output a key.
    Default is the first field.

    printFunction takes as input the list of values in f1 and f2 and
    returns a string to print on matches. Default is all f1 and f2.
    '''
    d={}
    outPairs = []
    everyFields1 = []
    if f1[-3:] == ".gz":
        file1 = gzip.open(f1, "r")
    else:
        file1 = open(f1, "r")
    for line in file1:
        fields1=line.rstrip("\n").split()
        key1=keyFunction1(fields1)
	everyFields1.append(fields1)
        if key1 in d:
            d[key1].append(fields1)
        else:
            d[key1] = [fields1]
    out=fout
    if f2[-3:] == ".gz":
        file2 = gzip.open(f2, "r")
    else:
        file2 = open(f2, "r")
    for line in file2:
        fields2=line.rstrip("\n").split()
        key2=keyFunction2(fields2)
        if key2 in d:
            for fields1 in d[key2]:
		outPairs.append([tuple(fields1), fields2])
    outD = {}
    for p1, p2 in outPairs:
    	if p1 in outD:
		outD[p1].append(p2)
	else:
		outD[p1] = [p2]
    #outD = dict(outPairs)
    for fields1 in map(tuple, everyFields1):
        if fields1 in outD:
		for fields2 in outD[fields1]:
			out.write(printFunction(list(fields1), fields2))
        else:
            out.write("\t".join(fields1) + "\n")


try:
	assert len(sys.argv) == 5
except AssertionError:
	print "need 4 arguments; offered %s: %s" % (len(sys.argv)-1, sys.argv[1:])

myMerge(sys.argv[1], sys.argv[2], sys.stdout, keyFunction1=lambda x: x[int(sys.argv[3])], keyFunction2=lambda y: y[int(sys.argv[4])], elsePrintFunction=lambda y: "\t".join(y) + "\n")
