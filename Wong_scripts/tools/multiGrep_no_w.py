import sys, gzip

def multiGrep(f):
    '''print lines in f that contain lines from stdin'''
    L = sys.stdin.readlines()
    L = map(lambda x: x.rstrip(), L)
    if f[-3:] == ".gz":
        file = gzip.open(f, "r")
    else:
        file = open(f, "r")
    for line in file:
	for key in L:
            if key in line:
		sys.stdout.write(line)

multiGrep(sys.argv[1])
