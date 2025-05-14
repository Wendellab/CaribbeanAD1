#!/usr/bin/python3

import sys
import gzip

if len(sys.argv) < 3:
    print("Usage: plink2treemix.py [gzipped input file] [gzipped output file]")
    print("ERROR: improper command line")
    sys.exit(1)

infile = gzip.open(sys.argv[1], "rt")  # 'rt' for text mode
outfile = gzip.open(sys.argv[2], "wt")  # 'wt' for text mode

pop2rs = dict()
rss = list()
rss2 = set()

line = infile.readline()
line = infile.readline()
while line:
    line = line.strip().split()
    rs = line[1]
    pop = line[2]
    mc = line[6]
    total = line[7]
    if rs not in rss2:
        rss.append(rs)
    rss2.add(rs)
    if pop not in pop2rs:
        pop2rs[pop] = dict()
    if rs not in pop2rs[pop]:
        pop2rs[pop][rs] = " ".join([mc, total])
    line = infile.readline()

pops = list(pop2rs.keys())
outfile.write(" ".join(pops) + "\n")

for rs in rss:
    line_output = []
    for pop in pops:
        tmp = pop2rs[pop][rs].split()
        c1 = int(tmp[0])
        c2 = int(tmp[1])
        c3 = c2 - c1
        line_output.append(",".join([str(c1), str(c3)]))
    outfile.write(" ".join(line_output) + "\n")

infile.close()
outfile.close()
