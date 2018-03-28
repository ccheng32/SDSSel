#!/usr/bin/python
import sys
import os
import networkx as nx
import random

NUM = 100000000
TMPNAME = '.tmp'
if len(sys.argv) != 4 + 1:
  print 'usage: ./py [num viewers] [num channels] [density ratio] [channel density ratio]'
  sys.exit(1)

nv = int(sys.argv[1])
nc = int(sys.argv[2])
ra = float(sys.argv[3])
rc = float(sys.argv[4])

# output file
f = open(TMPNAME, 'w')

print 'generating graph g'
G = nx.fast_gnp_random_graph(nv, ra)
edgelist = nx.generate_edgelist(G, data=False)

edgenum = 0
for line in edgelist:
  a = int(line.split()[0])
  b = int(line.split()[1])
  f.write(str(a + 1) + ' ' + str(b + 1) + '\n')
  edgenum += 1

# generate channels
chanedgenum = 0
for i in range(nv + 1, nv + nc + 1):
  for j in range(i + 1, nv + nc + 1):
    f.write(str(i) + ' ' + str(j) + ' '  + str(random.uniform(0.0, 1.0)) +  '\n')
    chanedgenum += 1

scedgenum = 0
for i in range(nv + 1, nv + nc + 1):
  for j in range(1, nv + 1):
    if random.randint(0, NUM) < NUM * rc:
      f.write(str(i) + ' ' + str(j) + '\n')
      scedgenum += 1
       
    
f.close()

ff = open('v' + str(nv) + '.' + 'c' + str(nc) + '.' + 'd' + str(ra) + 'rc' + str(rc), 'w')

# print metadata
ff.write(str(nv) + ' ' + str(nc) + ' ' + str(edgenum)+ ' ' + str(chanedgenum)+ ' ' + str(scedgenum)+ '\n')
with open(TMPNAME, 'r') as infile:
  for line in infile:
    ff.write(line)
os.remove(TMPNAME)
ff.close()


