#!/usr/bin/env python

from sys import argv

import os.path
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

print(argv[1])

try:
  logfn=argv[1]
  if(not os.path.isfile(logfn)):
    raise ValueError("Please specify an input file.")
  pass
except Exception as e:
  print(logfn + "Not found, please specify a valid input file.")
  raise

cmaps = ['Blues','Oranges','Purples','Greens','Reds','YlOrRd','Pastel1','Pastel2','Greys']

df = pd.read_csv(logfn)

print(df)

fig, ax = plt.subplots()

#ax = df.plot.scatter(x='o1_fitness', y='o2_fitness', label='Front')
#df.plot(kind='scatter', x='o1_fitness', y='o2_fitness', c='rank', cmap=cmaps[df['rank'][0]], ax=ax, colorbar=False)
#df.plot.scatter(x='o1_fitness', y='o2_fitness', colormap=cmaps[0], ax=ax)

ax.set_xlabel('Average Run Fitness')
ax.set_ylabel('Average Run Time')

ax.set_title('Front %d,%d' % (df['cycle'][0], df['rank'][0]))

intervals = df.groupby(['run','cycle','rank'])

ax = intervals.plot.scatter(x='o1_fitness', y='o2_fitness', c='rank')

for cid, interval in intervals:

  print(cid)
 
  for fid, front in interval.groupby(['rank']):
  
    print(fid)
  
  

#plt.show()
