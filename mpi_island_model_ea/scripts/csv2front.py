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

cmaps = ['Blues_r','Oranges_r','Purples_r','Greens_r','Reds_r','Greys_r']

df = pd.read_csv(logfn)

#print(df)

fig, ax = plt.subplots()

ax.set_xlabel('Average Run Fitness')
ax.set_ylabel('Average Run Time')

ax.set_title('Front %d,%d' % (df['cycle'][0], df['rank'][0]))

front_groups = df.groupby(['rank','cycle'])

#df.loc[ [('at', 1),('at', 3),('at', 5)], 'Dwell']

for r, c in front_groups:
  #fcolor=plt.cm.get_cmap('tab20')
  #print('[%s,%s] (%f,%f) => ' % (f[0], f[1], fronts.nth(['o1_fitness'], fronts['o2_fitness'])))
  #fronts.plot.line(ax=ax, x='o1_fitness', y='o2_fitness', color=fcolor(f[0]))
            
#ax.legend(loc='upper left', frameon=False)
#plt.show()
