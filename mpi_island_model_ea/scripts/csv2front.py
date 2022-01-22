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
ax.set_title('Meta EA Pareto Fronts at Generation ' + str(df['cycle'][0]))

intervals = df.groupby(['run','cycle','rank'])

imin = df['cycle'].min()
imax = df['cycle'].max()
ilen = df['cycle'].count()
isum = df['cycle'].sum()
  
for cid, interval in intervals:
  
  print('*** dataframe ***')
  print('id=%d,%d,%d | df_min=%d | df_max=%d | df_count=%d | df_sum=%d' % (cid[0], cid[1], cid[2], imin, imax, ilen, isum))
    
  cmin = interval['rank'].min()
  cmax = interval['rank'].max()
  clen = interval['rank'].count()
  csum = interval['rank'].sum()
  
  cidx = 1 - (cid[1] / (clen * cid[1]))
  cdist = np.linspace(0, clen, clen)
  cnorm = mpl.colors.Normalize(vmin=np.min(cdist), vmax=np.max(cdist))
  cmap = plt.cm.get_cmap(cmaps[cid[2]])
  ccmap = cmap(cnorm(cdist))
      
  for fid, front in interval.groupby(['rank']):
  
    print('*** interval ***')
    print('id=%d,%d,%d | interval_min=%d | max=%d | count=%d | sum=%d' % (cid[0], cid[1], cid[2], cmin, cmax, clen, csum))
  
    #test1 = interval['cycle'].agg([np.sum, np.mean, np.std])
    #print(test1)
    
    #test2 = interval['cycle'].count()
    #print(test2)
  
    #test3 = front.count()
    #print(test3)
    
    #test4 = intervals['cycle'].count()
    #print(test4)
  
    fmin = interval['rank'].min()
    fmax = interval['rank'].max()
    flen = interval['rank'].count()
    fdist = np.linspace(0, 1, flen)
    #fnorm = plt.Normalize(vmin=np.min(fdist), vmax=np.max(fdist))
    fmap = plt.cm.get_cmap(cmaps[fid])
    #fcmap = fmap(fnorm(fdist))
    #fccmap = cmap(cnorm(fdist))
    #fidx = fid / flen
    colormap = plt.cm.get_cmap(cmaps[cid[2]])
    colormap = colormap(np.linspace(0,1,int(8)))
    colormap = mpl.colors.ListedColormap(colormap[cid[2]:,:-1])

    print('*** front ***')
    #print('fid=%d | ranks=%d | min=%d | max=%d | index=%f | distance=' % (fid, flen, fmin, fmax, fidx))
  
    front.plot(ax=ax, kind='scatter', x='o1_fitness', y='o2_fitness', c='rank', cmap=colormap, zorder=1, alpha=1.0, colorbar=False)
    front.plot.line(x='o1_fitness', y='o2_fitness', ax=ax, style='-o', cmap=colormap, alpha=1.0, legend=False)
    
    #fstep = front['cycle'].max() / flen
    #fnorm = (front - front.mean()) / (front.max() - front.min())

#df.plot.scatter(x='o1_fitness', y='o2_fitness', c='rank', marker='o', linestyle='-', ax=ax, colormap=fmap, alpha=fidx)
#df.plot(kind='scatter', x='o1_fitness', y='o2_fitness', c='rank')
#df.plot.scatter(x='o1_fitness', y='o2_fitness', c=rank, ax=ax)

#plt.legend()
plt.show()
