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

runs = df.groupby(['run'])

for rid, run in runs:
  rstats = run.agg({'cycle': ['min', 'max']})
  for gid, gen in run.groupby(['cycle']):
    gstats = gen.agg({'rank': ['min', 'max'], 'o1_fitness': ['min', 'max'], 'o2_fitness': ['min', 'max']})
    #ax.text(gstats['o1_fitness'].min(),gstats['o2_fitness'].max(), 'Cycle %d' % (gid))
    for fid, fronts in gen.groupby(['rank']):
      fstats = fronts.agg({'o1_fitness': ['min', 'max'], 'o2_fitness': ['min', 'max']})
      #print(fronts)
      #print('cycle %d of %d, front: %d of %d' % (gid, rstats.max(), fid, gstats['rank'].max()))
      clrmap = plt.cm.get_cmap(cmaps[fid % len(cmaps)])
      dist = np.linspace(0,fid,int(gstats['rank'].max()))
      norm = plt.Normalize(vmin=fid, vmax=fid)
      hue = clrmap(norm(fid))
      for f, front in fronts.groupby(['cycle','rank']):
        if(f[1] == 0 and f[0]%2 == 0):
          print(front)
          front.plot.line(x='o1_fitness', y='o2_fitness', ax=ax, style='o-', color=hue, legend=False)
          #front.plot.line(x='o1_fitness', y='o2_fitness', ax=ax, style='o-', color=hue, legend=False)
          #ax.text(fstats['o1_fitness'].min(), fstats['o2_fitness'].max(), '(%d)' % (gid))
          #ax.annotate('(%d)' % gid, xy=(fstats['o1_fitness'].min(), fstats['o2_fitness'].max()), xycoords='data', xytext=(-60, -30), textcoords='offset points', arrowprops=dict(arrowstyle="->", color='black'))
          #front.plot(ax=ax, x='o1_fitness', y='o2_fitness', label='Front %d, Cycle %d' % (fid, gid) if gid == 1 else "")
      
#ax.legend(loc='upper left', frameon=False)
plt.show()


#print(gstats)
    #gen.plot.scatter(ax=ax, x='o1_fitness', y='o2_fitness', c='rank', colorbar=False)
##for rid, run in runs:
#  print(rid)
  

#print(intervals)

#print('###########################')
  
#for cid, interval in intervals:

  #cres = intervals.agg({'cycle': ['mean', 'min', 'max']})
  #max_cycle = interval.max()['cycle']
  #max_cycle = intervals.max()['cycle']
  #print(cres)
  #print('== cycle %d of %d ================' % (cid[1], max_cycle))
  #for fid, front in interval.groupby(['rank']):
   ## print('-- rank %d of %d ---------------' % (fid, max_rank))
  

#print('====================')
#print(interval.min())
#max_cycle = interval.max()['cycle']
#max_rank = interval.max()['cycle']
#interval_x, interval_y = interval.get_data()
#lbl = 'Front ' + str(cid[1]) + ',' + str(cid[2])
#colormap = plt.cm.get_cmap(cmaps[cid[2]])
##colormap = colormap(np.linspace(0,1,int(max_cycle)))
#colormap = mpl.colors.ListedColormap(colormap[cid[2]:,:-1])
#color = colormap(np.linspace(0,1,20))
#interval.plot.scatter(ax=ax, x='o1_fitness', y='o2_fitness', c='cycle', cmap=colormap, colorbar=False)
#interval.plot.line(ax=ax, x='o1_fitness', y='o2_fitness', cmap=colormap, label=lbl, legend=False)
#interval.plot(ax=ax, x='o1_fitness', y='o2_fitness', label=lbl)
#plt.annotate(lbl, xy=(interval['o1_fitness'],interval['o2_fitness']), xytext=(2,2), textcoords="offset points", ha='center')
#ax = intervals.plot.scatter(x='o1_fitness', y='o2_fitness', c='rank')


#for cid, interval in intervals:

#  print(cid)
 
#  for fid, front in interval.groupby(['rank']):
  
#    print(fid)

#plt.show()




#ax = df.plot.scatter(x='o1_fitness', y='o2_fitness', label='Front')
#df.plot(kind='scatter', x='o1_fitness', y='o2_fitness', c='rank', cmap=cmaps[df['rank'][0]], ax=ax, colorbar=False)
#df.plot.scatter(x='o1_fitness', y='o2_fitness', colormap=cmaps[0], ax=ax)

