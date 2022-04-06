#!/usr/bin/env python

from sys import argv

import os.path
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

try:
    logfn = argv[1]
    if (not os.path.isfile(logfn)):
        raise ValueError("Please specify an input file.")
    pass
except Exception as e:
    print(logfn + "Not found, please specify a valid input file.")
    raise

cmaps = ['Blues_r', 'Oranges_r', 'Purples_r', 'Greens_r', 'Reds_r', 'Greys_r']

df = pd.read_csv(logfn, usecols=['cycle', 'front', 'o1_fitness', 'o2_fitness'])


fig, ax = plt.subplots()

fronts = df[df['cycle'] == 1].groupby(['front'])
f1 = df[df['front'] == 1]

for f, front in fronts:
    print(front)
    plt.plot('o2_fitness', 'o1_fitness', '-o', data=front, label=f)

plt.title("CYCLE")
plt.legend(loc="upper left")


#plt.scatter('o1_fitness', 'o2_fitness', c='front', s='front', data=df)
#plt.plot('o1_fitness', 'o2_fitness', data=df)
#plt.plot.line('o1_fitness', 'o2_fitness', data=df)

plt.show()

