#!/usr/bin/env python

# Needs the python modules numpy and matplotlib

import sys, csv
import numpy as np; import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(sys.argv[1], sep='\t', header=False)
df.columns = ['Length', 'Count']

plt.figure()
df.plot(x='Length', y='Count', kind='bar')

plt.title('Length Distribution') 
plt.xlabel('Length (bp)')
plt.ylabel('Number of Sequences')
plt.savefig(sys.stdout, format='pdf')
