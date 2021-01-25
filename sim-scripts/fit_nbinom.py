#!/usr/bin/env python3

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
from scipy.special import comb
import matplotlib.pyplot as plt
from statsmodels.discrete.discrete_model import NegativeBinomial
import sys
import pandas as pd
import seaborn as sns

data_df = pd.read_csv(sys.argv[1], delimiter='\t')
data_df_wide = data_df.pivot(index='pos',columns='cell', values='gt')

data = data_df_wide.notna().sum(axis=1)

bin_heights, bin_borders, patches = plt.hist(data, bins=135, label='hist')
bin_centers = bin_borders[:-1] + np.diff(bin_borders) / 2

# nbinom = NegativeBinomial(bin_heights, bin_centers)
# results = nbinom.fit()
# print(results.summary())

def nbinom(x, n, p):
    '''
    x: data point
    n: positive parameter
    p: probability of success
    '''
    pdf = (p**n) * ((1-p)**x) * comb(n+x-1, n-1)
    return pdf

popt, popc = curve_fit(nbinom, bin_centers, bin_heights,p0=[2, 0.09])
print(popt)
fig, ax = plt.subplots()
sns.distplot(a=data, bins=135, hist=True, fit=stats.nbinom, label='hist', ax=ax)
x_s = np.linspace(min(data), max(data), num=1000)
y_s = nbinom(x_s, popt[0], popt[1])
ax.plot(x_s, y_s, label='curve_fit')
y_pre_fit = nbinom(x_s, 2, 0.09)
ax.plot(x_s, y_pre_fit, label='p0')
fig.suptitle('Curve fit for not NA in {}'.format(sys.argv[2]))
ax.legend()
fig.savefig("{}_notNA_hist.png".format(sys.argv[2]))
plt.close(fig)
