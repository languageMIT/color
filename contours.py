#!/usr/bin/python3
from __future__ import print_function
import sys

import pandas as pd
import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt

FILENAME = "output/special_focal_density.csv"
NUM_ROWS = 8
NUM_COLUMNS = 20
EPSILON = 10 ** -2

def plot_unique_hue_contours(d, **kwds):
    # Make data into a matrix
    nax = len(d['Language'].unique())
    fig, ax = plt.subplots(nrows=nax, sharex=True, sharey=True)
    for axis, lang in zip(ax, d['Language'].unique()):
        for term in d.loc[d['Language'] == lang, 'term'].unique():
            local_d = d[(d['Language'] == lang) & (d['term'] == term)]
            plot_single_term_contours(axis, local_d, term, **kwds)
            axis.set_title(lang)
    return fig, ax

def counts_to_grid_matrix(d):
    m = np.zeros((NUM_ROWS, NUM_COLUMNS))
    for i, row in enumerate(range(1, NUM_ROWS + 1)):
        for j, col in enumerate(range(1, NUM_COLUMNS + 1)):
            the_row = d[(d['Row'] == -row) & (d['Column'] == col)]
            if len(the_row) > 0:
                assert len(the_row) == 1
                m[i, j] = list(the_row['value'])[0]
    return np.flipud(m)

def get_levels(m, n):
    # levels such that the SUM of values above each level = x%
    total = m.sum()
    if hasattr(n, '__iter__'):
        breaks = np.array(n)
    else:
        breaks = np.linspace(0, 1, n + 2)[1:-1]

    breaks *= total
    breaks = list(breaks)
    
    def gen():
        for i in reversed(range(int(total))):
            subtracted = m - i
            total_above = np.sum(subtracted[subtracted > 0])
            yield i, total_above
    levels = []
    for i, total_above in gen():
        while breaks and total_above >= breaks[0]:
            levels.append(i)
            print(i, total_above, breaks[0], file=sys.stderr)
            breaks = breaks[1:]
            
    return list(reversed(levels))

def plot_single_term_contours(axis, d, term, n=[.05, .25, .5], smooth=10, edge=1/4):
    m = counts_to_grid_matrix(d)
    if term == 'green':
        term = '#5DFC0A'
    if smooth:
        m = scipy.ndimage.zoom(m * 5, smooth, mode='wrap')
    levels = get_levels(m, n)
    for i in range(len(levels)):
        if levels[i] == 0:
            levels[i] = EPSILON
    extent = (-edge, NUM_COLUMNS - 1 + edge, -edge, NUM_ROWS - 1 + edge)
    axis.contour(m, colors=term, levels=levels, linewidths=1, extent=extent)
    axis.yaxis.set_ticks(np.arange(NUM_ROWS))
    axis.yaxis.set_ticklabels(list(reversed("ABCDEFGH")))
    axis.yaxis.set_tick_params(direction='out')
    axis.xaxis.set_ticks(np.arange(NUM_COLUMNS))
    axis.xaxis.set_ticklabels(np.arange(NUM_COLUMNS) + 1)
    axis.xaxis.set_tick_params(direction='out')    

def main():
    d = pd.read_csv(FILENAME)
    plot_unique_hue_contours(d, smooth=10)
    plt.savefig("output/unique_hue_contours.pdf")

if __name__ == '__main__':
    main()
    
            
        
