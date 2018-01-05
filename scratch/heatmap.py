'''
Some stylized heatmap generators for displaying oligo interactions

-> borrowed heavily from http://stackoverflow.com/questions/14391959/heatmap-in-matplotlib-with-pcolor
'''
from primer3 import calcHeterodimerTm as chtm
import matplotlib.pyplot as plt
import numpy as np
import sys


MIN_INT = 0.001


def heatmap2D(data, row_labels, col_labels, tm_cutoff=None, binary=False,
              square_size=0.6):
    if tm_cutoff:
        data = np.copy(data)
        if binary:
            np.clip(data, a_min=tm_cutoff, a_max=tm_cutoff+MIN_INT, out=data)
        else:
            data = data - tm_cutoff
    fig, ax = plt.subplots(figsize=(data.shape[1] * square_size,
                                    data.shape[0] * square_size))
    heatmap = ax.pcolor(data, cmap=plt.cm.Reds, alpha=0.8)
    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
    ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)

    # want a more natural, table-like display
    ax.set_frame_on(False)
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(col_labels, minor=False, rotation=90)
    ax.set_yticklabels(row_labels, minor=False)
    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    return fig, ax, heatmap


def qcOligos(list1, list2, row_labels=None, col_labels=None, **kwargs):
    l1len = len(list1)
    l2len = len(list2)
    arr = np.zeros((l1len, l2len))
    col_labels = col_labels or ['2-%d' % i for i in range(l2len)]
    row_labels = row_labels or ['1-%d' % i for i in range(l1len)]
    for i in range(l1len):
        for j in range(l2len):
            arr[i][j] = chtm(list1[i], list2[j])
    return heatmap2D(arr, row_labels, col_labels, **kwargs)

