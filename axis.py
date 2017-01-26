#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import os
import sys
import argparse
import numpy as np
import pandas as pd
pd.set_option('display.max_rows', 20)
pd.set_option('display.width', None)


def make_init_axes(infile, step):
    genotypes = pd.read_table(infile)
    v = np.linspace(1.0, 0.0, 1.0 / step, endpoint=False)
    df = pd.DataFrame()
    for name in genotypes.columns:
        df[name] = v
    return df


def read_best(infile):
    df = pd.read_table(infile, comment='#')
    return df.iloc[df['loglik'].idxmax(), 1:]


def make_vicinity(center, step, grid_density):
    lower = center - step
    upper = center + step
    df = pd.DataFrame()
    for (name, l, u) in zip(center.index, lower, upper):
        v = np.linspace(u, l, grid_density, endpoint=True)
        df[name] = v[0.0 <= v]
    return df


parser = argparse.ArgumentParser()
parser.add_argument('-g', '--grid-density', type=int, default=11)
parser.add_argument('-d', '--step', type=float, default=0.1)
parser.add_argument('-i', '--init', action='store_true')
parser.add_argument('infile')
args = parser.parse_args()

if args.init:
    df = make_init_axes(args.infile, args.step)
    df.to_csv(sys.stdout, '\t', index=False, float_format='%g')
else:
    best = read_best(args.infile)
    print(best)
    df = make_vicinity(best, args.step, args.grid_density)
    print(df)
    (base, ext) = os.path.splitext(os.path.basename(args.infile))
    outfile = base + '.next-axes.tsv'
    print('writing: ' + outfile)
    df.to_csv(outfile, '\t', index=False, float_format='%g')
