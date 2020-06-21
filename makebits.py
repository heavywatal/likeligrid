#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import sys
import json
import numpy as np


def sample_values(n=50000):
    v = np.random.randint(3, 32, n)
    s1 = (4, 8, 16)
    s = np.fromiter(('{:05b}'.format(x) for x in v if x not in s1), 'U8')
    return s.tolist()


def make_dict():
    d = {"pathway": ["A", "B"], "annotation": ["00111", "11100"]}
    d['sample'] = sample_values()
    return d


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args()
    json.dump(make_dict(), args.outfile)


if __name__ == '__main__':
    main()
