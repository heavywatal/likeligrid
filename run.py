#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import wtl.options as wopt

import gzip
import json
import itertools

program = 'likeligrid'


def iter_args(infiles, range_s, concurrency, rest, epistasis, tp53):
    const = [program, '-j{}'.format(concurrency)] + rest
    axes = wopt.OrderedDict()
    axes['s'] = range_s
    for v in wopt.sequential(axes):
        args = wopt.make_args(v)
        for f in infiles:
            if epistasis:
                if '-g' in const:
                    npath = count_pathways_tsv(f)
                else:
                    npath = count_pathways_json(f)
                for x in itertools.combinations(range(npath), 2):
                    yield const + ['-e {} {}'.format(*x)] + args + [f]
            elif tp53:
                pair = tp53_pleiotropic_pair(f)
                yield const + ['-e {} {}'.format(*pair), '-p'] + args + [f]
            else:
                yield const + args + [f]


def read_pathways_tsv(infile):
    """Read result file"""
    with gzip.open(infile, 'rt') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            header = line.rstrip().split('\t')
            break
    header.pop(0)  # loglik
    return header


def count_pathways_tsv(infile):
    """Read result file"""
    header = read_pathways_tsv(infile)
    if header[-1] == 'pleiotropy':
        header.pop(-1)
    if ':' in header[-1]:
        header.pop(-1)
    return len(header)


def count_pathways_json(infile):
    """Read genotype file"""
    with gzip.open(infile, 'r') as fin:
        d = json.load(fin)
        return len(d['pathway'])


def tp53_pleiotropic_pair(infile):
    """Read result file"""
    columns = read_pathways_tsv(infile)
    return (columns.index('Cycle'), columns.index('Damage'))


def tp53_pleiotropic_pair_json(infile):
    """Read genotype file"""
    with gzip.open(infile, 'r') as fin:
        d = json.load(fin)
        p = d['pathway']
        return (p.index('Cycle'), p.index('Damage'))


def main():
    parser = wopt.ArgumentParser()
    parser.add_argument('-o', '--outdir', default='.stdout')
    parser.add_argument('--begin', type=int, default=4)
    parser.add_argument('--end', type=int, default=6)
    parser.add_argument('-e', '--epistasis', action='store_true')
    parser.add_argument('--tp53', action='store_true')
    parser.add_argument('infile', nargs='+')
    (args, rest) = parser.parse_known_args()

    print("cpu_count(): {}".format(wopt.cpu_count()))
    print('{} jobs * {} cores/job'.format(args.jobs, args.parallel))
    range_s = range(args.begin, args.end)
    it = iter_args(args.infile, range_s, args.parallel, rest,
                   args.epistasis, args.tp53)
    wopt.map_async(it, args.jobs, args.dry_run, outdir=args.outdir)
    print('End of ' + __file__)


if __name__ == '__main__':
    main()
