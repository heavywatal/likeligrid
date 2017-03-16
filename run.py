#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import wtl.options as wopt

program = 'likeligrid'


def iter(infiles, range_s):
    params = wopt.OrderedDict()
    params['s'] = range_s
    for p in wopt.sequential(params):
        for f in infiles:
            yield ' '.join([program, p, f])


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-j', '--jobs', type=int, default=wopt.cpu_count())
    parser.add_argument('--begin', type=int, default=2)
    parser.add_argument('--end', type=int, default=6)
    parser.add_argument('infile', nargs='+')
    (args, rest) = parser.parse_known_args()

    range_s = range(args.begin, args.end)
    wopt.shell.map(iter(args.infile, range_s), args.jobs, args.dry_run)
    print('End of ' + __file__)
