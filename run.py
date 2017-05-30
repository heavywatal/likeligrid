#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import wtl.options as wopt

program = 'likeligrid'


def iter_args(infiles, range_s, concurrency, rest):
    const = [program, '-j{}'.format(concurrency)] + rest
    axes = wopt.OrderedDict()
    axes['s'] = range_s
    for v in wopt.sequential(axes):
        args = wopt.make_args(v)
        for f in infiles:
            yield const + args + [f]


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-j', '--jobs', type=int, default=wopt.cpu_count())
    parser.add_argument('--begin', type=int, default=2)
    parser.add_argument('--end', type=int, default=6)
    parser.add_argument('infile', nargs='+')
    (args, rest) = parser.parse_known_args()

    range_s = range(args.begin, args.end)
    wopt.map_async(iter_args(args.infile, range_s, args.jobs, rest),
                   1, args.dry_run)
    print('End of ' + __file__)


if __name__ == '__main__':
    main()
