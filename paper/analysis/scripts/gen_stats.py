#!/usr/bin/env python3
"""Derive prose statistics for the frozen latest-release comparison."""
from _stats import Stats
from _scientific import INPUTS, add_stats
from _report import REPORT_INPUTS, add_report_stats


def main():
    stats = Stats()
    add_stats(stats)
    add_report_stats(stats)
    return stats.write(inputs=INPUTS + REPORT_INPUTS)


if __name__ == '__main__':
    raise SystemExit(main())
