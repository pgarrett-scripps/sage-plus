#!/usr/bin/env python3
"""Derive prose statistics for the frozen latest-release comparison."""
from _stats import Stats
from _scientific import INPUTS, add_stats
from _report import REPORT_INPUTS, add_report_stats
from _matched_fdp import MATCHED_INPUTS, add_matched_stats
from _mass_offset import MASS_OFFSET_INPUTS, add_mass_offset_stats


def main():
    stats = Stats()
    add_stats(stats)
    add_report_stats(stats)
    add_matched_stats(stats)
    add_mass_offset_stats(stats)
    return stats.write(inputs=INPUTS + REPORT_INPUTS + MATCHED_INPUTS + MASS_OFFSET_INPUTS)


if __name__ == '__main__':
    raise SystemExit(main())
