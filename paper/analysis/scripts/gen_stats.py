#!/usr/bin/env python3
"""Derive prose statistics for the frozen latest-release comparison."""
import json
from _stats import Stats
from _scientific import PAPER
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
    execution = json.loads((PAPER / 'analysis/data/execution.json').read_text())
    for engine in ('upstream', 'plus'):
        attempts = [row for row in execution['attempts'] if row['engine'] == engine]
        for key, value in (('attempts', len(attempts)), ('failures', sum(row['status'] != 'complete' for row in attempts))):
            stats.add(f'execution.{engine}.{key}', value, fmt=',', desc=f'{engine} comparative {key}', between=(0,10000))
    named = json.loads((PAPER / 'analysis/data/named-modifications.json').read_text())
    stats.add('named.checks', named['checks'], fmt=',', desc='Named attachment CLI checks passed', between=(0,1000))
    return stats.write(inputs=INPUTS + REPORT_INPUTS + MATCHED_INPUTS + MASS_OFFSET_INPUTS + ['analysis/data/execution.json', 'analysis/data/named-modifications.json'])


if __name__ == '__main__':
    raise SystemExit(main())
