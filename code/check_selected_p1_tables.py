#!/usr/bin/env python3
"""Exact arithmetic audit of four selected, already derived P1 tables.

This standard-library-only program does not reproduce FASTQ processing,
alignment, observation rules or biological truth. It checks finite-record
arithmetic conditional on the supplied accepted E/NE/U assignments.

Example from a publication candidate:
  python3 code/check_selected_p1_tables.py --tables-dir tables \
      --output results/SELECTED_P1_ARITHMETIC_CHECK.json

With no arguments, use ../tables relative to this script when present, else
../manuscript_p1/tables in the development tree. Individual table paths may
override --tables-dir. No table, manuscript or other code file is modified.
"""
from __future__ import annotations

import argparse
import ast
import copy
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import sys


TABLE_NAMES = {
    'S1': 'Table_S1_Primary_Run_Manifest.csv',
    'S7': 'Table_S7_Coding_Witness_Counts.csv',
    'S8': 'Table_S8_Coding_Witness_Contrasts.csv',
    'S10': 'Table_S10_Resolved_Label_Revocations.csv',
}
NATIVE = 'native_local_peptide_on_witnessed_observed_coding'
PROGRAMMED = 'CTT_full_SE2_on_witnessed_observed_coding'
EVENTS = {NATIVE, PROGRAMMED}
OBSERVER = 'triple_observed_coding_witness'
# Exact fraction comparisons are used for the scientific bounds. Decimal pp
# columns are display approximations; this tolerance is in percentage points.
PP_DISPLAY_TOLERANCE = Fraction(1, 10**11)


class CheckFailure(ValueError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise CheckFailure(message)


def integer(value: str, location: str) -> int:
    require(re.fullmatch(r'\d+', value) is not None,
            f'{location}: expected a nonnegative integer, got {value!r}')
    return int(value)


def fraction(value: str, location: str) -> Fraction:
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise CheckFailure(f'{location}: invalid rational {value!r}') from error


def displayed_pp(row: dict, expected_l: Fraction, expected_h: Fraction,
                 location: str) -> None:
    for key, expected in [('L_pp', expected_l), ('H_pp', expected_h)]:
        actual = fraction(row[key], location + '.' + key)
        require(abs(actual - 100*expected) <= PP_DISPLAY_TOLERANCE,
                f'{location}.{key}: display does not match the exact bound')


def counts(row: dict, location: str) -> dict[str, int]:
    result = {key: integer(row[key], location + '.' + key)
              for key in ('N', 'E', 'NE', 'U')}
    require(result['N'] > 0, location + ': N must be positive')
    require(result['E'] + result['NE'] + result['U'] == result['N'],
            location + ': E + NE + U differs from N')
    return result


def maximum_reduction(k: int, sources: list[tuple[str, Fraction, int]]) -> tuple[Fraction, dict[str, int]]:
    """Best lower-bound reduction for k revocations on accepted capacities.

    Each treated E -> U reduces the lower contrast by 1/N_t. Each control
    NE -> U reduces it by 1/N_c. An exchange from a smaller rate to an available
    larger rate cannot decrease the reduction, so descending rates are optimal.
    """
    require(k >= 0, 'Negative revocation budget')
    require(k <= sum(capacity for _, _, capacity in sources),
            'Revocation budget exceeds accepted-label capacities')
    remaining = k
    reduction = Fraction(0)
    allocation = {name: 0 for name, _, _ in sources}
    for name, rate, capacity in sorted(sources, key=lambda item: item[1], reverse=True):
        take = min(capacity, remaining)
        allocation[name] = take
        reduction += take*rate
        remaining -= take
    require(remaining == 0, 'Incomplete capacity allocation')
    return reduction, allocation


def minimum_revocations(lower: Fraction, sources: list[tuple[str, Fraction, int]]) -> int:
    require(lower > 0, 'Selected revocation table requires a strictly positive initial lower bound')
    remaining = lower
    total = 0
    for _, rate, capacity in sorted(sources, key=lambda item: item[1], reverse=True):
        available = capacity*rate
        if remaining <= available:
            ratio = remaining/rate
            needed = (ratio.numerator + ratio.denominator - 1)//ratio.denominator
            require(0 <= needed <= capacity, 'Invalid critical capacity')
            return total + needed
        remaining -= available
        total += capacity
    raise CheckFailure('Accepted capacities cannot remove strict positivity')


def check_tables(tables: dict[str, list[dict]]) -> dict:
    manifest = {}
    hashes = set()
    for index, row in enumerate(tables['S1'], 1):
        loc = f'S1 row {index}'
        run = row['run_accession']
        require(run not in manifest, loc + ': duplicate run accession')
        require(row['condition'] in {'untreated', 'mock', 'PE treated'}, loc + ': unexpected condition')
        replica = integer(row['source_replica'], loc + '.source_replica')
        require(replica in {1, 2, 3}, loc + ': unexpected source replica')
        raw = integer(row['raw_reads'], loc + '.raw_reads')
        require(raw > 0, loc + ': empty run')
        sha = row['sequence_quality_stream_sha256']
        require(re.fullmatch(r'[0-9a-f]{64}', sha) is not None, loc + ': malformed stream hash')
        require(sha not in hashes, loc + ': repeated sequence/quality stream hash')
        hashes.add(sha)
        manifest[run] = row
    require(len(manifest) == 9, 'S1: expected exactly nine selected primary runs')
    require(sum(int(row['raw_reads']) for row in manifest.values()) == 177006,
            'S1: selected raw-record total differs from 177006')
    expected_roles = {(condition, str(replica))
                      for condition in ('untreated', 'mock', 'PE treated')
                      for replica in (1, 2, 3)}
    require({(row['condition'], row['source_replica']) for row in manifest.values()} == expected_roles,
            'S1: incomplete condition/source-replica coverage')

    by_run_event = {}
    for index, row in enumerate(tables['S7'], 1):
        loc = f'S7 row {index}'
        run, event = row['run_accession'], row['event']
        require(run in manifest, loc + ': run missing from S1')
        require(event in EVENTS, loc + ': unexpected event')
        require(row['observer'] == OBSERVER, loc + ': unexpected observer')
        require((run, event) not in by_run_event, loc + ': duplicate run/event')
        source = manifest[run]
        require(row['condition'] == source['condition'] and row['source_replica'] == source['source_replica'],
                loc + ': condition/replica differs from S1')
        require(row['universe_sha256'] == source['sequence_quality_stream_sha256'],
                loc + ': universe hash differs from S1')
        c = counts(row, loc)
        require(c['N'] == int(source['raw_reads']), loc + ': N differs from S1 raw_reads')
        lower, upper = Fraction(c['E'], c['N']), Fraction(c['E']+c['U'], c['N'])
        require(fraction(row['L_exact'], loc) == lower, loc + ': incorrect lower bound')
        require(fraction(row['H_exact'], loc) == upper, loc + ': incorrect upper bound')
        displayed_pp(row, lower, upper, loc)
        by_run_event[run, event] = (row, c)
    require(set(by_run_event) == {(run, event) for run in manifest for event in EVENTS},
            'S7: incomplete nine-run/two-event coverage')
    native = [c for (run, event), (_, c) in by_run_event.items() if event == NATIVE]
    require(sum(c['N'] for c in native) == 177006, 'S7: native N total differs from 177006')
    require(sum(c['U'] for c in native) == 6281, 'S7: native U total differs from 6281')

    contrasts = {}
    for index, row in enumerate(tables['S8'], 1):
        loc = f'S8 row {index}'
        event = row['event']
        require(event in EVENTS and row['observer'] == OBSERVER, loc + ': event/observer mismatch')
        trow, t = by_run_event[row['treated'], event]
        crow, c = by_run_event[row['comparator'], event]
        require(trow['condition'] == 'PE treated', loc + ': treated run is not PE treated')
        require(crow['condition'] == row['control'] and row['control'] in {'mock', 'untreated'},
                loc + ': control condition mismatch')
        require(trow['source_replica'] == crow['source_replica'] == row['preparation'],
                loc + ': source-replica pairing mismatch')
        lower = Fraction(t['E'], t['N']) - Fraction(c['E']+c['U'], c['N'])
        upper = Fraction(t['E']+t['U'], t['N']) - Fraction(c['E'], c['N'])
        require(fraction(row['L_exact'], loc) == lower, loc + ': incorrect contrast lower')
        require(fraction(row['H_exact'], loc) == upper, loc + ': incorrect contrast upper')
        displayed_pp(row, lower, upper, loc)
        for key, expected in [('treated_counts', t), ('control_counts', c)]:
            try:
                embedded = ast.literal_eval(row[key])
            except (ValueError, SyntaxError) as error:
                raise CheckFailure(loc + ': invalid embedded count dictionary') from error
            require(type(embedded) is dict and set(embedded) == {'N', 'E', 'NE', 'U'},
                    loc + ': malformed embedded count dictionary')
            require(all(type(value) is int for value in embedded.values()) and embedded == expected,
                    loc + ': embedded counts differ from S7')
        key = (row['preparation'], row['control'], event)
        require(key not in contrasts, loc + ': duplicate contrast')
        contrasts[key] = (lower, upper, t, c)
    expected_contrasts = {(str(replica), control, event)
                          for replica in (1, 2, 3) for control in ('mock', 'untreated') for event in EVENTS}
    require(set(contrasts) == expected_contrasts, 'S8: incomplete selected contrasts')

    revocations = {}
    certificates = []
    for index, row in enumerate(tables['S10'], 1):
        loc = f'S10 row {index}'
        key = (row['preparation'], row['control'], row['event'])
        require(key in contrasts and key not in revocations, loc + ': missing or duplicate contrast')
        lower, upper, t, c = contrasts[key]
        require(fraction(row['L_exact'], loc) == lower and fraction(row['H_exact'], loc) == upper,
                loc + ': initial bounds differ from S8')
        displayed_pp(row, lower, upper, loc)
        sources = [('control_NE', Fraction(1, c['N']), c['NE']),
                   ('treated_E', Fraction(1, t['N']), t['E'])]
        k = integer(row['minimum_revocations'], loc + '.minimum_revocations')
        expected_k = minimum_revocations(lower, sources)
        require(k == expected_k, loc + ': claimed revocations are not minimal')
        allocation = {name: integer(row[name], loc + '.' + name) for name, _, _ in sources}
        require(sum(allocation.values()) == k, loc + ': selected withdrawals do not sum to K')
        require(all(allocation[name] <= capacity for name, _, capacity in sources),
                loc + ': selected withdrawal exceeds source capacity')
        reduction = sum((allocation[name]*rate for name, rate, _ in sources), Fraction(0))
        best_k, optimal_allocation = maximum_reduction(k, sources)
        best_previous, _ = maximum_reduction(k-1, sources)
        require(reduction == best_k, loc + ': supplied allocation does not attain optimal reduction')
        after, previous = lower-best_k, lower-best_previous
        require(after <= 0 < previous, loc + ': K/K-1 minimality certificate failed')
        require(fraction(row['lower_after_K_exact'], loc) == after, loc + ': wrong lower after K')
        require(fraction(row['lower_after_K_minus_1_exact'], loc) == previous,
                loc + ': wrong optimal lower after K-1')
        revocations[key] = k
        certificates.append({'preparation':key[0], 'control':key[1], 'event':key[2],
            'minimum_revocations':k, 'capacities':{'control_NE':c['NE'], 'treated_E':t['E']},
            'rates_exact':{'control_NE':str(Fraction(1,c['N'])), 'treated_E':str(Fraction(1,t['N']))},
            'supplied_allocation':allocation, 'one_optimal_allocation':optimal_allocation,
            'lower_after_K_exact':str(after), 'best_lower_after_K_minus_1_exact':str(previous),
            'K_feasible_and_K_minus_1_insufficient':True})
    require(set(revocations) == expected_contrasts, 'S10: incomplete revocation certificates')
    return {'manifest_runs':len(manifest), 'count_rows':len(by_run_event),
        'contrast_rows':len(contrasts), 'minimum_revocation_certificates':len(certificates),
        'native_totals':{key:sum(c[key] for c in native) for key in ('N','E','NE','U')},
        'exact_count_balance_and_bounds':True, 'S8_rederived_from_S7':True,
        'S10_rederived_from_S7_capacities':True, 'revocation_certificates':certificates}


def negative_controls(tables: dict[str, list[dict]]) -> list[dict]:
    tests = []
    for name, table, column in [('altered_N','S7','N'),
                                ('altered_exact_bound','S7','H_exact'),
                                ('inflated_minimum_revocations','S10','minimum_revocations')]:
        changed = copy.deepcopy(tables)
        old = changed[table][0][column]
        changed[table][0][column] = str(Fraction(old)+1) if column=='H_exact' else str(int(old)+1)
        try:
            check_tables(changed)
        except (CheckFailure, KeyError) as error:
            tests.append({'control':name, 'in_memory_only':True, 'rejected':True, 'reason':str(error)})
        else:
            raise CheckFailure('Negative control was wrongly accepted: '+name)
    return tests


def main(argv: list[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    default_tables = root/'tables'
    if not all((default_tables/name).is_file() for name in TABLE_NAMES.values()):
        default_tables = root/'manuscript_p1/tables'
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--tables-dir', type=Path, default=default_tables)
    parser.add_argument('--output', type=Path, default=root/'results/SELECTED_P1_ARITHMETIC_CHECK.json')
    for key in TABLE_NAMES:
        parser.add_argument('--table-'+key.lower(), type=Path)
    args = parser.parse_args(argv)
    inputs = []
    try:
        tables = {}
        for key, name in TABLE_NAMES.items():
            path = getattr(args, 'table_'+key.lower()) or args.tables_dir/name
            data = path.read_bytes()
            with path.open(encoding='utf-8-sig', newline='') as stream:
                tables[key] = list(csv.DictReader(stream))
            inputs.append({'table':key, 'path':str(path.resolve()), 'bytes':len(data),
                           'sha256':hashlib.sha256(data).hexdigest(), 'rows':len(tables[key])})
        result = check_tables(tables)
        controls = negative_controls(tables)
        output = {'status':'PASS', 'scope':'Exact arithmetic of selected derived tables; conditional on supplied accepted labels.',
                  'inputs':inputs, 'verification':result, 'negative_controls':controls,
                  'arithmetic':'fractions.Fraction; no floating-point bounds or minimality decisions',
                  'pp_display_tolerance_exact':str(PP_DISPLAY_TOLERANCE),
                  'not_claimed':{'FASTQ_processing_reproduced':False, 'alignment_reproduced':False,
                                 'observer_reexecuted':False, 'source_labels_independently_authenticated':False,
                                 'new_biological_observations':False, 'biological_truth_certified':False}}
        exit_code = 0
    except (OSError, CheckFailure, KeyError, csv.Error, ValueError) as error:
        output = {'status':'FAIL', 'error':str(error), 'inputs_read':inputs,
                  'scope':'Selected derived-table arithmetic only; no FASTQ or biological validation.'}
        exit_code = 1
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')
    summary = {key:value for key,value in output.items() if key in ('status','error','scope')}
    if exit_code == 0:
        summary.update({'runs':result['manifest_runs'], 'count_rows':result['count_rows'],
                        'contrasts':result['contrast_rows'], 'revocation_certificates':result['minimum_revocation_certificates'],
                        'native_N':result['native_totals']['N'], 'native_U':result['native_totals']['U'],
                        'negative_controls_rejected':len(controls)})
    summary['output'] = str(args.output.resolve())
    print(json.dumps(summary, ensure_ascii=False))
    return exit_code


if __name__ == '__main__':
    sys.exit(main())
