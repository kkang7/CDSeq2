"""Independently verify Tables S20-S22 against per-test input tables.

Usage: python3 verify_eqtl_accuracy.py OUTPUT_FOLDER
Only the Python standard library is required. No inference is rerun.
"""
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
tables = root / "supplementary_tables"


def read(name):
    with (tables / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def number(x):
    return float("nan") if x in ("NA", "", None) else float(x)


tests = [r for r in read("S8_eqtl_tests.csv")
         if r["engine"] in ("CDSeq2", "BayesPrism")]
for r in tests:
    for k in ("q", "effect", "injected_sign", "signal_ratio"):
        r[k] = number(r[k])
    r["replicate"] = int(r["replicate"])
    r["injected"] = r["injected"] == "TRUE"
    r["selected"] = False

selection_groups = defaultdict(list)
for r in tests:
    if r["engine"] == "CDSeq2":
        selection_groups[(r["scenario"], r["cell_type"])].append(r)
features = ["scenario_GAM6_score", "scenario_contribution_score",
            "scenario_expected_count_score"]
for rows in selection_groups.values():
    eligible = [r for r in rows if all(math.isfinite(number(r[k])) for k in features)]
    eligible.sort(key=lambda r: (-number(r[features[0]]), r["gene"]))
    for r in eligible[:int(0.4 * len(eligible) + 1e-9)]:
        r["selected"] = True

exc_selection = {(r["scenario"], r["gene"]): r["selected"] == "TRUE"
                 for r in read("S18_eqtl_before_after.csv") if r["engine"] == "CDSeq2"}
for r in tests:
    if r["engine"] == "CDSeq2" and r["cell_type"] == "Excit":
        assert r["selected"] == exc_selection[(r["scenario"], r["gene"])]

# Confirm a common truth universe before method-specific availability/selection.
truth = defaultdict(dict)
for r in tests:
    truth[r["engine"]][(r["scenario"], r["gene"], r["cell_type"])] = (
        r["injected"], r["injected_sign"])
assert truth["CDSeq2"] == truth["BayesPrism"]


def key(r):
    return (r["workflow"], float(r["signal_ratio"]), int(r["replicate"]), r["cell_type"])


def close(actual, expected, context):
    actual = number(actual)
    assert ((math.isnan(actual) and math.isnan(expected)) or
            math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-12)), (context, actual, expected)


def compute(rows, threshold, positive_total):
    calls = [r for r in rows if math.isfinite(r["q"]) and r["q"] <= threshold]
    tp = sum(r["injected"] and math.isfinite(r["effect"]) and
             (1 if r["effect"] > 0 else -1 if r["effect"] < 0 else 0) == r["injected_sign"]
             for r in calls)
    noninjected = sum(not r["injected"] for r in calls)
    fp = len(calls) - tp
    return dict(calls=len(calls), TP=tp, FP=fp, FP_noninjected=noninjected,
                FP_wrong_sign=fp-noninjected, FN=positive_total-tp,
                n_injected_all=positive_total,
                precision=tp/len(calls) if calls else math.nan,
                recall=tp/positive_total if positive_total else math.nan,
                F1=2*tp/(len(calls)+positive_total) if positive_total else math.nan,
                directional_FDP=fp/len(calls) if calls else 0,
                noninjected_FDP=noninjected/len(calls) if calls else 0)


def verify_metrics(row, expected):
    for field, value in expected.items():
        close(row[field], value, (key(row), field, row.get("q_threshold")))


curves = defaultdict(list)
for r in read("S20_precision_recall_curves.csv"):
    curves[key(r)].append(r)
fixed = {key(r): r for r in read("S21_precision_recall_F1.csv")}
matched = defaultdict(list)
for r in read("S22_matched_empirical_FDP.csv"):
    matched[key(r)].append(r)
assert set(curves) == set(fixed) == set(matched)
curve_count = matched_count = 0
for group, curve_rows in curves.items():
    workflow, strength, repetition, cell_type = group
    engine = "BayesPrism" if workflow == "BayesPrism standard" else "CDSeq2"
    universe = [r for r in tests if r["engine"] == engine and r["signal_ratio"] == strength
                and (repetition == 0 or r["replicate"] == repetition)
                and (cell_type == "ALL_TYPES" or r["cell_type"] == cell_type)]
    positive_total = sum(r["injected"] for r in universe)
    retained = [r for r in universe if r["selected"]] if workflow == "CDSeq2 GAM6" else universe
    availability = dict(n_tests_universe=len(universe), n_retained=len(retained),
                        n_available_before_selection=sum(math.isfinite(r["q"]) for r in universe),
                        n_available_retained=sum(math.isfinite(r["q"]) for r in retained))
    observed_q = sorted({r["q"] for r in retained if math.isfinite(r["q"])})
    exported_q = [number(r["q_threshold"]) for r in curve_rows if r["calls"] != "0"]
    assert len(observed_q) == len(exported_q)
    for a, b in zip(observed_q, exported_q):
        close(a, b, (group, "exhaustive thresholds"))
    assert len(curve_rows) == len(observed_q) + 1
    # Exported decimals may round by ~1e-15. Recover the original source cutoff
    # by its position, then count calls directly instead of using exported counts.
    attainable = []
    for row, cutoff in zip(curve_rows, [-math.inf] + observed_q):
        expected = compute(retained, cutoff, positive_total)
        verify_metrics(row, expected)
        for field, value in availability.items():
            close(row[field], value, (group, field))
        attainable.append((cutoff, expected))
        curve_count += 1
    verify_metrics(fixed[group], compute(retained, .05, positive_total))
    for row in matched[group]:
        cap = float(row["FDP_cap"])
        feasible = [(cutoff, metrics) for cutoff, metrics in attainable
                    if metrics["directional_FDP"] <= cap + 1e-12]
        chosen_q, expected = min(feasible, key=lambda pair: (-pair[1]["TP"], pair[1]["FP"], pair[0]))
        verify_metrics(row, expected)
        close(row["q_threshold"], chosen_q if math.isfinite(chosen_q) else math.nan,
              (group, "optimal threshold", cap))
        matched_count += 1

print(f"PASS: {len(curves)} workflow/strength/repetition/type groups")
print(f"PASS: {curve_count} exhaustive precision-recall rows, {len(fixed)} fixed-threshold rows")
print(f"PASS: {matched_count} matched-cap rows, maximum TP and all tie rules")
print("PASS: common injected truth, original GAM6 membership, full-target recall denominators")
print("PASS: all BayesPrism available tests retained without output selection")
