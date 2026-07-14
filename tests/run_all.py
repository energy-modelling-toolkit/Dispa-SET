#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Dispa-SET test-suite runner.

Runs the full Dispa-SET test suite (unit + integration + failure +
ultimate), prints per-test pass/fail status to stdout, and writes a
Markdown report at ``tests/output/TEST_REPORT.md`` with a results
summary table.

Usage
-----

    # From the repository root (with the dispaset2 conda env active):
    python tests/run_all.py

    # Or with a tighter top-level timeout (seconds):
    python tests/run_all.py --timeout 300

    # With per-test wall-clock timing recorded in the report:
    python tests/run_all.py --profile

The runner uses pytest under the hood. It enables ``--tb=short`` and
captures stdout/stderr to keep the output digestible. The resulting
Markdown report can be committed as part of the development workflow
to track regressions.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
TESTS_DIR = REPO_ROOT / "tests"
OUTPUT_DIR = TESTS_DIR / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
REPORT_PATH = OUTPUT_DIR / "TEST_REPORT.md"
JSON_PATH = OUTPUT_DIR / "test_report.json"

# Test groups (folder relative to tests/ -> human label)
GROUPS = [
    ("unit", "Unit tests"),
    ("integration", "Integration tests"),
    ("failure", "Failure / error-message tests"),
    ("ultimate", "Ultimate end-to-end test"),
]


def _parse_junit_xml(xml_path: Path) -> list[dict]:
    """Parse a JUnit XML file and return a list of per-test timing dicts.

    Each dict has keys: ``id`` (str), ``time`` (float), ``status`` (str).
    """
    records: list[dict] = []
    if not xml_path.exists():
        return records
    try:
        tree = ET.parse(xml_path)
    except ET.ParseError:
        return records
    root = tree.getroot()
    # JUnit XML can have <testsuites><testsuite>…</testsuite></testsuites>
    # or a bare <testsuite> at the root.
    testsuites = root.findall(".//testcase")
    for tc in testsuites:
        classname = tc.get("classname", "")
        name = tc.get("name", "")
        test_id = f"{classname}::{name}" if classname else name
        elapsed = float(tc.get("time", 0.0))
        if tc.find("failure") is not None or tc.find("error") is not None:
            status = "FAILED"
        elif tc.find("skipped") is not None:
            status = "SKIPPED"
        else:
            status = "PASSED"
        records.append({"id": test_id, "time": elapsed, "status": status})
    return records

def run_pytest(target: Path, timeout: int, profile: bool = False) -> dict:
    """Run pytest on ``target`` and return a dict with the outcome.

    When *profile* is True a JUnit XML file is generated alongside the run
    and per-test wall-clock times are parsed from it.
    """
    junit_path = OUTPUT_DIR / f"_junit_{target.name}.xml"
    cmd = [
        sys.executable, "-m", "pytest", str(target),
        "--tb=short", "-q",
        "--maxfail=200",
        f"--timeout={timeout}",
    ]
    if profile:
        cmd += [f"--junit-xml={junit_path}"]
    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd, cwd=str(REPO_ROOT), capture_output=True, text=True,
            timeout=timeout + 30,
        )
        stdout = proc.stdout
        stderr = proc.stderr
        rc = proc.returncode
    except subprocess.TimeoutExpired as e:
        stdout = e.stdout or ""
        stderr = (e.stderr or "") + "\n[runner] timeout reached"
        rc = -1
    elapsed = time.time() - t0

    summary = parse_pytest_summary(stdout)
    summary.update(
        target=str(target.relative_to(REPO_ROOT)),
        returncode=rc,
        elapsed=elapsed,
        stdout=stdout,
        stderr=stderr,
    )
    if profile:
        summary["per_test_times"] = _parse_junit_xml(junit_path)
    return summary


def parse_pytest_summary(stdout: str) -> dict:
    """Extract per-test pass/fail/skip counts and short test list from
    a pytest stdout. Robust to slightly different pytest versions."""
    import re
    summary = {"passed": 0, "failed": 0, "skipped": 0, "errors": 0,
               "tests": []}
    pattern = re.compile(r"(\d+)\s+(passed|failed|skipped|error|errors|warning|warnings)")
    for line in stdout.splitlines():
        line_strip = line.strip()
        # Per-test PASSED/FAILED/SKIPPED in verbose pytest output.
        if "::" in line and any(tok in line for tok in
                                 (" PASSED", " FAILED", " SKIPPED",
                                  " ERROR")):
            parts = line.split()
            test_id = parts[0]
            status = "?"
            for tok in ("PASSED", "FAILED", "SKIPPED", "ERROR"):
                if tok in line:
                    status = tok
                    break
            summary["tests"].append({"id": test_id, "status": status})
        # Final summary line of pytest: e.g.
        #   "12 passed, 17 warnings in 0.81s"
        #   "===== 1 failed, 11 passed in 0.6s ====="
        if " in " in line_strip and re.search(r"\d+\s+(passed|failed|skipped|error)", line_strip):
            for n_str, kind in pattern.findall(line_strip):
                try:
                    n = int(n_str)
                except ValueError:
                    continue
                if kind == "passed":
                    summary["passed"] = n
                elif kind == "failed":
                    summary["failed"] = n
                elif kind == "skipped":
                    summary["skipped"] = n
                elif kind in ("error", "errors"):
                    summary["errors"] = n
    return summary


def render_md(group_results: list[tuple[str, str, dict]],
              total_elapsed: float,
              profile: bool = False) -> str:
    lines: list[str] = []
    lines.append("# Dispa-SET test-suite report\n")
    lines.append(f"_Total wall-clock runtime: **{total_elapsed:.1f} s**._\n")
    lines.append("")
    lines.append("## Summary by test group\n")
    lines.append("| Group | Path | Passed | Failed | Skipped | Errors "
                 "| Time (s) | Status |")
    lines.append("|---|---|---:|---:|---:|---:|---:|---|")
    overall_failed = 0
    for label, target, res in group_results:
        status_icon = "PASS" if res["failed"] == 0 and res["errors"] == 0 \
            else "FAIL"
        if res["returncode"] == 5:
            status_icon = "NO TESTS"
        overall_failed += res["failed"] + res["errors"]
        lines.append(
            f"| {label} | `{res['target']}` | {res['passed']} "
            f"| {res['failed']} | {res['skipped']} | {res['errors']} "
            f"| {res['elapsed']:.1f} | {status_icon} |"
        )
    lines.append("")

    # Per-test timing table (only when --profile is used)
    if profile:
        all_timings: list[dict] = []
        for label, target, res in group_results:
            for rec in res.get("per_test_times", []):
                all_timings.append({**rec, "group": label})
        if all_timings:
            all_timings.sort(key=lambda r: r["time"], reverse=True)
            lines.append("## Per-test wall-clock times (slowest first)\n")
            lines.append("| Group | Test | Time (s) | Status |")
            lines.append("|---|---|---:|---|")
            for rec in all_timings:
                lines.append(
                    f"| {rec['group']} | `{rec['id']}` "
                    f"| {rec['time']:.3f} | {rec['status']} |"
                )
            lines.append("")

    lines.append("## Details by test group\n")
    for label, target, res in group_results:
        lines.append(f"### {label} (`{res['target']}`)\n")
        lines.append("```")
        tail = "\n".join(res["stdout"].splitlines()[-50:])
        lines.append(tail or "<no stdout>")
        lines.append("```")
        if res["stderr"].strip():
            lines.append("\n<details><summary>stderr</summary>\n")
            lines.append("```")
            lines.append("\n".join(res["stderr"].splitlines()[-30:]))
            lines.append("```")
            lines.append("</details>\n")
        lines.append("")

    lines.append("---")
    if overall_failed == 0:
        lines.append("All tests passed.")
    else:
        lines.append(
            f"**{overall_failed} test failure(s) / error(s) detected.**"
        )

    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--timeout", type=int, default=300,
        help="Per-test pytest timeout in seconds (default 300)."
    )
    parser.add_argument(
        "--only", nargs="+", default=None,
        help="Run only the named groups (e.g. unit integration)."
    )
    parser.add_argument(
        "--profile", action="store_true", default=False,
        help=(
            "Record the wall-clock time of each individual test and include "
            "a per-test timing table in the Markdown report. Useful for "
            "detecting regressions in test cost."
        ),
    )
    args = parser.parse_args()

    t0 = time.time()
    group_results = []
    for folder, label in GROUPS:
        if args.only and folder not in args.only:
            continue
        target = TESTS_DIR / folder
        if not target.exists():
            continue
        print(f"\n=== Running {label} ({target.relative_to(REPO_ROOT)}) ===")
        res = run_pytest(target, timeout=args.timeout, profile=args.profile)
        print(
            f"   passed={res['passed']} failed={res['failed']} "
            f"skipped={res['skipped']} errors={res['errors']} "
            f"({res['elapsed']:.1f}s)"
        )
        group_results.append((label, str(target), res))
    total_elapsed = time.time() - t0

    md = render_md(group_results, total_elapsed, profile=args.profile)
    REPORT_PATH.write_text(md, encoding="utf-8")
    JSON_PATH.write_text(
        json.dumps(
            {"groups": [
                {"label": l, "target": t,
                 "summary": {k: v for k, v in r.items()
                             if k in ("passed", "failed", "skipped",
                                      "errors", "elapsed", "returncode",
                                      "target")}}
                for (l, t, r) in group_results
            ], "total_elapsed": total_elapsed},
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"\nWrote report to {REPORT_PATH}")
    print(f"Wrote JSON summary to {JSON_PATH}")
    overall_failed = sum(r["failed"] + r["errors"]
                         for _, _, r in group_results)
    return 0 if overall_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
