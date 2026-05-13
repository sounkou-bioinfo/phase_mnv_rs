#!/usr/bin/env python3
"""Pairwise phase-concordance reproducibility matrix over VCF/BCF manifests."""
from __future__ import annotations

import argparse
import csv
import itertools
import math
import os
import re
import shutil
import statistics
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run phase_compare for all within-method pairs in a manifest and "
            "summarize replicate reproducibility."
        )
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="TSV with columns method, run, vcf and optional sample",
    )
    parser.add_argument("--out-dir", required=True, help="Output directory")
    parser.add_argument(
        "--phase-compare",
        default="target/release/phase_compare",
        help="phase_compare binary (default: target/release/phase_compare)",
    )
    parser.add_argument("--sample", help="Sample name used for every comparison")
    parser.add_argument(
        "--threads", type=int, default=1, help="htslib reader threads (default: 1)"
    )
    parser.add_argument(
        "--only-snvs",
        action="store_true",
        help="Pass --only-snvs to phase_compare",
    )
    parser.add_argument(
        "--write-pairs",
        action="store_true",
        help="Also write phase_compare --pair-tsv for each pair",
    )
    parser.add_argument(
        "--limit-pairs",
        type=int,
        help="Optional maximum number of pairs per method, useful for smoke runs",
    )
    return parser.parse_args()


def safe_name(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return safe or "unnamed"


def read_manifest(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise SystemExit(f"manifest not found: {path}")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"method", "run", "vcf"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                f"manifest missing required column(s): {', '.join(sorted(missing))}"
            )
        rows = []
        for row_no, row in enumerate(reader, start=2):
            if None in row:
                raise SystemExit(f"manifest row {row_no} has more fields than the header")
            row = {key: (value or "").strip() for key, value in row.items()}
            for key in required:
                if not row[key]:
                    raise SystemExit(f"manifest row {row_no} has empty {key}")
            vcf = Path(row["vcf"])
            if not vcf.is_absolute():
                vcf = (path.parent / vcf).resolve()
            if not vcf.exists():
                raise SystemExit(f"manifest row {row_no} VCF/BCF not found: {row['vcf']}")
            row["vcf"] = str(vcf)
            rows.append(row)
    return rows


def command_for_pair(
    args: argparse.Namespace,
    row1: dict[str, str],
    row2: dict[str, str],
    pair_tsv: Path | None,
) -> list[str]:
    cmd = [args.phase_compare, "--threads", str(args.threads)]
    if args.only_snvs:
        cmd.append("--only-snvs")
    if pair_tsv is not None:
        cmd.extend(["--pair-tsv", str(pair_tsv)])
    if args.sample:
        cmd.extend(["--sample", args.sample])
    else:
        sample = row1.get("sample", "") or row2.get("sample", "")
        if sample:
            cmd.extend(["--sample", sample])
    cmd.extend([row1["vcf"], row2["vcf"]])
    return cmd


def total_row(summary_text: str) -> tuple[list[str], list[str]]:
    lines = [line for line in summary_text.splitlines() if line.strip()]
    if not lines:
        raise RuntimeError("phase_compare produced empty summary")
    header = lines[0].split("\t")
    for line in lines[1:]:
        fields = line.split("\t")
        if fields and fields[0] == "TOTAL":
            if len(fields) != len(header):
                raise RuntimeError(
                    "phase_compare TOTAL row has "
                    f"{len(fields)} columns but header has {len(header)}"
                )
            return header, fields
    raise RuntimeError("phase_compare summary lacks TOTAL row")


def parse_number(value: str) -> float | None:
    if value in {"", "NA", "nan", "NaN"}:
        return None
    try:
        out = float(value)
    except ValueError:
        return None
    if math.isnan(out):
        return None
    return out


def enforce_sample_scope(
    args: argparse.Namespace, by_method: dict[str, list[dict[str, str]]]
) -> None:
    if args.sample:
        return
    for method, rows in by_method.items():
        samples = sorted({row.get("sample", "") for row in rows if row.get("sample", "")})
        if len(samples) > 1:
            raise SystemExit(
                "manifest has multiple samples for method "
                f"'{method}'; split methods by sample or pass --sample explicitly"
            )


def write_method_summary(path: Path, long_rows: list[dict[str, str]], metrics: Iterable[str]) -> None:
    by_method: dict[str, list[dict[str, str]]] = {}
    for row in long_rows:
        by_method.setdefault(row["method"], []).append(row)

    with path.open("w", newline="") as handle:
        fieldnames = ["method", "pairs", "metric", "n", "min", "median", "mean", "max", "sum"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for method in sorted(by_method):
            rows = by_method[method]
            for metric in metrics:
                values = [
                    value
                    for value in (parse_number(row.get(metric, "")) for row in rows)
                    if value is not None
                ]
                if not values:
                    continue
                writer.writerow(
                    {
                        "method": method,
                        "pairs": len(rows),
                        "metric": metric,
                        "n": len(values),
                        "min": f"{min(values):.6g}",
                        "median": f"{statistics.median(values):.6g}",
                        "mean": f"{statistics.mean(values):.6g}",
                        "max": f"{max(values):.6g}",
                        "sum": f"{sum(values):.6g}",
                    }
                )


def remove_path(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    elif path.exists() or path.is_symlink():
        path.unlink()


def replace_output_dir(work_dir: Path, out_dir: Path) -> None:
    backup = out_dir.with_name(f".{out_dir.name}.backup.{os.getpid()}")
    if backup.exists() or backup.is_symlink():
        remove_path(backup)
    had_previous = out_dir.exists() or out_dir.is_symlink()
    if had_previous:
        out_dir.rename(backup)
    try:
        work_dir.rename(out_dir)
    except OSError:
        if had_previous and not (out_dir.exists() or out_dir.is_symlink()):
            backup.rename(out_dir)
        raise
    if had_previous:
        remove_path(backup)


def main() -> None:
    args = parse_args()
    if args.threads < 1:
        raise SystemExit("--threads must be >= 1")
    if args.limit_pairs is not None and args.limit_pairs < 0:
        raise SystemExit("--limit-pairs must be >= 0")
    manifest = Path(args.manifest).resolve()
    out_dir = Path(args.out_dir).resolve()

    rows = read_manifest(manifest)
    by_method: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_method.setdefault(row["method"], []).append(row)
    enforce_sample_scope(args, by_method)

    pair_plan: list[tuple[str, dict[str, str], dict[str, str]]] = []
    for method in sorted(by_method):
        method_rows = sorted(by_method[method], key=lambda row: row["run"])
        pairs = itertools.combinations(method_rows, 2)
        if args.limit_pairs is not None:
            pairs = itertools.islice(pairs, args.limit_pairs)
        pair_plan.extend((method, row1, row2) for row1, row2 in pairs)
    if not pair_plan:
        raise SystemExit("manifest has no within-method pairs to compare")

    out_dir.parent.mkdir(parents=True, exist_ok=True)
    work_dir = Path(
        tempfile.mkdtemp(prefix=f".{out_dir.name}.tmp.", dir=out_dir.parent)
    )
    try:
        pairwise_dir = work_dir / "pairwise"
        pair_tsv_dir = work_dir / "pairs"
        pairwise_dir.mkdir(parents=True, exist_ok=True)
        if args.write_pairs:
            pair_tsv_dir.mkdir(parents=True, exist_ok=True)

        long_rows: list[dict[str, str]] = []
        summary_header: list[str] | None = None

        for pair_no, (method, row1, row2) in enumerate(pair_plan, start=1):
            stem = (
                f"{pair_no:06d}.{safe_name(method)}."
                f"{safe_name(row1['run'])}_vs_{safe_name(row2['run'])}"
            )
            summary_path = pairwise_dir / f"{stem}.summary.tsv"
            pair_tsv = pair_tsv_dir / f"{stem}.pairs.tsv" if args.write_pairs else None
            cmd = command_for_pair(args, row1, row2, pair_tsv)
            try:
                result = subprocess.run(cmd, text=True, capture_output=True, check=True)
            except FileNotFoundError:
                raise SystemExit(f"phase_compare not found: {args.phase_compare}") from None
            except OSError as exc:
                reason = exc.strerror or str(exc)
                raise SystemExit(
                    f"cannot execute phase_compare '{args.phase_compare}': {reason}"
                ) from None
            except subprocess.CalledProcessError as exc:
                (pairwise_dir / f"{stem}.stdout.txt").write_text(exc.stdout or "")
                (pairwise_dir / f"{stem}.stderr.txt").write_text(exc.stderr or "")
                raise SystemExit(
                    "phase_compare failed for "
                    f"method='{method}' run1='{row1['run']}' run2='{row2['run']}' "
                    f"with exit code {exc.returncode}"
                ) from None
            summary_path.write_text(result.stdout)
            if result.stderr:
                (pairwise_dir / f"{stem}.stderr.txt").write_text(result.stderr)
            try:
                header, fields = total_row(result.stdout)
            except RuntimeError as err:
                raise SystemExit(
                    "invalid phase_compare summary for "
                    f"method='{method}' run1='{row1['run']}' run2='{row2['run']}': "
                    f"{err}"
                ) from None
            if summary_header is None:
                summary_header = header
            elif header != summary_header:
                raise SystemExit(
                    "phase_compare summary header changed for "
                    f"method='{method}' run1='{row1['run']}' run2='{row2['run']}'"
                )
            long_row = {"method": method, "run1": row1["run"], "run2": row2["run"]}
            long_row.update(dict(zip(header, fields)))
            long_rows.append(long_row)

        if summary_header is None:
            raise SystemExit("manifest has no within-method pairs to compare")

        long_path = work_dir / "pairwise_long.tsv"
        with long_path.open("w", newline="") as handle:
            fieldnames = ["method", "run1", "run2"] + summary_header
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(long_rows)

        default_metrics = [
            "common_het",
            "assessed_pairs",
            "phase_match_pairs",
            "switch_errors",
            "switch_rate",
            "blockwise_hamming",
            "blockwise_hamming_rate",
        ]
        missing_metrics = [metric for metric in default_metrics if metric not in summary_header]
        if missing_metrics:
            raise SystemExit(
                "phase_compare summary missing expected metric(s): "
                + ", ".join(missing_metrics)
            )
        write_method_summary(work_dir / "summary_by_method.tsv", long_rows, default_metrics)
        replace_output_dir(work_dir, out_dir)
    finally:
        if work_dir.exists():
            remove_path(work_dir)
    print(f"wrote {len(pair_plan)} pairwise comparison(s) to {args.out_dir}")


if __name__ == "__main__":
    main()
