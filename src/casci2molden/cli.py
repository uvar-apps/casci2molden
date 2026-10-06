"""Command-line interface for casci2molden."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .converter import convert_tree


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="casci2molden",
        description=(
            "Convert PySCF CASCI.log files (with sibling .xyz and RHF/ROHF.log) "
            "into Molden files for EigenVista."
        ),
        epilog=(
            "Each CASCI.log needs a *.xyz (geometry) and a RHF.log or ROHF.log "
            "(charge/spin) in the same directory or its parent. Output names are "
            "the log's directory path relative to ROOT, joined with '_'."
        ),
    )
    ap.add_argument(
        "root",
        nargs="?",
        default=".",
        help="directory to scan recursively for CASCI.log files (default: cwd)",
    )
    ap.add_argument(
        "-o",
        "--out",
        default="molden",
        help="output directory for the .molden files (default: ./molden)",
    )
    ap.add_argument(
        "--pattern",
        default="CASCI.log",
        help="filename pattern to scan for (default: CASCI.log)",
    )
    ap.add_argument(
        "--strip",
        action="append",
        default=[],
        metavar="COMPONENT",
        help="drop this path component from output names (repeatable)",
    )
    ap.add_argument(
        "--report",
        default=None,
        metavar="PATH",
        help="path for the JSON report (default: OUT/conversion_report.json)",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="list what would be written without writing anything",
    )
    ap.add_argument(
        "--version",
        action="version",
        version=f"casci2molden {__version__}",
    )
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    root = Path(args.root)
    if not root.exists():
        print(f"root does not exist: {root}", file=sys.stderr)
        return 2

    reports, failures = convert_tree(
        root,
        Path(args.out),
        pattern=args.pattern,
        strip=tuple(args.strip),
        dry=args.dry_run,
    )

    if not reports and not failures:
        print(f"No {args.pattern} files found under {root}", file=sys.stderr)
        return 1

    if not args.dry_run:
        report_path = (
            Path(args.report)
            if args.report
            else Path(args.out) / "conversion_report.json"
        )
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(reports, indent=2) + "\n")
        print(f"Wrote {len(reports)} molden file(s); report: {report_path}")
    else:
        print(f"[dry] {len(reports)} file(s) would be written")

    for src, err in failures:
        print(f"FAILED {src}: {err}", file=sys.stderr)

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
