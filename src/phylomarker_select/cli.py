"""Interfaz de linea de ordenes de phylomarker-select."""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from .busco import discover_busco_runs, validate_runs
from .metadata import load_metadata
from .pipeline import run_pipeline

LOGGER = logging.getLogger("phylomarker-select")


def configure_logging(verbose: bool) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )


def validate_command(
    args: argparse.Namespace,
) -> None:
    metadata = load_metadata(
        Path(args.metadata),
        args.sample_id_column,
    )

    runs = discover_busco_runs(
        Path(args.busco_directory)
    )

    validated, warnings = validate_runs(
        runs,
        metadata,
        args.sample_id_column,
    )

    output_directory = Path(
        args.output
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    runs.to_csv(
        output_directory
        / "discovered_runs.tsv",
        sep="\t",
        index=False,
    )

    validated.to_csv(
        output_directory
        / "validated_runs.tsv",
        sep="\t",
        index=False,
    )

    warnings.to_csv(
        output_directory
        / "warnings.tsv",
        sep="\t",
        index=False,
    )

    print(
        f"Detected BUSCO runs: {len(runs)}"
    )

    print(
        f"Validated BUSCO runs: {len(validated)}"
    )

    print(
        f"Warnings or errors: {len(warnings)}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="phylomarker-select",
        description=(
            "Evolution-aware selection of "
            "phylogenomic marker panels from BUSCO outputs."
        ),
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable detailed logging.",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    run_parser = subparsers.add_parser(
        "run",
        help="Run the complete workflow.",
    )

    run_parser.add_argument(
        "--config",
        required=True,
        type=Path,
        help="YAML configuration file.",
    )

    validate_parser = subparsers.add_parser(
        "validate",
        help="Validate BUSCO runs and metadata.",
    )

    validate_parser.add_argument(
        "--busco-directory",
        required=True,
    )

    validate_parser.add_argument(
        "--metadata",
        required=True,
    )

    validate_parser.add_argument(
        "--sample-id-column",
        default="sample_ID",
    )

    validate_parser.add_argument(
        "--output",
        default="validation_results",
    )

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    configure_logging(
        args.verbose
    )

    try:
        if args.command == "run":
            run_pipeline(
                args.config
            )
        elif args.command == "validate":
            validate_command(
                args
            )
        else:
            parser.error(
                f"Unsupported command: {args.command}"
            )
    except Exception as error:
        LOGGER.error("%s", error)

        if args.verbose:
            raise

        sys.exit(1)


if __name__ == "__main__":
    main()
