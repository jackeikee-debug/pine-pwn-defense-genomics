#!/usr/bin/env python3
"""Placeholder CLI for building gene family count matrices."""

import argparse


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthogroups", required=True)
    parser.add_argument("--output", required=True)
    parser.parse_args()
    raise SystemExit("Not implemented yet. See project Task 5.")


if __name__ == "__main__":
    main()

