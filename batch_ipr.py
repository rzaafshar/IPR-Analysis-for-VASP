#!/usr/bin/env python3

from pathlib import Path
import argparse
import csv
import gc
import re
import sys

import numpy as np
from pymatgen.io.vasp.outputs import Chgcar


def normalized_ipr(rho):
    """
    Dimensionless normalized real-space IPR:

        IPR* = N_grid * sum(rho^2) / sum(rho)^2
    """
    rho = np.asarray(rho, dtype=np.float64)

    density_sum = np.sum(rho)

    if abs(density_sum) < 1.0e-14:
        return np.nan

    return float(
        rho.size
        * np.sum(rho**2)
        / density_sum**2
    )


def read_ibands(incar_path):
    """Read IBAND values from INCAR."""
    if not incar_path.is_file():
        return []

    text = incar_path.read_text(
        encoding="utf-8",
        errors="ignore",
    )

    match = re.search(
        r"(?im)^\s*IBAND\s*=\s*([^!#\n]+)",
        text,
    )

    if match is None:
        return []

    return [
        int(number)
        for number in re.findall(
            r"\d+",
            match.group(1),
        )
    ]


def determine_band(parchg_path, ibands):
    """
    Determine the band number from the PARCHG filename
    or from IBAND in INCAR.
    """
    suffix_numbers = []

    for token in parchg_path.name.split(".")[1:]:
        if token.isdigit():
            suffix_numbers.append(int(token))

    # Prefer a filename number that also occurs in IBAND.
    for number in suffix_numbers:
        if number in ibands:
            return str(number)

    # For an exact file named PARCHG with one IBAND.
    if len(ibands) == 1:
        return str(ibands[0])

    # Multiple bands combined in one PARCHG.
    if len(ibands) > 1:
        return ";".join(str(number) for number in ibands)

    # Last fallback: first numeric filename suffix.
    if suffix_numbers:
        return str(suffix_numbers[0])

    return ""


def format_number(value):
    if value is None:
        return ""

    if not np.isfinite(value):
        return ""

    return f"{value:.8f}"


def find_parchg_files(root):
    """
    Find PARCHG and PARCHG.* files only inside
    directories named 'parch'.
    """
    files = []

    for path in root.rglob("PARCHG*"):
        if not path.is_file():
            continue

        if path.parent.name != "parch":
            continue

        files.append(path)

    return sorted(files)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate spin-resolved normalized real-space "
            "IPRs from all existing VASP PARCHG files."
        )
    )

    parser.add_argument(
        "root",
        nargs="?",
        default=".",
        help="Root directory to search. Default: current directory",
    )

    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output CSV filename",
    )

    args = parser.parse_args()

    root = Path(args.root).resolve()

    if args.output is None:
        output_file = root / "ipr_results.csv"
    else:
        output_file = Path(args.output).resolve()

    parchg_files = find_parchg_files(root)

    if not parchg_files:
        print(
            "ERROR: No PARCHG files were found inside "
            "directories named 'parch'.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Found {len(parchg_files)} PARCHG files.")

    columns = [
        "material",
        "defect",
        "parchg_file",
        "selected_band",
        "grid_x",
        "grid_y",
        "grid_z",
        "n_grid",
        "ipr_total",
        "ipr_spin_up",
        "ipr_spin_down",
        "ipr_max",
        "higher_ipr_spin",
        "spin_ipr_ratio",
        "sum_spin_up",
        "sum_spin_down",
        "status",
    ]

    with output_file.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as csv_handle:

        writer = csv.DictWriter(
            csv_handle,
            fieldnames=columns,
        )

        writer.writeheader()

        for index, parchg_path in enumerate(
            parchg_files,
            start=1,
        ):
            parch_directory = parchg_path.parent
            defect_directory = parch_directory.parent

            defect_name = defect_directory.name
            material_name = defect_directory.parent.name

            ibands = read_ibands(
                parch_directory / "INCAR"
            )

            selected_band = determine_band(
                parchg_path,
                ibands,
            )

            relative_path = str(
                parchg_path.relative_to(root)
            )

            print(
                f"[{index}/{len(parchg_files)}] "
                f"{material_name}/{defect_name}/"
                f"{parchg_path.name}"
            )

            row = {
                column: ""
                for column in columns
            }

            row.update(
                {
                    "material": material_name,
                    "defect": defect_name,
                    "parchg_file": relative_path,
                    "selected_band": selected_band,
                }
            )

            rho_total = None
            rho_magnetic = None
            rho_up = None
            rho_down = None
            parchg = None

            try:
                parchg = Chgcar.from_file(
                    str(parchg_path)
                )

                rho_total = np.asarray(
                    parchg.data["total"],
                    dtype=np.float64,
                )

                if rho_total.ndim != 3:
                    raise ValueError(
                        "The PARCHG grid is not three-dimensional."
                    )

                grid_x, grid_y, grid_z = (
                    rho_total.shape
                )

                ipr_total = normalized_ipr(
                    rho_total
                )

                row.update(
                    {
                        "grid_x": grid_x,
                        "grid_y": grid_y,
                        "grid_z": grid_z,
                        "n_grid": rho_total.size,
                        "ipr_total": format_number(
                            ipr_total
                        ),
                    }
                )

                if "diff" in parchg.data:
                    rho_magnetic = np.asarray(
                        parchg.data["diff"],
                        dtype=np.float64,
                    )

                    if rho_magnetic.shape != rho_total.shape:
                        raise ValueError(
                            "Total and magnetic grids "
                            "have different dimensions."
                        )

                    rho_up = 0.5 * (
                        rho_total + rho_magnetic
                    )

                    rho_down = 0.5 * (
                        rho_total - rho_magnetic
                    )

                    ipr_up = normalized_ipr(rho_up)
                    ipr_down = normalized_ipr(rho_down)

                    if ipr_up >= ipr_down:
                        higher_spin = "up"
                        ipr_max = ipr_up
                    else:
                        higher_spin = "down"
                        ipr_max = ipr_down

                    smaller_ipr = min(
                        ipr_up,
                        ipr_down,
                    )

                    if smaller_ipr > 0:
                        spin_ratio = (
                            ipr_max / smaller_ipr
                        )
                    else:
                        spin_ratio = np.nan

                    row.update(
                        {
                            "ipr_spin_up": format_number(
                                ipr_up
                            ),
                            "ipr_spin_down": format_number(
                                ipr_down
                            ),
                            "ipr_max": format_number(
                                ipr_max
                            ),
                            "higher_ipr_spin": higher_spin,
                            "spin_ipr_ratio": format_number(
                                spin_ratio
                            ),
                            "sum_spin_up": format_number(
                                np.sum(rho_up)
                            ),
                            "sum_spin_down": format_number(
                                np.sum(rho_down)
                            ),
                        }
                    )

                else:
                    row.update(
                        {
                            "ipr_max": format_number(
                                ipr_total
                            ),
                            "higher_ipr_spin": (
                                "not_spin_resolved"
                            ),
                        }
                    )

                row["status"] = "OK"

            except Exception as error:
                row["status"] = (
                    f"ERROR: {type(error).__name__}: {error}"
                )

                print(
                    f"    ERROR: {error}",
                    file=sys.stderr,
                )

            writer.writerow(row)
            csv_handle.flush()

            del parchg
            del rho_total
            del rho_magnetic
            del rho_up
            del rho_down

            gc.collect()

    print()
    print(f"Completed: {output_file}")


if __name__ == "__main__":
    main()
