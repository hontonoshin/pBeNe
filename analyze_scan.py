#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


LOG_DIR = Path("logs_scan")
OUT_DIR = Path("analysis_scan")
OUT_DIR.mkdir(exist_ok=True)


CHANNEL_LABELS = {
    "B11_pn": "11B(p,n)11C",
    "B10_pn": "10B(p,n)10C",
    "N14_pn": "14N(p,n)14O",
    "B11_an": "11B(a,n)14N",
    "B10_an": "10B(a,n)13N",
    "N14_an": "14N(a,n)17F",
}


def parse_energy_from_filename(name: str) -> float:
    m = re.search(r"pBN_([0-9]+(?:p[0-9]+)?)MeV", name)
    if not m:
        raise ValueError(f"Could not parse energy from filename: {name}")
    return float(m.group(1).replace("p", "."))


def parse_int(s: str) -> int:
    return int(s.replace(",", "").strip())


def clean_line(line: str) -> str:
    # Remove box-drawing characters and normalize spaces.
    for ch in "║╠╚╔═╗╝":
        line = line.replace(ch, " ")
    return " ".join(line.split())


def parse_channel_line(line: str, label: str):
    """
    Robustly parse a channel row after box characters are removed.

    Expected cleaned examples:
        11B(p,n)11C 825 8.250e-04 78.7%
        11B(a,n)14N 0 0.000e+00 0.0%
    """
    line = clean_line(line)
    if label not in line:
        return None

    rest = line.split(label, 1)[1].strip()
    m = re.search(r"([0-9,]+)\s+([0-9.eE+\-]+)\s+([0-9.eE+\-]+)%", rest)
    if not m:
        return None

    return {
        "N": parse_int(m.group(1)),
        "yield": float(m.group(2)),
        "frac_percent": float(m.group(3)),
    }


def parse_total_line(text: str, label: str):
    for raw in text.splitlines():
        line = clean_line(raw)
        if label not in line:
            continue

        rest = line.split(label, 1)[1].strip()
        m = re.search(r"([0-9,]+)\s+([0-9.eE+\-]+)\s+([0-9.eE+\-]+)%", rest)
        if m:
            return {
                "N": parse_int(m.group(1)),
                "yield": float(m.group(2)),
                "frac_percent": float(m.group(3)),
            }
    return None


def parse_log(path: Path) -> dict:
    text = path.read_text(errors="ignore")

    result = {
        "energy_MeV": parse_energy_from_filename(path.name),
        "log_file": str(path),
        "primary_protons": np.nan,
        "total_neutrons": 0,
        "alpha_secondaries": 0,
    }

    m = re.search(r"Primary protons\s*:\s*([0-9,]+)", text)
    if m:
        result["primary_protons"] = parse_int(m.group(1))

    m = re.search(r"Total neutrons\s*:\s*([0-9,]+)", text)
    if m:
        result["total_neutrons"] = parse_int(m.group(1))

    m = re.search(r"Alpha secondaries\s*:\s*([0-9,]+)", text)
    if m:
        result["alpha_secondaries"] = parse_int(m.group(1))

    # Initialize channel columns.
    for key in CHANNEL_LABELS:
        result[f"{key}_N"] = 0
        result[f"{key}_yield"] = 0.0
        result[f"{key}_frac_percent"] = 0.0

    # Parse channel rows line-by-line.
    for raw in text.splitlines():
        for key, label in CHANNEL_LABELS.items():
            parsed = parse_channel_line(raw, label)
            if parsed is not None:
                result[f"{key}_N"] = parsed["N"]
                result[f"{key}_yield"] = parsed["yield"]
                result[f"{key}_frac_percent"] = parsed["frac_percent"]

    # Parse group totals.
    pn_total = parse_total_line(text, "Total p,n")
    an_total = parse_total_line(text, "Total alpha,n")

    if pn_total is None:
        pn_N = result["B11_pn_N"] + result["B10_pn_N"] + result["N14_pn_N"]
        pn_yield = pn_N / result["primary_protons"] if result["primary_protons"] else 0.0
        pn_frac = 100.0 * pn_N / result["total_neutrons"] if result["total_neutrons"] else 0.0
    else:
        pn_N, pn_yield, pn_frac = pn_total["N"], pn_total["yield"], pn_total["frac_percent"]

    if an_total is None:
        an_N = result["B11_an_N"] + result["B10_an_N"] + result["N14_an_N"]
        an_yield = an_N / result["primary_protons"] if result["primary_protons"] else 0.0
        an_frac = 100.0 * an_N / result["total_neutrons"] if result["total_neutrons"] else 0.0
    else:
        an_N, an_yield, an_frac = an_total["N"], an_total["yield"], an_total["frac_percent"]

    result["total_pn_N"] = pn_N
    result["total_pn_yield"] = pn_yield
    result["total_pn_frac_percent"] = pn_frac

    result["total_an_N"] = an_N
    result["total_an_yield"] = an_yield
    result["total_an_frac_percent"] = an_frac

    total = result["total_neutrons"]
    result["N14_pn_fraction_of_total_percent"] = 100.0 * result["N14_pn_N"] / total if total else 0.0
    result["N14_pn_fraction_of_pn_percent"] = 100.0 * result["N14_pn_N"] / pn_N if pn_N else 0.0
    result["alpha_n_fraction_of_total_percent"] = 100.0 * an_N / total if total else 0.0

    return result


def save_plot(df: pd.DataFrame, ycols: list[str], ylabel: str, filename: str, *, logy: bool = False) -> None:
    plt.figure(figsize=(8, 5))

    for col in ycols:
        if col not in df.columns:
            continue
        plt.plot(df["energy_MeV"], df[col], marker="o", label=col)

    plt.xlabel("Proton energy [MeV]")
    plt.ylabel(ylabel)
    plt.grid(True, alpha=0.3)
    if logy:
        plt.yscale("log")
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT_DIR / filename, dpi=300)
    plt.close()


def main() -> None:
    logs = sorted(LOG_DIR.glob("pBN_*MeV_*p.log"))
    if not logs:
        raise SystemExit(f"No log files found in {LOG_DIR.resolve()}")

    df = pd.DataFrame(parse_log(p) for p in logs)
    df = df.sort_values("energy_MeV").reset_index(drop=True)

    csv_path = OUT_DIR / "scan_summary.csv"
    df.to_csv(csv_path, index=False)

    cols = [
        "energy_MeV",
        "primary_protons",
        "total_neutrons",
        "total_pn_N",
        "total_an_N",
        "B11_pn_N",
        "B10_pn_N",
        "N14_pn_N",
        "B11_an_N",
        "B10_an_N",
        "N14_an_N",
        "N14_pn_fraction_of_total_percent",
        "alpha_n_fraction_of_total_percent",
    ]
    print(df[cols].to_string(index=False))

    save_plot(
        df,
        ["total_pn_yield", "total_an_yield"],
        "Yield [neutrons / primary proton]",
        "yield_pn_vs_alpha_n.png",
    )

    save_plot(
        df,
        ["B11_pn_yield", "B10_pn_yield", "N14_pn_yield"],
        "Yield [neutrons / primary proton]",
        "primary_pn_channel_yields.png",
    )

    save_plot(
        df,
        ["B11_pn_yield", "B10_pn_yield", "N14_pn_yield"],
        "Yield [neutrons / primary proton]",
        "primary_pn_channel_yields_log.png",
        logy=True,
    )

    save_plot(
        df,
        ["B11_an_yield", "B10_an_yield", "N14_an_yield"],
        "Yield [neutrons / primary proton]",
        "secondary_alpha_n_channel_yields.png",
    )

    save_plot(
        df,
        ["N14_pn_fraction_of_total_percent", "alpha_n_fraction_of_total_percent"],
        "Fraction of total neutrons [%]",
        "fractions_nitrogen_and_alpha_n.png",
    )

    save_plot(
        df,
        ["N14_pn_fraction_of_pn_percent"],
        "Fraction of primary p,n neutrons [%]",
        "nitrogen_fraction_of_primary_pn.png",
    )

    print()
    print(f"Wrote: {csv_path}")
    print(f"Wrote plots to: {OUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
