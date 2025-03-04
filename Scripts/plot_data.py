"""Plot simulation data output."""

import argparse
import pathlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
from cycler import cycler
from functools import reduce

plt.style.use(pathlib.Path(__file__).parent.resolve() / "project.mplstyle")

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Plot simulation data output"
    )
    parser.add_argument(
        "-f",
        "--fnames",
        help="Files to plot",
        required=True,
        nargs="+",
        type=str,
    )
    parser.add_argument(
        "-l",
        "--labels",
        help="Labels",
        required=True,
        nargs="+",
        type=str,
    )
    args = parser.parse_args()

    if not len(args.labels) == len(args.fnames):
        raise AssertionError("Need same number of labels and file names")

    types = ["message", "processed_message", "conjugate_message", "undefined_message"]
    for lbl, fname in zip(args.labels, args.fnames):
        df = pd.read_csv(fname)

        plt.figure("gvt")
        plt.plot(df.step, df.gvt,label=lbl)

        plt.figure("rate")
        plt.plot(df.step, df.avg_rate,label=lbl)

        plt.figure(f"rlbk")
        for nrlbk in range(1, 20):
            if df[f"rollback_{nrlbk}"].sum() > 0:
                plt.semilogy(df.step, df[f"rollback_{nrlbk}"], label = f"nrlbk = {nrlbk}")
            
        for typ in types:
            plt.figure(f"{typ}s")
            plt.fill_between(df.step, df[f"min_{typ}s"], df[f"max_{typ}s"], alpha=0.5)

    pname = "profile_data.pdf"
    with PdfPages(pname) as pdf:
        plt.figure("gvt")
        plt.xlabel(r"$s~[-]$")
        plt.ylabel(r"$g~[s]$")
        legend = plt.legend(loc="lower right")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("rate")
        plt.xlabel(r"$s~[-]$")
        plt.ylabel(r"$r~[\#/s]$")
        legend = plt.legend(loc="lower right")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure(f"rlbk")
        plt.xlabel(r"$s~[-]$")
        plt.ylabel(r"$n~[\#]$")
        legend = plt.legend(loc="lower right")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        for typ in types:
            plt.figure(f"{typ}s")
            plt.xlabel(r"$s~[-]$")
            plt.ylabel(r"$m~[\#]$")
            plt.tight_layout()
            pdf.savefig(dpi=300)
