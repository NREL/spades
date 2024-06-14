"""Plot tiny profile output in a log file."""

import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
from cycler import cycler
from functools import reduce

default_cycler = cycler(
    color=[
        "#EE2E2F",
        "#008C48",
        "#185AA9",
        "#F47D23",
        "#662C91",
        "#A21D21",
        "#B43894",
        "#010202",
    ]
) + cycler(
    linestyle=[
        "-",
        "--",
        ":",
        "-.",
        (0, (1, 10)),
        (0, (1, 1)),
        (0, (5, 10)),
        (0, (5, 1)),
    ]
)

plt.rc("text", usetex=True)
plt.rc("lines", linewidth=2)
plt.rc("axes", prop_cycle=default_cycler)

if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Plot tiny profile output in a log file"
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

    df_lst = []
    for fname in args.fnames:
        lst = []
        with open(fname, "r") as f:
            logging = False
            for line in f:
                if line.startswith("Name") and ("NCalls" in line):
                    line = next(f)
                    line = next(f)
                    logging = True
                if logging and line.startswith("-----------"):
                    logging = False
                    break
                if logging:
                    line = line.split()
                    lst.append(
                        {
                            "function": line[0].replace("(", "").replace(")", ""),
                            f"average-{fname}": float(line[3]),
                        }
                    )

        df_lst.append(pd.DataFrame(lst))
    df = reduce(
        lambda left, right: pd.merge(left, right, on=["function"], how="outer"), df_lst
    )
    df = df.sort_values(by=f"average-{args.fnames[0]}", ascending=False).head(10)
    df.function = df.function.str.replace("_", "\_", regex=False)

    pname = "plots.pdf"
    plt.figure("timing", figsize=(14, 6))
    ax = plt.gca()
    ind = np.arange(len(df))
    width = 0.25
    for k, fname in enumerate(args.fnames):
        ax.barh(
            ind + k * width,
            df[f"average-{fname}"],
            width,
            align="center",
            label=args.labels[k],
        )
    ax.set(yticks=ind + width, yticklabels=df.function, ylim=[2 * width - 1, len(df)])
    ax.invert_yaxis()

    # Save the plots
    with PdfPages(pname) as pdf:
        plt.figure("timing")
        plt.xlabel(r"Time $[s]$", fontsize=22, fontweight="bold")
        # plt.ylabel(r"$\bar{u} (x=0)$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=22, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=22, fontweight="bold")
        legend = ax.legend(loc="lower right", fontsize=22)
        plt.tight_layout()
        pdf.savefig(dpi=300)
