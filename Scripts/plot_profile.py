"""Plot tiny profile output in a log file."""

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
    renames = {
        "::encoded_sort::": "::sort::",
        "::nonencoded_sort::": "::sort::",
    }
    for fname in args.fnames:
        lst = []
        with open(fname, "r") as f:
            logging = False
            for line in f:
                if line.startswith("Time spent in evolve():"):
                    total_time = float(line.split()[-1])
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

        lst.append(
            {
                "function": "spades::Total",
                f"average-{fname}": total_time,
            }
        )
        df = pd.DataFrame(lst)
        df.function = df.function.apply(
            lambda x: reduce(lambda s, kv: s.replace(*kv), renames.items(), x)
        )
        df_lst.append(df)

    # Make sure all the df have the same reported functions, fill missing ones with 0
    all_functions = list(set([x for df in df_lst for x in df["function"].to_list()]))
    for df in df_lst:
        for x in all_functions:
            if x not in df.function.values:
                dct = {col: 0.0 for col in df.columns if "function" not in col}
                dct["function"] = x
                df.loc[len(df)] = dct

    # Make a single df
    df = reduce(
        lambda left, right: pd.merge(left, right, on=["function"], how="inner"), df_lst
    )
    df["average"] = df[[x for x in df.columns if "function" not in x]].mean(axis=1)
    df = df.sort_values(by="average", ascending=False).head(10)
    df.function = df.function.str.replace("_", "\_", regex=False)

    pname = "profile_plots.pdf"
    plt.figure("timing", figsize=(14, 6))
    ax = plt.gca()
    ind = np.arange(len(df))
    width = 0.8 / (len(args.fnames))
    offset = 0.5* (len(args.fnames)-1) * width
    for k, fname in enumerate(args.fnames):
        ax.barh(
            ind - offset + k * width,
            df[f"average-{fname}"],
            width,
            align="center",
            label=args.labels[k],
        )
    ax.set(yticks=ind, yticklabels=df.function, ylim=[2 * width - 1, len(df)])
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
