"""Plot tiny profile output in a log file."""

import argparse
import pathlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import numpy as np
from cycler import cycler
from functools import reduce
import itertools
import re

plt.style.use(pathlib.Path(__file__).parent.resolve() / "project.mplstyle")
marker_shapes = ("s", "d", "o", "p", "h")
markers = itertools.cycle(marker_shapes)

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
    extra_dct = {}
    entity_regex = r"^  \d+ entities$"
    for fname in args.fnames:
        lst = []
        with open(fname, "r") as f:
            logging = False
            for line in f:
                if "MPI processes" in line and "MPI initialized with" in line:
                    nranks = int(line.split()[3])
                if re.match(entity_regex, line):
                    nentities = int(line.split()[0])
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
        extra_dct[f"average-{fname}"] = {"nranks": nranks, "nentities": nentities}

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
    df.function = df.function.str.replace("_", "\\_", regex=False)

    # Group communication and computation functions and make a new df
    excluded_function = "spades::Total"
    fdf = df[df["function"] != excluded_function]
    communication_functions = [
        "amrex::communicateParticlesFinish",
        "ParticleCopyPlan::doHandShake",
        "Redistribute_partition",
        "amrex::communicateParticlesStart",
        "amrex::unpackRemotes",
        "ParticleCopyPlan::buildMPIStart",
        "ParticleCopyPlan::build",
    ]
    grouped_data = {
        "function": ["communication", "computation", "total", "nranks", "nentities"]
    }
    for col in df.select_dtypes(include="number").columns:
        communication_sum = fdf[fdf["function"].isin(communication_functions)][
            col
        ].sum()
        computation_sum = fdf[~fdf["function"].isin(communication_functions)][col].sum()
        grouped_data[col] = [
            communication_sum,
            computation_sum,
            communication_sum + computation_sum,
            extra_dct[col]["nranks"],
            extra_dct[col]["nentities"],
        ]
    grouped_df = pd.DataFrame(grouped_data).set_index("function").T

    norm = grouped_df.iloc[0]
    norm_cols = ["total", "communication", "computation"]
    for col in norm_cols:
        grouped_df[f"norm-{col}"] = grouped_df[col] / norm[col]
    grouped_df["entities_per_rank"] = grouped_df.nentities / grouped_df.nranks
    theory_idx = 0
    grouped_df["theory_entities"] = (
        grouped_df.total.iloc[theory_idx]
        * grouped_df.entities_per_rank
        / grouped_df.entities_per_rank.iloc[theory_idx]
    )

    # sort and keep the top consuming functions for the original df
    df["average"] = df[[x for x in df.columns if "function" not in x]].mean(axis=1)
    df = df.sort_values(by="average", ascending=False).head(8)

    # normalize the data by the first entry
    norm = df.loc[df.function == "spades::Total"][f"average-{args.fnames[0]}"].values[0]
    norm_cols = df.columns.difference(["function"])
    for col in norm_cols:
        df[f"norm-{col}"] = df[col] / norm

    # print_fname = "spades::Total"
    print_fname = "spades::SpadesParticleContainer::sort::sort"
    print(df[df.function == print_fname].T)

    pname = "profile_plots.pdf"
    plt.figure("timing", figsize=(14, 6))
    ax = plt.gca()
    ind = np.arange(len(df))
    width = 0.8 / (len(args.fnames))
    offset = 0.5 * (len(args.fnames) - 1) * width
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

    plt.figure("norm-timing", figsize=(14, 6))
    ax = plt.gca()
    ind = np.arange(len(df))
    width = 0.8 / (len(args.fnames))
    offset = 0.5 * (len(args.fnames) - 1) * width
    for k, fname in enumerate(args.fnames):
        ax.barh(
            ind - offset + k * width,
            df[f"norm-average-{fname}"],
            width,
            align="center",
            label=args.labels[k],
        )
    ax.set(yticks=ind, yticklabels=df.function, ylim=[2 * width - 1, len(df)])
    ax.invert_yaxis()

    plt.figure("scaling")
    plt.semilogx(
        grouped_df.nranks,
        grouped_df.communication,
        label="Communication",
        marker=next(markers),
    )
    plt.semilogx(
        grouped_df.nranks,
        grouped_df.computation,
        label="Computation",
        marker=next(markers),
    )
    plt.semilogx(
        grouped_df.nranks, grouped_df.total, label="Total", marker=next(markers)
    )

    plt.figure("norm-scaling")
    markers = itertools.cycle(marker_shapes)
    plt.semilogx(
        grouped_df.nranks,
        grouped_df["norm-communication"],
        label="Communication",
        marker=next(markers),
    )
    plt.semilogx(
        grouped_df.nranks,
        grouped_df["norm-computation"],
        label="Computation",
        marker=next(markers),
    )
    plt.semilogx(
        grouped_df.nranks, grouped_df["norm-total"], label="Total", marker=next(markers)
    )

    plt.figure("scaling-entities")
    plt.loglog(
        grouped_df.entities_per_rank,
        grouped_df.communication,
        label="Communication",
        marker=next(markers),
    )
    plt.loglog(
        grouped_df.entities_per_rank,
        grouped_df.computation,
        label="Computation",
        marker=next(markers),
    )
    plt.loglog(
        grouped_df.entities_per_rank,
        grouped_df.total,
        label="Total",
        marker=next(markers),
    )
    plt.loglog(
        grouped_df.entities_per_rank,
        grouped_df.theory_entities,
        label="Perfect scaling",
        color="k",
        ls="-",
    )

    # Save the plots
    with PdfPages(pname) as pdf:
        plt.figure("timing")
        ax = plt.gca()
        plt.xlabel(r"$t~[s]$", fontsize=22, fontweight="bold")
        # plt.ylabel(r"$\bar{u} (x=0)$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=22, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=22, fontweight="bold")
        legend = ax.legend(loc="lower right", fontsize=22)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("norm-timing")
        ax = plt.gca()
        plt.xlabel(r"$t / \tau~[-]$", fontsize=22, fontweight="bold")
        # plt.ylabel(r"$\bar{u} (x=0)$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=22, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=22, fontweight="bold")
        legend = ax.legend(loc="lower right", fontsize=22)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("scaling")
        plt.xlabel(f"Number of ranks")
        plt.ylabel(r"$t~[s]$")
        legend = plt.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("norm-scaling")
        plt.xlabel(f"Number of ranks")
        plt.ylabel(r"$t~[-]$")
        legend = plt.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("scaling-entities")
        plt.xlabel(f"Entities per rank")
        plt.ylabel(r"$t~[s]$")
        legend = plt.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)
