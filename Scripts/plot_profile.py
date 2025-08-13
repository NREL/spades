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
    mean_data_lst = []
    renames = {
        "::encoded_sort::": "::sort::",
        "::nonencoded_sort::": "::sort::",
    }
    types = ["message", "processed_message", "conjugate_message", "undefined_message"]
    n_mean_steps = 100
    for lbl, fname in zip(args.labels, args.fnames):
        lst = []
        with open(fname, "r") as f:
            logging = False
            for line in f:
                if "MPI processes" in line and "MPI initialized with" in line:
                    nranks = int(line.split()[3])
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

        # Temporal data
        dname = fname.replace("phold", "data").replace("log", "csv")
        data_df = pd.read_csv(dname)
        mean_data_df = data_df.tail(n_mean_steps).mean()
        mean_data_df["final_gvt"] = data_df.gvt.iloc[-1]
        mean_data_df["final_step"] = data_df.step.iloc[-1]
        mean_data_df["nranks"] = nranks
        mean_data_df["function"] = f"average-{fname}"
        mean_data_lst.append(mean_data_df)

        plt.figure("gvt")
        plt.plot(data_df.step, data_df.gvt, label=lbl)

        plt.figure("rate")
        plt.plot(data_df.step, data_df.avg_rate, label=lbl)

        plt.figure(f"rlbk")
        for nrlbk in range(1, 20):
            if data_df[f"rollback_{nrlbk}"].sum() > 0:
                plt.semilogy(
                    data_df.step, data_df[f"rollback_{nrlbk}"], label=f"nrlbk = {nrlbk}"
                )

        for typ in types:
            plt.figure(f"{typ}s")
            plt.plot(data_df.step, data_df[f"{typ}s"] / data_df.lps)
            plt.fill_between(
                data_df.step, data_df[f"min_{typ}s"], data_df[f"max_{typ}s"], alpha=0.5
            )

    mean_data_df = pd.DataFrame(mean_data_lst)
    print(mean_data_df)

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
        "ParticleContainer::RedistributeMPI",
        "ParticleContainer::RedistributeCPU",
    ]
    grouped_data = {
        "function": [
            "communication",
            "computation",
            "total",
            "nranks",
            "nentities",
            "rate",
            "gvt",
            "avg_time"
        ]
    }
    for col in df.select_dtypes(include="number").columns:
        communication_sum = fdf[fdf["function"].isin(communication_functions)][
            col
        ].sum()
        computation_sum = fdf[~fdf["function"].isin(communication_functions)][col].sum()
        nranks = mean_data_df.loc[mean_data_df["function"] == col, "nranks"].iloc[0]
        nentities = mean_data_df.loc[mean_data_df["function"] == col, "entities"].iloc[
            0
        ]
        rate = mean_data_df.loc[mean_data_df["function"] == col, "avg_rate"].iloc[0]
        gvt = mean_data_df.loc[mean_data_df["function"] == col, "final_gvt"].iloc[0]
        avg_time = mean_data_df.loc[mean_data_df["function"] == col, "avg_time"].iloc[0]
        grouped_data[col] = [
            communication_sum,
            computation_sum,
            communication_sum + computation_sum,
            nranks,
            nentities,
            rate,
            gvt,
            avg_time,
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
    print(grouped_df)

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

    plt.figure("scaling-time-gvt")
    markers = itertools.cycle(marker_shapes)
    plt.semilogx(
        grouped_df.nranks,
        grouped_df.communication / grouped_df.gvt,
        label="Communication",
        marker=next(markers),
    )
    plt.semilogx(
        grouped_df.nranks,
        grouped_df.computation / grouped_df.gvt,
        label="Computation",
        marker=next(markers),
    )
    plt.semilogx(
        grouped_df.nranks,
        grouped_df.total / grouped_df.gvt,
        label="Total",
        marker=next(markers),
    )

    plt.figure("scaling-efficiency")
    markers = itertools.cycle(marker_shapes)
    plt.semilogx(
        grouped_df.nranks,
        grouped_df.communication.iloc[0] / grouped_df.communication * 100,
        label="Communication",
        marker=next(markers),
    )
    plt.semilogx(
        grouped_df.nranks,
        grouped_df.computation.iloc[0] / grouped_df.computation * 100,
        label="Computation",
        marker=next(markers),
    )
    plt.semilogx(
        grouped_df.nranks,
        grouped_df.total.iloc[0] / grouped_df.total * 100,
        label="Total",
        marker=next(markers),
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
    markers = itertools.cycle(marker_shapes)
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
        zorder=0,
    )

    plt.figure("scaling-entities-gvt")
    markers = itertools.cycle(marker_shapes)
    plt.loglog(
        grouped_df.entities_per_rank,
        grouped_df.communication / grouped_df.gvt,
        label="Communication",
        marker=next(markers),
    )
    plt.loglog(
        grouped_df.entities_per_rank,
        grouped_df.computation / grouped_df.gvt,
        label="Computation",
        marker=next(markers),
    )
    plt.loglog(
        grouped_df.entities_per_rank,
        grouped_df.total / grouped_df.gvt,
        label="Total",
        marker=next(markers),
    )
    plt.loglog(
        grouped_df.entities_per_rank,
        grouped_df.theory_entities / grouped_df.gvt,
        label="Perfect scaling",
        color="k",
        ls="-",
        zorder=0,
    )

    plt.figure("scaling-rate-efficiency")
    markers = itertools.cycle(marker_shapes)
    plt.semilogx(
        grouped_df.nranks,
        # grouped_df.rate / grouped_df.nentities / ( grouped_df.rate / grouped_df.nentities).iloc[0],
        (grouped_df.rate / grouped_df.nentities)
        / (grouped_df.rate / grouped_df.nentities).iloc[0]
        * 100,
        label="Rate",
        marker=next(markers),
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

        plt.figure("scaling-time-gvt")
        plt.xlabel(f"Number of ranks")
        plt.ylabel(r"$t / g~[-]$")
        legend = plt.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("scaling-efficiency")
        plt.axhline(y=100, color="k", ls="-", label="Perfect scaling", zorder=0)
        plt.xlabel(f"Number of ranks")
        plt.ylabel(r"Parallel efficiency$~[\%]$")
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

        plt.figure("scaling-entities-gvt")
        plt.xlabel(f"Entities per rank")
        plt.ylabel(r"$t / g~[-]$")
        legend = plt.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("scaling-rate-efficiency")
        plt.axhline(y=100, color="k", ls="-", label="Perfect scaling", zorder=0)
        plt.xlabel(f"Number of ranks")
        plt.ylabel(r"Parallel efficiency$~[\%]$")
        legend = plt.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

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
        # plt.ylim([0, 3e6])
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure(f"rlbk")
        plt.xlabel(r"$s~[-]$")
        plt.ylabel(r"$n~[\#]$")
        legend = plt.legend(loc="upper right")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        for typ in types:
            plt.figure(f"{typ}s")
            plt.xlabel(r"$s~[-]$")
            plt.ylabel(f"$m_{typ[0]}~[\\#]$")
            plt.tight_layout()
            pdf.savefig(dpi=300)
