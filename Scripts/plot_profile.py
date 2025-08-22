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
plt.rcParams.update({"figure.max_open_warning": 0})
marker_shapes = ("s", "d", "o", "p", "h")
markers = itertools.cycle(marker_shapes)


def legend_without_duplicates(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique_labels = dict(zip(labels, handles))

    handles = list(unique_labels.values())
    labels = list(unique_labels.keys())

    # Move a label to the end
    last_label = "perfect scaling"
    if last_label in labels:
        last_index = labels.index(last_label)
        handles.append(handles.pop(last_index))
        labels.append(labels.pop(last_index))

    return plt.legend(handles, labels, loc="best")


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
    parser.add_argument(
        "--breakdown",
        action="store_true",
        help="Show communication/computation breakdown",
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
        arch = "Intel CPU node"
        basename = pathlib.Path(fname).name
        with open(fname, "r") as f:
            logging = False
            for line in f:
                if "MPI processes" in line and "MPI initialized with" in line:
                    nranks = int(line.split()[3])
                if line.startswith("Time spent in evolve():"):
                    total_time = float(line.split()[-1])
                if line.startswith("CUDA initialized with"):
                    arch = "NVIDIA GPU"
                if line.startswith("HIP initialized with"):
                    arch = "AMD GPU"
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
        mean_data_df["arch"] = arch
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
            "entities",
            "lps",
            "avg_rate",
            "gvt",
            "avg_time",
            "arch",
            "rollback_1",
            "rollback_2",
            "rollback_3",
            "rollback_4",
        ]
    }
    for col in df.select_dtypes(include="number").columns:
        communication_sum = fdf[fdf["function"].isin(communication_functions)][
            col
        ].sum()
        computation_sum = fdf[~fdf["function"].isin(communication_functions)][col].sum()

        lst = []
        for name in grouped_data["function"]:
            if name not in ["communication", "computation", "total"]:
                lst.append(
                    mean_data_df.loc[mean_data_df["function"] == col, name].iloc[0]
                )
        grouped_data[col] = [
            communication_sum,
            computation_sum,
            communication_sum + computation_sum,
        ] + lst
    grouped_df = pd.DataFrame(grouped_data).set_index("function").T

    norm = grouped_df.iloc[0]
    norm_cols = ["total", "communication", "computation"]
    for col in norm_cols:
        grouped_df[f"norm-{col}"] = grouped_df[col] / norm[col]
    grouped_df["entities_per_rank"] = grouped_df.entities / grouped_df.nranks
    grouped_df["entities_per_lp"] = grouped_df.entities / grouped_df.lps

    # sort and keep the top consuming functions for the original df
    df["average"] = df[[x for x in df.columns if "function" not in x]].mean(axis=1)
    df = df.sort_values(by="average", ascending=False).head(8)

    # get and sort the fnames/labels based on the number of ranks
    fnames_df = pd.DataFrame(
        {
            "fname": grouped_df.index,
            "nranks": grouped_df["nranks"],
            "label": args.labels,
        }
    ).reset_index(drop=True)
    fnames_df["fname"] = fnames_df["fname"].str.replace("average-", "")
    fnames_df.sort_values(by=["nranks"], inplace=True)

    # normalize the data by the first entry
    norm = df.loc[df.function == "spades::Total"][
        f"average-{fnames_df["fname"].iloc[0]}"
    ].values[0]
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
    for k, fname in enumerate(fnames_df["fname"]):
        ax.barh(
            ind - offset + k * width,
            df[f"average-{fname}"],
            width,
            align="center",
            label=fnames_df["label"].iloc[k],
        )
    ax.set(yticks=ind, yticklabels=df.function, ylim=[2 * width - 1, len(df)])
    ax.invert_yaxis()

    plt.figure("norm-timing", figsize=(14, 6))
    ax = plt.gca()
    ind = np.arange(len(df))
    width = 0.8 / (len(args.fnames))
    offset = 0.5 * (len(args.fnames) - 1) * width
    for k, fname in enumerate(fnames_df["fname"]):
        ax.barh(
            ind - offset + k * width,
            df[f"norm-average-{fname}"],
            width,
            align="center",
            label=fnames_df["label"].iloc[k],
        )
    ax.set(yticks=ind, yticklabels=df.function, ylim=[2 * width - 1, len(df)])
    ax.invert_yaxis()

    markers = itertools.cycle(marker_shapes)
    for cnt, arch in enumerate(grouped_df["arch"].unique()):
        mkr = marker_shapes[cnt]
        sg = grouped_df[grouped_df["arch"] == arch].copy()
        theory_idx = 0
        sg.sort_values(by=["nranks"], inplace=True, ignore_index=True, ascending=True)
        sg["theory_ranks"] = (
            sg.total.iloc[theory_idx] * sg.nranks.iloc[theory_idx] / sg.nranks
        )
        sg.sort_values(
            by=["entities_per_rank"], inplace=True, ignore_index=True, ascending=False
        )
        sg["theory_entities"] = (
            sg.total.iloc[theory_idx]
            * sg.entities_per_rank
            / sg.entities_per_rank.iloc[theory_idx]
        )

        plt.figure("scaling")
        sg.sort_values(by=["nranks"], inplace=True)
        if args.breakdown:
            plt.loglog(
                sg.nranks,
                sg.communication,
                label=f"{arch} communication",
                marker=mkr,
            )
            plt.loglog(
                sg.nranks,
                sg.computation,
                label=f"{arch} computation",
                marker=mkr,
            )
        plt.loglog(sg.nranks, sg.total, label=f"{arch}", marker=next(markers))
        plt.loglog(
            sg.nranks,
            sg.theory_ranks,
            label="perfect scaling",
            color="k",
            ls="-",
            zorder=0,
        )

        plt.figure("scaling-time-gvt")
        sg.sort_values(by=["nranks"], inplace=True)
        if args.breakdown:
            plt.loglog(
                sg.nranks,
                sg.gvt / sg.communication,
                label=f"{arch} communication",
                marker=mkr,
            )
            plt.loglog(
                sg.nranks,
                sg.gvt / sg.computation,
                label=f"{arch} computation",
                marker=mkr,
            )
        plt.loglog(
            sg.nranks,
            sg.gvt / sg.total,
            label=f"{arch}",
            marker=mkr,
        )
        plt.loglog(
            sg.nranks,
            sg.gvt / sg.theory_ranks,
            label="perfect scaling",
            color="k",
            ls="-",
            zorder=0,
        )

        plt.figure("scaling-efficiency")
        markers = itertools.cycle(marker_shapes)
        sg.sort_values(by=["nranks"], inplace=True)
        if args.breakdown:
            plt.semilogx(
                sg.nranks,
                sg.communication.iloc[0] / sg.communication * 100,
                label=f"{arch} communication",
                marker=mkr,
            )
            plt.semilogx(
                sg.nranks,
                sg.computation.iloc[0] / sg.computation * 100,
                label=f"{arch} computation",
                marker=mkr,
            )
        plt.semilogx(
            sg.nranks,
            sg.total.iloc[0] / sg.total * 100,
            label=f"{arch}",
            marker=mkr,
        )

        plt.figure("norm-scaling")
        markers = itertools.cycle(marker_shapes)
        sg.sort_values(by=["nranks"], inplace=True)
        if args.breakdown:
            plt.semilogx(
                sg.nranks,
                sg["norm-communication"],
                label=f"{arch} communication",
                marker=mkr,
            )
            plt.semilogx(
                sg.nranks,
                sg["norm-computation"],
                label=f"{arch} computation",
                marker=mkr,
            )
        plt.semilogx(sg.nranks, sg["norm-total"], label=f"{arch}", marker=mkr)

        plt.figure("scaling-entities")
        markers = itertools.cycle(marker_shapes)
        sg.sort_values(by=["entities_per_rank"], inplace=True)
        if args.breakdown:
            plt.loglog(
                sg.entities_per_rank,
                sg.communication,
                label=f"{arch} communication",
                marker=mkr,
            )
            plt.loglog(
                sg.entities_per_rank,
                sg.computation,
                label=f"{arch} computation",
                marker=mkr,
            )
        plt.loglog(
            sg.entities_per_rank,
            sg.total,
            label=f"{arch}",
            marker=mkr,
        )
        plt.loglog(
            sg.entities_per_rank,
            sg.theory_entities,
            label="perfect scaling",
            color="k",
            ls="-",
            zorder=0,
        )

        plt.figure("scaling-entities-gvt")
        markers = itertools.cycle(marker_shapes)
        sg.sort_values(by=["entities_per_rank"], inplace=True)
        if args.breakdown:
            plt.loglog(
                sg.entities_per_rank,
                sg.gvt / sg.communication,
                label=f"{arch} communication",
                marker=mkr,
            )
            plt.loglog(
                sg.entities_per_rank,
                sg.gvt / sg.computation,
                label=f"{arch} computation",
                marker=mkr,
            )
        plt.loglog(
            sg.entities_per_rank,
            sg.gvt / sg.total,
            label=f"{arch}",
            marker=mkr,
        )
        plt.loglog(
            sg.entities_per_rank,
            sg.gvt / sg.theory_entities,
            label="perfect scaling",
            color="k",
            ls="-",
            zorder=0,
        )

        plt.figure("scaling-entities-rate")
        markers = itertools.cycle(marker_shapes)
        sg.sort_values(by=["entities_per_rank"], inplace=True)
        plt.loglog(
            sg.entities_per_rank,
            sg.avg_rate,
            label=f"{arch}",
            marker=mkr,
        )
        # plt.loglog(
        #     sg.entities_per_rank,
        #     sg.gvt / sg.theory_entities,
        #     label="perfect scaling",
        #     color="k",
        #     ls="-",
        #     zorder=0,
        # )

        plt.figure("scaling-rate-efficiency")
        markers = itertools.cycle(marker_shapes)
        sg.sort_values(by=["nranks"], inplace=True)
        plt.semilogx(
            sg.nranks,
            # sg.avg_rate / sg.entities / ( sg.avg_rate / sg.entities).iloc[0],
            (sg.avg_rate / sg.entities) / (sg.avg_rate / sg.entities).iloc[0] * 100,
            label=f"{arch}",
            marker=mkr,
        )

        plt.figure("lp-entity-time")
        markers = itertools.cycle(marker_shapes)
        sg.sort_values(by=["entities_per_lp"], inplace=True)
        plt.semilogx(
            sg["entities_per_lp"],
            sg.gvt / sg.total,
            label=f"{arch}",
            marker=mkr,
        )

        plt.figure("lp-entity-rate")
        markers = itertools.cycle(marker_shapes)
        sg.sort_values(by=["entities_per_lp"], inplace=True)
        plt.semilogx(
            sg["entities_per_lp"],
            sg.avg_rate,
            label=f"{arch}",
            marker=mkr,
        )

        for rlb in range(1, 5):
            plt.figure(f"rollback_{rlb}")
            markers = itertools.cycle(marker_shapes)
            sg.sort_values(by=["entities_per_lp"], inplace=True)
            plt.semilogx(
                sg["entities_per_lp"],
                sg[f"rollback_{rlb}"],
                label=f"{arch}",
                marker=mkr,
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
        legend = legend_without_duplicates(plt.gca())
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("scaling-time-gvt")
        plt.xlabel(f"Number of ranks")
        plt.ylabel(r"$g / t~[-]$")
        legend = legend_without_duplicates(plt.gca())
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("scaling-efficiency")
        plt.axhline(y=100, color="k", ls="-", label="Perfect scaling", zorder=0)
        plt.xlabel(f"Number of ranks")
        plt.ylabel(r"Parallel efficiency$~[\%]$")
        legend = legend_without_duplicates(plt.gca())
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("norm-scaling")
        plt.xlabel(f"Number of ranks")
        plt.ylabel(r"$t~[-]$")
        legend = legend_without_duplicates(plt.gca())
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("scaling-entities")
        plt.xlabel(f"Entities per rank")
        plt.ylabel(r"$t~[s]$")
        legend = legend_without_duplicates(plt.gca())
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("scaling-entities-gvt")
        plt.xlabel(f"Entities per rank")
        plt.ylabel(r"$g / t~[-]$")
        legend = legend_without_duplicates(plt.gca())
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("scaling-entities-rate")
        plt.xlabel(f"Entities per rank")
        plt.ylabel(r"$r~[\#/s]$")
        legend = legend_without_duplicates(plt.gca())
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("scaling-rate-efficiency")
        plt.axhline(y=100, color="k", ls="-", label="Perfect scaling", zorder=0)
        plt.xlabel(f"Number of ranks")
        plt.ylabel(r"Parallel efficiency$~[\%]$")
        legend = legend_without_duplicates(plt.gca())
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("lp-entity-time")
        plt.xlabel(r"$n_e / n_{lp}$")
        plt.ylabel(r"$g / t~[-]$")
        legend = legend_without_duplicates(plt.gca())
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("lp-entity-rate")
        plt.xlabel(r"$n_e / n_{lp}$")
        plt.ylabel(r"$r~[\#/s]$")
        legend = legend_without_duplicates(plt.gca())
        plt.tight_layout()
        pdf.savefig(dpi=300)

        for rlb in range(1, 5):
            plt.figure(f"rollback_{rlb}")
            plt.xlabel(r"$n_e / n_{lp}$")
            plt.ylabel(f"{rlb} rollback $~[\#]$")
            legend = legend_without_duplicates(plt.gca())
            plt.tight_layout()
            pdf.savefig(dpi=300)

        plt.figure("gvt")
        plt.xlabel(r"$s~[-]$")
        plt.ylabel(r"$g~[s]$")
        legend = legend_without_duplicates(plt.gca())
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("rate")
        plt.xlabel(r"$s~[-]$")
        plt.ylabel(r"$r~[\#/s]$")
        legend = legend_without_duplicates(plt.gca())
        # plt.ylim([0, 3e6])
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure(f"rlbk")
        plt.xlabel(r"$s~[-]$")
        plt.ylabel(r"$n~[\#]$")
        legend = legend_without_duplicates(plt.gca())
        plt.tight_layout()
        pdf.savefig(dpi=300)

        for typ in types:
            plt.figure(f"{typ}s")
            plt.xlabel(r"$s~[-]$")
            plt.ylabel(f"$m_{typ[0]}~[\\#]$")
            plt.tight_layout()
            pdf.savefig(dpi=300)
