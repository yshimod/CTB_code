import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

plt.rcParams["font.family"] = "Helvetica"
plt.rcParams["mathtext.fontset"] = "stixsans"
plt.rcParams["pdf.fonttype"] = 42

sns.set_theme(
    style="whitegrid",
    rc={"font.family": "Helvetica", "mathtext.fontset": "stixsans", "pdf.fonttype": 42},
)


df = pd.read_stata("base_estimated.dta")

df["type"] = (df["labnumber"] // 1000).astype(int)
df["type_sd"] = df["type"] % 10
df["type_delta"] = ((df["type"] // 10) % 10).astype(int)
df["type_rho"] = ((df["type"] // 100) % 10).astype(int)
df["type_beta"] = (df["type"] // 1000).astype(int)

df = df[df["type_rho"] > 2]


gt_delta = np.round(np.linspace(0.9912, 1.0025, 10), 4).tolist()

gt_delta_lower = [
    (-i - 1) * (gt_delta[-1] - gt_delta[0]) / 9 + gt_delta[0] for i in range(10)
]
gt_delta_upper = [
    (i + 10) * (gt_delta[-1] - gt_delta[0]) / 9 + gt_delta[0] for i in range(10)
]

gt_beta = np.round(np.linspace(0.85, 1.12, 10), 2).tolist()

gt_beta_lower = [
    (-i - 1) * (gt_beta[-1] - gt_beta[0]) / 9 + gt_beta[0] for i in range(10)
]
gt_beta_upper = [
    (i + 10) * (gt_beta[-1] - gt_beta[0]) / 9 + gt_beta[0] for i in range(10)
]

gt_lnsigma = np.linspace(-2, 5, 10).tolist()

gt_lnsigma_lower = [(-i - 1) * 7 / 9 - 2 for i in range(5)]
gt_lnsigma_upper = [(i + 10) * 7 / 9 - 2 for i in range(5)]


def figure1():
    plt.clf()

    _, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    tmpplotdata = (
        df.groupby(["type_sd", "type_delta"])["neqone_delta"]
        .mean()
        .values.reshape(5, 10)
    )

    for i in range(5):
        ax1.plot(
            np.arange(0, 7),
            tmpplotdata[i, 0:7],
            marker="o",
            color=plt.get_cmap("crest_r", 5)(i / 5),
            linestyle="-",
            lw=3,
        )
        ax1.plot(
            np.arange(8, 10),
            tmpplotdata[i, 8:10],
            marker="o",
            color=plt.get_cmap("crest_r", 5)(i / 5),
            linestyle="-",
            lw=3,
        )

    tmpplotdata = (
        df.groupby(["type_sd", "type_beta"])["neqone_beta"].mean().values.reshape(5, 10)
    )

    for i in range(5):
        ax2.plot(
            np.arange(0, 5),
            tmpplotdata[i, 0:5],
            marker="o",
            color=plt.get_cmap("crest_r", 5)(i / 5),
            linestyle="-",
            lw=3,
        )
        ax2.plot(
            np.arange(6, 10),
            tmpplotdata[i, 6:10],
            marker="o",
            color=plt.get_cmap("crest_r", 5)(i / 5),
            linestyle="-",
            lw=3,
            label="$s$ = {:.2f}".format([0.01, 0.05, 0.10, 0.15, 0.20][i]),
        )

    ax1.set_xlabel("Ground-truth values of $\delta$", fontsize=16)
    ax2.set_xlabel("Ground-truth values of $\\beta$", fontsize=16, labelpad=14)

    ax1.set_ylabel("The rate of successful rejection", fontsize=16)

    ax1.set_ylim([-0.05, 1.05])
    ax2.set_ylim([-0.05, 1.05])

    ax1.set_xticks(list(range(10)))
    ax1.set_xticklabels(
        ["${:.4f}$".format(v) for v in gt_delta],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=16,
    )

    ax2.set_xticks(list(range(10)))
    ax2.set_xticklabels(
        ["${:.2f}$".format(v) for v in gt_beta],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=16,
    )

    ax1.set_yticks([v * 2 / 10 for v in range(6)])
    ax1.set_yticklabels(["${:.1f}$".format(v * 2 / 10) for v in range(6)], fontsize=16)
    ax2.set_yticks([v * 2 / 10 for v in range(6)])
    ax2.set_yticklabels([])

    leg = ax2.legend(
        bbox_to_anchor=(1.02, 0.5), loc="center left", borderaxespad=0, fontsize=16
    )
    leg.get_frame().set_edgecolor("grey")
    leg.get_frame().set_linewidth(0.4)

    ax1.text(
        0,
        1.1,
        "a",
        verticalalignment="top",
        transform=ax1.transAxes,
        fontname="Arial",
        fontweight="bold",
    ).set_fontsize(20)
    ax2.text(
        0,
        1.1,
        "b",
        verticalalignment="top",
        transform=ax2.transAxes,
        fontname="Arial",
        fontweight="bold",
    ).set_fontsize(20)

    ax1.grid(True)
    ax2.grid(True)

    plt.subplots_adjust(left=0.07, right=0.83, bottom=0.20, top=0.92, wspace=0.05)

    plt.savefig("Fig1.pdf")


def figure2():
    plt.clf()

    _, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    for i, v in enumerate(gt_delta):
        ax1.plot(
            [i - 0.52, i + 0.52],
            [v, v],
            linestyle="-",
            linewidth=1,
            c="red",
            zorder=0.6,
        )

    for i, v in enumerate(gt_beta):
        plt.plot(
            [i - 0.52, i + 0.52],
            [v, v],
            linestyle="-",
            linewidth=1,
            c="red",
            zorder=0.6,
        )

    sns.boxplot(
        data=df.query("type_sd == 1 | type_sd == 3 | type_sd == 5"),
        x="delta",
        y="deltai",
        hue="type_sd",
        showfliers=False,
        whis=[5, 95],
        linewidth=0.8,
        palette=sns.color_palette("tab10")[:3],
        ax=ax1,
    )

    sns.boxplot(
        data=df.query("type_sd == 1 | type_sd == 3 | type_sd == 5"),
        x="beta",
        y="betai",
        hue="type_sd",
        showfliers=False,
        whis=[5, 95],
        linewidth=0.8,
        palette=sns.color_palette("tab10")[:3],
        ax=ax2,
    )

    ytickslist1 = gt_delta_lower[::-1] + gt_delta + gt_delta_upper
    ax1.set_yticks(ytickslist1)
    ax1.set_yticklabels(ytickslist1, fontsize=16)
    ax1.yaxis.set_major_formatter(plt.FormatStrFormatter("$%.4f$"))
    ax1.set_xticks(list(range(10)))
    ax1.set_xticklabels(
        ["${:.4f}$".format(v) for v in gt_delta],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=16,
    )

    ytickslist2 = gt_beta_lower[::-1] + gt_beta + gt_beta_upper
    ax2.set_yticks(ytickslist2)
    ax2.set_yticklabels(
        ["${:.2f}$".format(v) if i % 2 else "" for i, v in enumerate(ytickslist2)],
        fontsize=16,
    )
    ax2.set_xticks(list(range(10)))
    ax2.set_xticklabels(
        ["${:.2f}$".format(v) for v in gt_beta],
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=16,
    )

    ax1.set_xlim([-0.5, 9.5])
    ax1.set_ylim([0.987, 1.0045])

    ax2.set_xlim([-0.5, 9.5])
    ax2.set_ylim([0.69, 1.34])

    ax1.grid(axis="y", linewidth=0.5)

    ax2.grid(axis="y", linewidth=0.5)

    ax1.set_xlabel("Ground-truth values of $\delta$", fontsize=16)
    ax1.set_ylabel("Estimates of $\delta$", fontsize=16)

    ax2.set_xlabel("Ground-truth values of $\\beta$", fontsize=16, labelpad=14)
    ax2.set_ylabel("Estimates of $\\beta$", fontsize=16)

    ax1.get_legend().remove()

    handler, label = ax2.get_legend_handles_labels()
    ax2.get_legend().remove()

    legloc_x = 1.05
    legloc_y = 0.35

    ax2.text(
        legloc_x,
        legloc_y,
        "$s=0.01$",
        fontsize=16,
        va="bottom",
        ha="left",
        transform=ax2.transAxes,
        rotation=90,
    )
    ax2.text(
        legloc_x + 0.06,
        legloc_y,
        "$s=0.10$",
        fontsize=16,
        va="bottom",
        ha="left",
        transform=ax2.transAxes,
        rotation=90,
    )
    ax2.text(
        legloc_x + 0.12,
        legloc_y,
        "$s=0.20$",
        fontsize=16,
        va="bottom",
        ha="left",
        transform=ax2.transAxes,
        rotation=90,
    )

    ax2.add_patch(
        patches.Rectangle(
            xy=(legloc_x + 0.005, legloc_y + 0.2),
            width=0.035,
            height=0.1,
            fill=True,
            ec=handler[0].get_edgecolor(),
            color=handler[0].get_facecolor(),
            transform=ax2.transAxes,
            clip_on=False,
        )
    )
    ax2.add_patch(
        patches.Rectangle(
            xy=(legloc_x + 0.06 + 0.005, legloc_y + 0.2),
            width=0.035,
            height=0.1,
            fill=True,
            ec=handler[1].get_edgecolor(),
            color=handler[1].get_facecolor(),
            transform=ax2.transAxes,
            clip_on=False,
        )
    )
    ax2.add_patch(
        patches.Rectangle(
            xy=(legloc_x + 0.12 + 0.005, legloc_y + 0.2),
            width=0.035,
            height=0.1,
            fill=True,
            ec=handler[2].get_edgecolor(),
            color=handler[2].get_facecolor(),
            transform=ax2.transAxes,
            clip_on=False,
        )
    )

    ax2.add_patch(
        patches.FancyBboxPatch(
            xy=(legloc_x - 0.005, legloc_y - 0.02),
            width=0.185,
            height=0.35,
            boxstyle="round, pad=.01",
            lw=0.4,
            fill=False,
            ec="grey",
            transform=ax2.transAxes,
            clip_on=False,
        )
    )

    ax1.text(
        -0.1,
        1.1,
        "a",
        verticalalignment="top",
        transform=ax1.transAxes,
        fontname="Arial",
        fontweight="bold",
    ).set_fontsize(20)
    ax2.text(
        -0.1,
        1.1,
        "b",
        verticalalignment="top",
        transform=ax2.transAxes,
        fontname="Arial",
        fontweight="bold",
    ).set_fontsize(20)

    plt.subplots_adjust(left=0.095, right=0.91, bottom=0.17, top=0.92, wspace=0.25)

    plt.savefig("Fig2.pdf")


def figure3():
    plt.clf()

    plt.figure(figsize=(6, 5))

    for i, v in enumerate(gt_lnsigma[3:]):
        plt.plot(
            [i - 0.52, i + 0.52],
            [v, v],
            linestyle="-",
            linewidth=1,
            c="red",
            zorder=0.6,
        )

    plt.plot(
        [-0.5, 6.5], [5.5, 5.5], linestyle="dashed", linewidth=1, c="k", zorder=0.6
    )
    plt.plot(
        [-0.5, 6.5], [-2.5, -2.5], linestyle="dashed", linewidth=1, c="k", zorder=0.6
    )

    sns.boxplot(
        data=df.query("type_sd == 1 | type_sd == 3 | type_sd == 5"),
        x="lnsigma",
        y="lnsigmai",
        hue="type_sd",
        showfliers=False,
        whis=[5, 95],
        linewidth=0.8,
        palette=sns.color_palette("tab10")[:3],
    )

    plt.gca().set_yticks(gt_lnsigma_lower[::-1] + gt_lnsigma + gt_lnsigma_upper)
    plt.gca().set_yticklabels(
        gt_lnsigma_lower + gt_lnsigma + gt_lnsigma_upper, fontsize=16
    )
    plt.gca().yaxis.set_major_formatter(plt.FormatStrFormatter("$%.2f$"))

    plt.gca().set_xticks(list(range(7)))
    plt.gca().set_xticklabels(
        ["${:.2f}$".format(v) for v in gt_lnsigma[3:]], fontsize=16
    )

    plt.xlim(-0.5, 6.5)
    plt.ylim(-2.7, 5.7)

    plt.grid(axis="y", linewidth=0.5)

    plt.xlabel("Ground-truth values of $\ln\sigma$", fontsize=16)
    plt.ylabel("Estimates of $\ln\sigma$", fontsize=16)

    handler, _ = plt.gca().get_legend_handles_labels()
    plt.gca().get_legend().remove()

    legloc_x = 1.05
    legloc_y = 0.35

    plt.gca().text(
        legloc_x,
        legloc_y,
        "$s=0.01$",
        fontsize=16,
        va="bottom",
        ha="left",
        transform=plt.gca().transAxes,
        rotation=90,
    )
    plt.gca().text(
        legloc_x + 0.06,
        legloc_y,
        "$s=0.10$",
        fontsize=16,
        va="bottom",
        ha="left",
        transform=plt.gca().transAxes,
        rotation=90,
    )
    plt.gca().text(
        legloc_x + 0.12,
        legloc_y,
        "$s=0.20$",
        fontsize=16,
        va="bottom",
        ha="left",
        transform=plt.gca().transAxes,
        rotation=90,
    )

    plt.gca().add_patch(
        patches.Rectangle(
            xy=(legloc_x + 0.005, legloc_y + 0.2),
            width=0.035,
            height=0.1,
            fill=True,
            ec=handler[0].get_edgecolor(),
            color=handler[0].get_facecolor(),
            transform=plt.gca().transAxes,
            clip_on=False,
            zorder=0.9,
        )
    )
    plt.gca().add_patch(
        patches.Rectangle(
            xy=(legloc_x + 0.06 + 0.005, legloc_y + 0.2),
            width=0.035,
            height=0.1,
            fill=True,
            ec=handler[1].get_edgecolor(),
            color=handler[1].get_facecolor(),
            transform=plt.gca().transAxes,
            clip_on=False,
            zorder=0.9,
        )
    )
    plt.gca().add_patch(
        patches.Rectangle(
            xy=(legloc_x + 0.12 + 0.005, legloc_y + 0.2),
            width=0.035,
            height=0.1,
            fill=True,
            ec=handler[2].get_edgecolor(),
            color=handler[2].get_facecolor(),
            transform=plt.gca().transAxes,
            clip_on=False,
            zorder=0.9,
        )
    )

    plt.gca().add_patch(
        patches.FancyBboxPatch(
            xy=(legloc_x - 0.005, legloc_y - 0.02),
            width=0.185,
            height=0.35,
            boxstyle="round, pad=.01",
            lw=0.4,
            fill=True,
            fc="white",
            ec="grey",
            transform=plt.gca().transAxes,
            clip_on=False,
            zorder=0.85,
        )
    )

    plt.subplots_adjust(left=0.18, right=0.83, bottom=0.12, top=0.99)

    plt.savefig("Fig3.pdf")


##### main #####
figure1()
figure2()
figure3()
