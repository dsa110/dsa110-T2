import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

try:
    import seaborn as sns
except ImportError:
    print("no seaborn found. continuing...")
import dsautils.dsa_syslog as dsl

logger = dsl.DsaSyslogger()
logger.subsystem("software")
logger.app("T2")


def plot_clustered(clusterer, clsnr, snrs, data, tab, cols, plot_dir="./"):
    """
    Plot cluster probabilities in color and greyness.

    Parameters
    ----------
    clusterer : from cluster_data().
    clsnr: from get_peak().
    tab, data, snrs: from parse_candsfile().
    cols : parameters to plot.
    plot_dir : where to save output plots.
    """
    # grey for unclustered noise, pure color for well clustered points.
    palette = sns.color_palette()
    cluster_colors = [
        sns.desaturate(palette[col], sat) if col >= 0 else (0.5, 0.5, 0.5)
        for col, sat in zip(clusterer.labels_, clusterer.probabilities_)
    ]

    for i in range(len(cols)):
        fig, ax = plt.subplots()
        ax.cla()
        ax.scatter(data[:, i], snrs, s=3, c=cluster_colors)
        ax.set_xlabel(cols[i])
        ax.set_ylabel("snr")
        ax.set_title("cluster cols:" + str(cols))
        fig.savefig(plot_dir + "snr_" + str(cols[i]) + ".pdf")
        fig.clf()
        plt.close("all")

        for j in range(len(cols)):
            if j > i:
                fig, ax = plt.subplots()
                ax.cla()
                ax.scatter(data[:, i], data[:, j], c=cluster_colors)
                ax.set_xlabel(cols[i])
                ax.set_ylabel(cols[j])

                for k in range(len(clsnr)):
                    imaxsnr = clsnr[k][0]
                    maxsnr = int(clsnr[k][1])
                    ax.scatter(
                        data[:, i][imaxsnr],
                        data[:, j][imaxsnr],
                        s=maxsnr,
                        c="k",
                        marker="*",
                    )
                    ax.text(
                        data[:, i][imaxsnr], data[:, j][imaxsnr], str(maxsnr)
                    )

                ax.set_title("cluster cols:" + str(cols))
                fig.savefig(
                    plot_dir
                    + "cluster_prob_"
                    + cols[i]
                    + "_"
                    + cols[j]
                    + ".pdf"
                )
                fig.clf()
                plt.close("all")


# Below are functions to plot giants.
def plot_dm_hist(
    tab, nbins=30, plot_dir="./", multibeam=False, data_name=None
):
    """
    plot the giants DM histogram

    Parameters
    ----------
    tab : ascii table
        MBHeimdall giants output full table from parse_candsfile or parse_socket.
    nbins : int
        Bin number. The default is 30.
    plot_dir : optional.  Where to save plots. The default is "./".
    multibeam : bool, optional. The default is False.
        If True, plot a histogram for each beam.
    data_name : bool, optional. The default is None.

    Returns
    -------
    Save the plot in plot_dir.
    """
    fig, ax = plt.subplots()
    ax.cla()  # Clear the current axes
    dm_min = 1.0
    # tab['ibeam'] = tab['ibeam'].astype(int)
    if multibeam:
        for beam in range(
            tab["ibeam"].max() + 1
        ):  # due to zero based indexing
            if beam in tab["ibeam"]:
                cands = tab[tab["ibeam"] == beam]
                logbins = np.logspace(
                    np.log10(dm_min), np.log10(cands["dm"].max()), nbins + 1
                )
                vals, edges = np.histogram(cands["dm"], bins=logbins)
                ax.step(
                    edges, np.append(vals, 0.0), where="post", label=str(beam)
                )
                ax.legend(loc=9, ncol=4, fontsize=8)
    else:
        logbins = np.logspace(
            np.log10(dm_min), np.log10(tab["dm"].max()), nbins + 1
        )
        ax.hist(tab["dm"], bins=logbins, histtype="step", label=data_name)

    ### Set global plot window parameters
    ax.set_xlim(dm_min, tab["dm"].max())
    ax.set_xscale("log")
    ax.set_yscale("log")

    # Set ticks and grid
    ax.grid(axis="x", which="both", linewidth=0.2)
    ax.set_axisbelow(True)
    formatter = plt.FormatStrFormatter("%i")
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)

    ax.set_xlabel(r"$\rm DM\;(pc\;cm^{-3})$", size=12)
    ax.set_ylabel("$\\rm Giants count$", size=12)
    ax.set_title("giants dm")
    fig.savefig(
        plot_dir + "giants_dm_hist_multibeam_" + str(multibeam) + ".pdf"
    )
    fig.clf()
    plt.close("all")


def plot_dm_snr(ax, ax_cbar, tab, tsamp=1048e-6):
    """


    Parameters
    ----------
    ax_cbar : TYPE
        DESCRIPTION.
    tab : ascii table
        MBHeimdall giants output full table from parse_candsfile or parse_socket.
    tsamp : optional. The default is 1048e-6.

    Returns
    -------
    None.

    """
    # fig, ax = plt.subplots()
    ax.cla()
    ax_cbar.cla()

    max_snr = 100.0

    ax.set_xlim(6.0, max_snr)
    ax.set_xscale("log")
    ax.set_xlabel("$\\rm SNR$", size=12)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.tick_params(axis="y", labelleft=False)

    # Set ticks and grid
    ax.grid(axis="y", which="both", linewidth=0.2)
    ax.set_axisbelow(True)
    xformat = plt.FormatStrFormatter("%i")
    ax.xaxis.set_major_formatter(xformat)

    xticks = [6.0, 8.0, 10.0, 20.0, 40.0, 100.0]
    xticklabel = ["%i" % tick for tick in xticks]
    ax.xaxis.set_ticks(xticks)
    ax.xaxis.set_ticklabels(xticklabel)
    ax.xaxis.set_minor_locator(plt.NullLocator())

    if tab["snr"].max() > max_snr:
        logger.warning(
            "The SNR of a candidate is higher than the maximum SNRs plotted!"
        )

    colormap = ax.scatter(
        tab["snr"],
        tab["dm"],
        s=tab["ibox"],
        c=tab["ibox"],
        cmap="viridis",
        vmin=0.0,
        vmax=12,
        marker="o",
        picker=2,
        alpha=0.5,
    )
    # Add colorbar
    try:
        colormap
    except NameError:
        ax_cbar.axis("off")
        pass

    ax_cbar.axis("on")
    cbar = plt.colorbar(
        colormap,
        cax=ax_cbar,
        use_gridspec=True,
        label=r"$\rm Boxcar\ width\;(index)$",
    )
    #fix for matplotlib warning: UserWarning: FixedFormatter should only be used together with FixedLocator
    ibox_vals = np.arange(0, 13)
    cbar.set_ticks(ibox_vals)
    widths_ms = (2.0 ** ibox_vals) * tsamp * 1000.0
    ticklabels = [f"{val:.3f}" for val in widths_ms]
    cbar.set_ticklabels(ticklabels)
    
    cbar.set_alpha(1)
    cbar.draw_all()


def plot_time_dm(
    ax,
    tab,
    duration=None,
    snr_cut=8,
    snr_thr=25.0,
    mps=30.0,
    multibeam=False,
    axrange=True,
    axlabel=True,
):
    """Plot time versus DM."""

    ax.cla()
    max_point_size = mps  # diameter in points

    if axrange:
        if duration != None:
            ax.set_xlim(0.0, duration)  # duration in seconds
        else:
            ax.set_xlim(0.0, tab["mjds"].max() * 1.05)

        ax.set_ylim(1.0, tab["dm"].max() * 1.05)
        ax.set_yscale("log")

    if axlabel == True:
        ax.set_xlabel(r"$\rm Time\;(sec)$", size=12)
        ax.set_ylabel(r"$\rm DM\;(pc\;cm^{-3})$", size=12)
    else:
        ax.tick_params(axis="x", labelbottom=False)

    ax.grid(axis="y", which="both", linewidth=0.2)
    ax.set_axisbelow(True)
    yformat = plt.FormatStrFormatter("%i")
    ax.yaxis.set_major_formatter(yformat)

    ax.scatter(
        tab["mjds"],
        tab["dm"],
        s=(max_point_size * tab["snr"] / tab["snr"].mean()) ** 2,
        c=tab["ibox"],
        cmap="viridis",
        vmin=0.0,
        vmax=12,
        marker="o",
        picker=5,
        alpha=0.5,
    )


def plot_beam_time(tab, plot_dir="./"):
    fig, ax = plt.subplots()
    ax.cla()
    colormap = ax.scatter(
        tab["mjds"],
        tab["ibeam"],
        s=(tab["snr"] / tab["snr"].min()) * 5,
        c=tab["ibox"],
        marker="o",
        alpha=0.5,
    )
    fig.colorbar(
        colormap,
        ax=ax,
        use_gridspec=True,
        label=r"$\rm Boxcar\ width\;(index)$",
    )
    ax.set_xlabel("$\\rm mjd (s)$", size=12)
    ax.set_ylabel("$\\rm beam number$", size=12)
    ax.set_title("giants snr for each beam")
    fig.savefig(plot_dir + "giants_beam_time.pdf")


def plot_giants(tab, plot_dir="./"):
    """
    Plot the un-clustere heimdall output giants.out
    Parameters
    ----------
    tab : full table from parse_candsfile.
    plot_dir : optional. The default is "./".
    """

    plot_dm_hist(
        tab, nbins=30, plot_dir=plot_dir, multibeam=False, data_name=None
    )
    plot_dm_hist(
        tab, nbins=30, plot_dir=plot_dir, multibeam=True, data_name=None
    )
    plot_beam_time(tab, plot_dir=plot_dir)

    # subplot or just save plots in each function?
    fig2, ax2 = plt.subplots(1, 2, sharey=True)
    plot_time_dm(
        ax2[0],
        tab,
        duration=None,
        snr_cut=6.5,
        snr_thr=25.0,
        mps=3.0,
        multibeam=False,
        axrange=True,
        axlabel=True,
    )
    divider = make_axes_locatable(ax2[1])
    cbar_ax = divider.append_axes("right", size="5%", pad=0.1)
    plot_dm_snr(ax2[1], cbar_ax, tab)
    fig2.savefig(plot_dir + "giants_dm_time_snr.pdf")
    # fig2.clf()
    plt.close("all")
