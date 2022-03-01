"""
Draw a heatmap from the analysis results directories. The heatmap values are
the AUC values for specific SNP cut-off and measured variable.

Command line arguments are the directories to be used for processing.

adapted from matplotlib heatmap example.
"""

import collections
import pathlib

import argh
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from natsort import natsorted


def heatmap(data, row_labels, col_labels, ax=None, cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right", rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im,
    data=None,
    valfmt="{x:.2f}",
    textcolors=("black", "white"),
    threshold=None,
    **textkw
):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.0

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


@argh.arg("-d", "--dirs", nargs="+", type=str)
def main(dirs=[]):
    # parse data from dirs/k*/roc.csv
    rows = list()
    cols = list()
    measure = collections.defaultdict(list)
    for d in dirs:
        d = pathlib.Path(d)
        k_dirs = list(natsorted(d.glob("k*"), key=str))
        print(cols)
        print(k_dirs)

        for k_dir in k_dirs:
            cols.append(k_dir)
            roc_data = [
                line.split(",") for line in open(k_dir / "roc.csv").readlines()[1:]
            ]
            print(roc_data)
            if not rows:
                for roc_data_line in roc_data:
                    rows.append(roc_data_line[1])
            for roc_data_line in roc_data:
                measure[roc_data_line[1]].append(float(roc_data_line[2]))

    measure = np.array(list(measure.values()))

    print(cols)
    print(len(rows))
    print(len(cols))
    print(len(measure))
    print(len(measure[0]))

    # plot heatmap
    fig, ax = plt.subplots()

    im, cbar = heatmap(
        measure, rows, cols, vmin=0.9, ax=ax, cmap="Reds", cbarlabel="ROC curve AUC"
    )
    annotate_heatmap(im, valfmt="{x:.3f}")

    fig.tight_layout()
    fig.set_size_inches(50, 50)
    plt.savefig("auc_heatmap.png", dpi=200)


if __name__ == "__main__":
    argh.dispatch_command(main)
