import os
import collections
import pathlib
import sys

import argh
import pandas
import matplotlib.pyplot as plt


def go(input_dir):
    res = collections.defaultdict(list)
    graph_output_dir = pathlib.Path(input_dir) / "aucgraphs"
    graph_output_dir.mkdir()

    for k, d in enumerate(os.popen(f"ls -trhd {input_dir}/k*").readlines()):
        csvfn = pathlib.Path(d.strip()) / "roc.csv"
        df = pandas.read_csv(csvfn)
        keys = df["col_name"].tolist()
        print(df)
        for i, key in enumerate(keys):
            key_auc = df.at[i, "auc"]
            res[key].append(key_auc)

    for key, Y in res.items():
        plt.figure()
        plt.plot(Y)
        plt.ylim([0.0, 1.05])
        plt.title(f"ROC auc for {key}")
        plt.xlabel("cutoff")
        plt.ylabel("auc")
        plt.savefig(f"{graph_output_dir}/aucgraph_{key}.png")


if __name__ == "__main__":
    argh.dispatch_command(go)
