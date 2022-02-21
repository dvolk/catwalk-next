import collections
import math
import statistics
import pathlib

import argh
import igraph
import matplotlib.pyplot as plt
import pandas
import requests
import scipy.stats
import sklearn.metrics


def get_sample_list():
    return requests.get("http://localhost:5000/list_samples").json()


def get_neighbours(sample_name, snp_distance):
    return requests.get(
        f"http://localhost:5000/neighbours/{sample_name}/{snp_distance}"
    ).json()


def get_pairwise_distances(sample_name_list):
    return requests.post(
        "http://localhost:5000/get_pairwise_distances", json=sample_name_list
    ).json()


def count_elem_lessthan(xs, p):
    return len([x for x in xs if x <= p])


def make_csv_and_graphs(cutoff_distance, outprefix):
    sample_names = get_sample_list()
    sample_data = collections.defaultdict(dict)

    percentiles = [5, 25, 50, 75, 95]
    neighbour_counts_at = range(0, 16)

    pathlib.Path(outprefix).mkdir(exist_ok=True)
    pathlib.Path(f"{outprefix}/igraph_pngs").mkdir(exist_ok=True)

    for sample_index, sample_name in enumerate(sample_names):
        if "+" in sample_name and ":" not in sample_name:
            sample_data[sample_name]["mixed"] = True
        else:
            sample_data[sample_name]["mixed"] = False

        #
        # part 1. analysis of neighbours of sample
        #
        sample_neighbours_data = get_neighbours(sample_name, cutoff_distance)

        sample_data[sample_name]["sample_name"] = sample_name
        # add number of neighbours
        sample_data[sample_name]["1_count_neighbours"] = len(sample_neighbours_data)

        sample_neighbour_distances = [int(x[1]) for x in sample_neighbours_data]
        if len(sample_neighbour_distances) > 1:
            sample_data[sample_name]["1_mean"] = statistics.mean(
                sample_neighbour_distances
            )
            sample_data[sample_name]["1_sd"] = statistics.stdev(
                sample_neighbour_distances
            )
            for percentile in percentiles:
                r = scipy.stats.scoreatpercentile(
                    sample_neighbour_distances, percentile
                )
                sample_data[sample_name][f"1_percentile_{percentile}"] = r

            for neighbour_count_at in neighbour_counts_at:
                sample_data[sample_name][
                    f"1_neighbours_snp_{neighbour_count_at}"
                ] = count_elem_lessthan(sample_neighbour_distances, neighbour_count_at)

        #
        # part 2. analysis of sample neighbour pair matrix
        #
        sample_neighbour_names = [x[0] for x in sample_neighbours_data]
        sample_neighbour_pairwise_distance_data = get_pairwise_distances(
            sample_neighbour_names + [sample_name]
        )
        sample_neighbour_pairwise_distance_data = [
            (x[0], x[1], int(x[2])) for x in sample_neighbour_pairwise_distance_data
        ]

        sample_neighbour_pairwise_edges = [
            (x[0], x[1]) for x in sample_neighbour_pairwise_distance_data
        ]
        sample_neighbour_pairwise_distances = [
            int(x[2]) for x in sample_neighbour_pairwise_distance_data
        ]

        sample_data[sample_name]["2_count_neighbour_matrix"] = len(
            sample_neighbour_pairwise_distances
        )
        if len(sample_neighbour_pairwise_distances) >= 2:
            sample_data[sample_name]["2_mean"] = statistics.mean(
                sample_neighbour_pairwise_distances
            )
            sample_data[sample_name]["2_sd"] = statistics.stdev(
                sample_neighbour_pairwise_distances
            )
            for percentile in percentiles:
                r = scipy.stats.scoreatpercentile(
                    sample_neighbour_pairwise_distances, percentile
                )
                sample_data[sample_name][f"2_percentile_{percentile}"] = r

            for neighbour_count_at in neighbour_counts_at:
                sample_data[sample_name][
                    f"2_neighbours_snp_{neighbour_count_at}"
                ] = count_elem_lessthan(
                    sample_neighbour_pairwise_distances, neighbour_count_at
                )

            all_pairwise_names = sample_neighbour_names + [sample_name]

            g = igraph.Graph()
            g.add_vertices(len(all_pairwise_names))
            vertex_indices = dict()
            for i in range(len(g.vs)):
                g.vs[i]["id"] = i
                vertex_indices[all_pairwise_names[i]] = i
                g.vs[i]["label"] = all_pairwise_names[i]

            edges_relabeled = list()

            for v1, v2 in sample_neighbour_pairwise_edges:
                edges_relabeled.append((vertex_indices[v1], vertex_indices[v2]))
            g.add_edges(edges_relabeled)
            weights = [x + 0.01 for x in sample_neighbour_pairwise_distances]
            g.es["weight"] = weights
            # g.es["label"] = weights
            g.es["curved"] = False

            # rewire_trials = len(sample_neighbour_pairwise_distance_data)
            # g.rewire(n=rewire_trials)

            # how does 3sqrt scale to 0-1?
            vertex_betweenness = igraph.rescale(
                g.betweenness(weights=g.es["weight"]),
                clamp=True,
                scale=lambda x: math.pow(x, 1 / 3),
            )
            sample_data[sample_name][
                "2_sample_vertex_betweenness"
            ] = vertex_betweenness[-1]

            edge_betweenness = igraph.rescale(
                g.edge_betweenness(weights=g.es["weight"]),
                clamp=True,
                scale=lambda x: math.pow(x, 1 / 2),
            )

            vertex_size = igraph.rescale(vertex_betweenness, (10, 100))
            # print(vertex_size)
            print(sample_index, sample_name, "writing graph")
            igraph.plot(
                g,
                f"{outprefix}/igraph_pngs/{sample_index}.png",
                layout=g.layout_lgl(),
                bbox=(1000, 1000),
                vertex_size=vertex_size,
                vertex_frame_width=0.2,
                edge_color="grey",
                edge_width=igraph.rescale(edge_betweenness, (0.5, 5.0)),
            )

    df = pandas.DataFrame(sample_data.values())
    df.to_csv(f"{outprefix}/mixanalysis.csv")


def make_roc_curve_graphs(outprefix):
    df = pandas.read_csv(f"{outprefix}/mixanalysis-cleaned.csv")
    y_true = list(map(int, df["mixed"].values))
    pathlib.Path(outprefix).mkdir(exist_ok=True)
    pathlib.Path(f"{outprefix}/roc_curves").mkdir(exist_ok=True)
    out = collections.defaultdict(dict)
    for col_name in df.keys()[3:]:
        dat = df[col_name].values
        fpr, tpr, _ = sklearn.metrics.roc_curve(y_true, dat)
        score = sklearn.metrics.auc(fpr, tpr)
        cor = scipy.stats.pearsonr(y_true, dat)
        out[col_name]["col_name"] = col_name
        out[col_name]["auc"] = score
        out[col_name]["pearsonr"] = cor[0]
        print(f"{col_name}, aoc: {score}")
        plt.figure()
        plt.title(f"ROC curve for {col_name}")
        plt.plot(
            fpr,
            tpr,
            label=f"ROC curve for {col_name} (area={score:.2f}, cor={cor[0]:.2f})",
            lw=2,
        )
        plt.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.legend(loc="lower right")
        plt.savefig(f"{outprefix}/roc_curves/roc_{col_name}.png")

    df = pandas.DataFrame(out.values())
    df.to_csv(f"{outprefix}/roc.csv")


if __name__ == "__main__":
    argh.dispatch_commands([make_csv_and_graphs, make_roc_curve_graphs])
