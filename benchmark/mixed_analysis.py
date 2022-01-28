import collections
import statistics
import json
import math
import random

import argh
import requests
import scipy.stats

import igraph
import pandas


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


def count_neighbours(xs, p):
    count = 0
    for x in xs:
        if x == p:
            count += 1
    return count


def do_analysis(cutoff_distance):
    sample_names = get_sample_list()
    sample_data = collections.defaultdict(dict)

    percentiles = [5, 25, 50, 75, 95]
    neighbour_counts_at = range(0, 16)

    for sample_index, sample_name in enumerate(sample_names):
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
                ] = count_neighbours(sample_neighbour_distances, neighbour_count_at)

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
        if len(sample_neighbour_pairwise_distances) > 3:
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
                ] = count_neighbours(
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
            sample_data[sample_name]["2_sample_edge_betweenness"] = edge_betweenness[-1]

            vertex_size = igraph.rescale(vertex_betweenness, (10, 100))
            # print(vertex_size)
            print(sample_index, sample_name, "writing graph")
            igraph.plot(
                g,
                f"igraph_pngs/{sample_index}.png",
                layout=g.layout_lgl(),
                bbox=(1000, 1000),
                vertex_size=vertex_size,
                vertex_frame_width=0.2,
                edge_color="grey",
                edge_width=igraph.rescale(edge_betweenness, (0.5, 5.0)),
            )

    xs = [[k] + list(v.values()) for k, v in sample_data.items()]

    df = pandas.DataFrame(sample_data.values())
    df.to_csv("mixanalysis.csv")


if __name__ == "__main__":
    argh.dispatch_command(do_analysis)
