import collections
import statistics

import argh
import requests
import scipy.stats

import igraph
from rich import print


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

    for i, sample_name in enumerate(sample_names[0:11]):
        #
        # part 1. analysis of neighbours of sample
        #
        sample_neighbours_data = get_neighbours(sample_name, cutoff_distance)

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
            (x[0], x[1], int(x[2]) + 0.1)
            for x in sample_neighbour_pairwise_distance_data
        ]

        sample_neighbour_pairwise_edges = [
            (x[0], x[1]) for x in sample_neighbour_pairwise_distance_data
        ]
        sample_neighbour_pairwise_distances = [
            x[2] for x in sample_neighbour_pairwise_distance_data
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

            # sample_data[sample_name][
            # "2_debug_distance"
            # ] = sample_neighbour_pairwise_distance_data

            g = igraph.Graph()
            g.add_vertices(sample_neighbour_names + [sample_name])
            g.add_edges(sample_neighbour_pairwise_edges)
            g.es["weight"] = sample_neighbour_distances

            g = igraph.Graph().TupleList(
                sample_neighbour_pairwise_distance_data, weights=True
            )

            rewire_trials = 10 * len(sample_neighbour_pairwise_distance_data)
            g.rewire(n=rewire_trials)

            sample_data[sample_name]["2_centrality"] = g.eigenvector_centrality(
                weights=sample_neighbour_pairwise_distances
            )
            g.layout_auto()
            igraph.plot(g, f"{i}.png")

    print(sample_data)


if __name__ == "__main__":
    argh.dispatch_command(do_analysis)
