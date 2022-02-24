"""
Betweenness benchmark.

Requires a catwalk instance loaded with samples.

Prints csv to stdout.

# example results
#        neighbours_count  vertex_betweenness       time_s
# count       1050.000000         1050.000000  1050.000000
# mean          23.371429            0.238162     0.007021
# std           36.197598            0.286489     0.022450
# min            0.000000            0.000000     0.002069
# 25%            6.000000            0.000000     0.002233
# 50%           12.000000            0.000000     0.002419
# 75%           25.000000            0.500000     0.003010
# max          265.000000            1.000000     0.304404
"""

import argh
import igraph
import requests
import time
import math


def do_one_sample(sample_name, snp_distance):
    def get_neighbours(sample_name, snp_distance):
        return requests.get(
            f"http://localhost:5000/neighbours/{sample_name}/{snp_distance}"
        ).json()

    def get_pairwise_distances(sample_name_list):
        return requests.post(
            "http://localhost:5000/get_pairwise_distances", json=sample_name_list
        ).json()

    neighbours_data = get_neighbours(sample_name, snp_distance)
    all_pairwise_names = [x[0] for x in neighbours_data] + [sample_name]
    sample_neighbour_pairwise_distance_data = get_pairwise_distances(all_pairwise_names)
    sample_neighbour_pairwise_distances = [
        int(x[2]) for x in sample_neighbour_pairwise_distance_data
    ]
    sample_neighbour_pairwise_edges = [
        (x[0], x[1]) for x in sample_neighbour_pairwise_distance_data
    ]

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
    g.es["curved"] = False

    vertex_betweenness = igraph.rescale(
        g.betweenness(weights=g.es["weight"]),
        clamp=True,
        scale=lambda x: math.pow(x, 1 / 3),
    )
    sample_vertex_betweenness = vertex_betweenness[-1]

    return [sample_name, str(len(neighbours_data)), str(sample_vertex_betweenness)]


def main():
    def get_sample_list():
        return requests.get("http://localhost:5000/list_samples").json()

    samples = get_sample_list()

    rows = list()
    for sample in samples:
        start = time.time()
        row = do_one_sample(sample, 3)
        end = time.time()
        rows.append(row + [str(end - start)])

    print("sample_name,neighbours_count,vertex_betweenness,time_s")
    for row in rows:
        print(",".join(row))


if __name__ == "__main__":
    argh.dispatch_command(main)
