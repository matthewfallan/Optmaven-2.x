import os
from collections import OrderedDict

import numpy as np

import standards


def kmeans(data, k, tolerance=standards.DefaultKmeansTolerance, max_iterations=standards.DefaultKmeansMaxIterations):
    if k > len(data):
        raise ValueError("Cannot set k to {}, greater than number of data ({}).".format(k, len(data)))
    centroids = random_centroids(data, k)  
    movement = None
    iteration = 1
    while iteration <= max_iterations and (movement is None or movement > tolerance):
        clusters = update_clusters(data, centroids)
        movement = update_centroids(clusters, centroids)
        iteration += 1
    #FIXME did I correctly implement this error calculation?
    max_mse = np.max([np.mean((np.array(map(to_array, cluster)) - centroid)**2) for cluster, centroid in zip(clusters, centroids)])
    return clusters, max_mse
   

def update_clusters(data, centroids):
    k = len(centroids)
    clusters = [list() for i in range(k)]
    for item in data:
        nearest_cluster_index = np.argmin([np.linalg.norm(to_array(item) - centroid) for centroid in centroids])
        clusters[nearest_cluster_index].append(item)

    # If any cluster is empty then assign one point
    # from data set randomly so as to not have empty
    # clusters and 0 means.
    for i in range(k):
        cluster_i = clusters[i]
        if len(cluster_i) == 0:
            for j in np.random.permutation(range(k)):
                cluster_j = clusters[j]
                n_j = len(cluster_j)
                if n_j > 1:
                    cluster_i.append(cluster_j.pop(np.random.randint(n_j)))
    return clusters


def update_centroids(clusters, centroids):
    if len(clusters) != len(centroids):
        raise ValueError("Clusters and Centroids have different lengths.")
    square_distance = 0.0
    for index, cluster in enumerate(clusters):
        displacement = np.mean(map(to_array, cluster), axis=0) - centroids[index]
        square_distance += np.dot(displacement, displacement)
        centroids[index] += displacement
    return np.sqrt(square_distance)


def random_centroids(data, k):
    centroids = map(to_array, np.random.choice(data, size=k, replace=False).tolist())
    return centroids


def to_array(iterable):
    if isinstance(iterable, np.ndarray):
        return iterable
    else:
        return np.array(list(iterable))


def optimal_kmeans(data, threshold=standards.DefaultKmeansOptKThreshold, tolerance=standards.DefaultKmeansTolerance, max_iterations=standards.DefaultKmeansMaxIterations):
    print "KMEANS"
    k = 1
    clusters_opt = None
    mse_1 = None
    mse_k_lowest = None
    while k <= len(data) and (mse_k_lowest is None or mse_k_lowest / mse_1 > threshold):
        clusters, mse = kmeans(data, k, tolerance, max_iterations)
        print k, mse
        if k == 1:
            mse_1 = mse
        if mse_k_lowest is None or mse < mse_k_lowest:
            mse_k_lowest = mse
            clusters_opt = clusters
        k += 1
    print "FINISHED", len(clusters_opt)
    return clusters_opt
