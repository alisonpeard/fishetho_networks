import scipy
import networkx as nx
import numpy as np
from sklearn.cluster import KMeans


def sort_dict(dic):
    """Helper function to sort dictionary."""
    return {key: value for key, value in sorted(dic.items())}


def sponge_partition(G, k=2):
    """
    Cluster a signed network using SPONGE method. Isolated nodes are assigned there own clusters.
    
    Parameters:
    ----------
    G : networkx.Graph
    k : int, default=2
        The number of clusters to search for.

    References:
    ----------
    ..[1] Mihai Cucuringu, Peter Davies, Aldo Glielmo, and Hemant Tyagi. Sponge:
        A generalized eigenproblem for clustering signed networks. 2019.
        
    """
    
    G = G.copy()
    
    # lists of nodes
    nodes = [*G.nodes]
    isolated = [*nx.isolates(G)]
    nodes = [node for node in nodes if node not in isolated]
    
    G.remove_nodes_from(isolated)
    
    val, vect  = np.linalg.eigh(sponge(G))
    positions = vect[:, 0: k - 1]
    partition = KMeans(n_clusters=k).fit(positions)
    partition_dict = {node: cluster for node, cluster in zip(nodes, partition.labels_)}
    
    # add back isolated nodes
    for i, node in enumerate(isolated):
        partition_dict[node] = - (i + 1)
        
    partition_dict = sort_dict(partition_dict)
    return partition_dict

    
def sponge(G):
    """
    Calculate the SPONGE matrix of an undirected, weighted, signed graph.

    The SPONGE matrix is a numpy matrix with elements representing the effective
    conductance between pairs of nodes in the input graph. Original code from
    shazia.babul@new.ox.ac.uk.

    Parameters
    ----------
    G : networkx.Graph
        An undirected, weighted, signed graph.

    Returns
    -------
    numpy.matrix
        A matrix of effective conductances between pairs of nodes in the input graph.
    """
    
    selected_edges = [(u,v) for u,v in G.edges if G[u][v]['weight'] < 0]

    G_neg = nx.Graph()
    G_neg.add_nodes_from(G)
    G_neg.add_edges_from(selected_edges)
    degree_neg = [val for (node, val) in G_neg.degree()]
    D_neg = np.diag(degree_neg)
    A_neg = nx.to_numpy_array(G_neg)
    L_neg = D_neg - A_neg

    
    selected_edges = [(u,v) for u,v in G.edges if G[u][v]['weight'] > 0]
    G_pos = nx.Graph()
    G_pos.add_nodes_from(G)
    G_pos.add_edges_from(selected_edges)
    degree_pos = [val for (node, val) in G_pos.degree()]
    D_pos = np.diag(degree_pos)
    A_pos = nx.to_numpy_array(G_pos)
    L_pos = D_pos - A_pos
    
    
    N = L_neg + D_pos
    D, V = scipy.linalg.eigh(N)
    Bs = (V * np.sqrt(D)) @ V.T
    M1 = np.linalg.inv(Bs)
    M2 = L_pos + D_neg
    
    sponge = np.matmul(M1, M2)
    sponge = np.matmul(sponge, M1)
    
    return sponge