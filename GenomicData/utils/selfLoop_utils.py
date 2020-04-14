import numpy as np

def remove_selfLoop(edge_index):
    assert edge_index.shape[0] == 2
    i, j = edge_index
    duplicated_index = (i == j)
    return edge_index[:, ~duplicated_index]


def add_selfLoop(edge_index):
    assert edge_index.shape[0] == 2
    node_set = np.unique(edge_index)
    edge_index = np.concatenate([edge_index, np.vstack([node_set]*2)], axis=1)
    return edge_index