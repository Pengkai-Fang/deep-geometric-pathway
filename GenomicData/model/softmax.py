from torch_scatter import scatter_max, scatter_add

def maybe_num_nodes(index, num_nodes=None):
    return index.max().item() + 1 if num_nodes is None else num_nodes

def softmax(src, index, num_nodes=None):
    r"""Computes a sparsely evaluated softmax.
    Given a value tensor :attr:`src`, this function first groups the values
    along the first dimension based on the indices specified in :attr:`index`,
    and then proceeds to compute the softmax individually for each group.

    Args:
        src (Tensor): The source tensor.
        index (LongTensor): The indices of elements for applying the softmax.
        num_nodes (int, optional): The number of nodes, *i.e.*
            :obj:`max_val + 1` of :attr:`index`. (default: :obj:`None`)

    :rtype: :class:`Tensor`
    """

    num_nodes = maybe_num_nodes(index, num_nodes)

    dim = len(src.shape) - 2
    out = src - scatter_max(src, index, dim=dim, dim_size=num_nodes)[0].index_select(dim, index)
    out = out.exp()
    out = out / (
        scatter_add(out, index, dim=dim, dim_size=num_nodes).index_select(dim, index) + 1e-16)

    return out

