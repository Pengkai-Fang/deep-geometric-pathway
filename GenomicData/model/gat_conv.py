import torch
from torch.nn import Parameter
import torch.nn.functional as F
from .message_passing import MessagePassing
from torch_geometric.utils import remove_self_loops, add_self_loops
from .softmax import softmax

from .inits import glorot, zeros


class GATConv(MessagePassing):
    def __init__(self, in_channels, out_channels, heads=1, concat=True,
                 negative_slope=0.2, dropout=0, bias=True, **kwargs):
        super(GATConv, self).__init__(aggr='mean', **kwargs)

        self.num_obs = None
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.heads = heads
        self.concat = concat
        self.negative_slope = negative_slope
        self.dropout = dropout

        self.weight = Parameter(torch.Tensor(in_channels,
                                             heads * out_channels))
        self.att = Parameter(torch.Tensor(1, heads, 2 * out_channels))

        if bias and concat:
            self.bias = Parameter(torch.Tensor(heads * out_channels))
        elif bias and not concat:
            self.bias = Parameter(torch.Tensor(out_channels))
        else:
            self.register_parameter('bias', None)

        self.reset_parameters()

    def reset_parameters(self):
        glorot(self.weight)
        glorot(self.att)
        zeros(self.bias)

    def forward(self, x, edge_index, size=None):
        """"""
        if size is None and torch.is_tensor(x):
            edge_index, _ = remove_self_loops(edge_index)
            edge_index, _ = add_self_loops(edge_index,
                                           num_nodes=x.size(self.node_dim))

        if torch.is_tensor(x):
            x = torch.matmul(x, self.weight)
        else:
            x = (None if x[0] is None else torch.matmul(x[0], self.weight),
                 None if x[1] is None else torch.matmul(x[1], self.weight))

        return self.propagate(edge_index, size=size, x=x)

    def message(self, edge_index_i, x_i, x_j, size_i):
        # Compute attention coefficients.
        if (len(x_j.shape) == 3):
            self.num_obs = x_j.shape[0]
            x_j = x_j.view(self.num_obs, -1, self.heads, self.out_channels)
        else:
            x_j = x_j.view(-1, self.heads, self.out_channels)
        if x_i is None:
            alpha = (x_j * self.att[:, :, self.out_channels:]).sum(dim=-1)
        else:
            if (len(x_i.shape) == 3):
                x_i = x_i.view(self.num_obs, -1, self.heads, self.out_channels)
            else:
                x_i = x_i.view(-1, self.heads, self.out_channels)
            alpha = (torch.cat([x_i, x_j], dim=-1) * self.att).sum(dim=-1)

        alpha = F.leaky_relu(alpha, self.negative_slope)
        alpha = softmax(alpha, edge_index_i, size_i)

        # Sample attention coefficients stochastically.
        alpha = F.dropout(alpha, p=self.dropout, training=self.training)

        if (len(x_j.shape) == 4):
            alpha = alpha.view(self.num_obs, -1, self.heads, 1)
        else:
            alpha = alpha.view(-1, self.heads, 1)
        return x_j * alpha

    def update(self, aggr_out):
        if self.concat is True:
            if self.num_obs:
                aggr_out = aggr_out.view(self.num_obs, -1, self.heads * self.out_channels)
            else:
                aggr_out = aggr_out.view(-1, self.heads * self.out_channels)
        else:
            aggr_out = aggr_out.mean(dim=1)

        if self.bias is not None:
            aggr_out = aggr_out + self.bias
        return aggr_out

    def __repr__(self):
        return '{}({}, {}, heads={})'.format(self.__class__.__name__,
                                             self.in_channels,
                                             self.out_channels, self.heads)
