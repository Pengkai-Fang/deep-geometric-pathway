import pandas as pd
import torch
import torch.nn as nn
import numpy as np
import scipy.io as sio
import torch.nn.functional as F
import torch.optim as optim
import math
from torch.autograd import Variable
from torch.nn import Parameter
from torch_scatter import scatter_add
from torch_geometric.nn.conv import MessagePassing
from torch_geometric.utils import add_remaining_self_loops


def reset(nn):
    def _reset(item):
        if hasattr(item, 'reset_parameters'):
            item.reset_parameters()

    if nn is not None:
        if hasattr(nn, 'children') and len(list(nn.children())) > 0:
            for item in nn.children():
                _reset(item)
        else:
            _reset(nn)
            
def glorot(tensor):
    if tensor is not None:
        stdv = math.sqrt(6.0 / (tensor.size(-2) + tensor.size(-1)))
        tensor.data.uniform_(-stdv, stdv)

def zeros(tensor):
    if tensor is not None:
        tensor.data.fill_(0)
        
        
dataroot = './rna-genome-data/'
data_x = pd.read_csv(dataroot + 'availPPidata.csv')
data_edge_index = pd.read_csv(dataroot + 'availedge.csv')

sourceroot = './NewCleanDataSetBRCA/'
rna_data = sio.loadmat(sourceroot + 'Data_RNASeq2.mat')
cancerfreedata = rna_data['geneRNASeq2RawMatrix0']
cancerobtdata = rna_data['geneRNASeq2RawMatrix1']
target = np.concatenate([np.zeros((cancerfreedata.shape[1])),np.ones((cancerobtdata.shape[1]))]).astype(int)
data_x = data_x.values
data = np.concatenate([data_x, target.reshape(1,-1)],axis=0)
data = data.T


class GCNConv(MessagePassing):
    def __init__(self, in_channels, out_channels, improved=False, cached=False,
                 bias=True, normalize=True, **kwargs):
        super(GCNConv, self).__init__(aggr='add', **kwargs)

        self.in_channels = in_channels
        self.out_channels = out_channels
        self.improved = improved
        self.cached = cached
        self.normalize = normalize
        self.node_dim = 1

        self.weight = Parameter(torch.Tensor(in_channels, out_channels))

        if bias:
            self.bias = Parameter(torch.Tensor(out_channels))
        else:
            self.register_parameter('bias', None)

        self.reset_parameters()
        
    def reset_parameters(self):
        glorot(self.weight)
        zeros(self.bias)
        self.cached_result = None
        self.cached_num_edges = None

    @staticmethod
    def norm(edge_index, num_nodes, edge_weight=None, improved=False,
             dtype=None):
        if edge_weight is None:
            edge_weight = torch.ones((edge_index.size(1), ), dtype=dtype,
                                     device=edge_index.device)

        fill_value = 1 if not improved else 2
        edge_index, edge_weight = add_remaining_self_loops(
            edge_index, edge_weight, fill_value, num_nodes)

        row, col = edge_index
        deg = scatter_add(edge_weight, row, dim=0, dim_size=num_nodes)
        deg_inv_sqrt = deg.pow(-0.5)
        deg_inv_sqrt[deg_inv_sqrt == float('inf')] = 0

        return edge_index, deg_inv_sqrt[row] * edge_weight * deg_inv_sqrt[col]

    def forward(self, x, edge_index, edge_weight=None):
        """"""
        x = torch.matmul(x, self.weight)

        if self.cached and self.cached_result is not None:
            if edge_index.size(1) != self.cached_num_edges:
                raise RuntimeError(
                    'Cached {} number of edges, but found {}. Please '
                    'disable the caching behavior of this layer by removing '
                    'the `cached=True` argument in its constructor.'.format(
                        self.cached_num_edges, edge_index.size(1)))

        if not self.cached or self.cached_result is None:
            self.cached_num_edges = edge_index.size(1)
            if self.normalize:
                edge_index, norm = self.norm(edge_index, x.size(self.node_dim),
                                             edge_weight, self.improved,
                                             x.dtype)
            else:
                norm = edge_weight
                
            self.cached_result = edge_index, norm

        edge_index, norm = self.cached_result

        return self.propagate(edge_index, x=x, norm=norm)


    def message(self, x_j, norm):
        return norm.view(-1, 1) * x_j if norm is not None else x_j

    def update(self, aggr_out):
        if self.bias is not None:
            aggr_out = aggr_out + self.bias
        return aggr_out

    def __repr__(self):
        return '{}({}, {})'.format(self.__class__.__name__, self.in_channels,
                                   self.out_channels)
                                   
class net(nn.Module):
    def __init__(self, in_channel, out_channel, node_num):
        super(net, self).__init__()
        self.node_num = node_num
        self.out_channel = out_channel
        self.layer1 = GCNConv(in_channel, 32, bias=True, normalize=True)
        self.layer2 = GCNConv(32, 64, bias=True, normalize=True)
        self.layer3 = GCNConv(64, out_channel, bias=True,normalize=True)
        self.lin1 = nn.Linear(out_channel * self.node_num, 1)
        
    def forward(self, x, edge_index):
        x = x.float()
        x = F.relu(self.layer1(x, edge_index),True)
        x = F.relu(self.layer2(x, edge_index),True)
        x = F.relu(self.layer3(x, edge_index),True)
        x = x.view(-1, self.out_channel*self.node_num)
        x = self.lin1(x)
        x = F.sigmoid(x)
        return x
        
def train_val_split(opdata, val_ratio=0.05, test_ratio=0.1, shuffle=False):
    """
    make the seperation of dataset, the ratio is 0.1 of dataset to be the 
    """
    num_obs = opdata.shape[0]
    n_v = int(math.floor(val_ratio * num_obs))
    n_t = int(math.floor(test_ratio * num_obs))
    if shuffle:
        opdata = np.random.shuffle(opdata)
    valdata = opdata[:n_v,:]
    testdata = opdata[n_v:n_v+n_t,:]
    traindata = opdata[n_v+n_t:,:]
    return traindata, valdata, testdata
    
    
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
traindata, valdata, testdata = train_val_split(data)
node_num = data.shape[1]
model = net(1, 64, node_num).to(device)
criterion = nn.BCELoss()
optimizer = optim.Adam(model.parameters(), lr=1e-4)


traindata_t = torch.from_numpy(traindata).unsqueeze(-1)
valdata_t = torch.from_numpy(valdata).unsqueeze(-1)
testdata_t = torch.from_numpy(testdata).unsqueeze(-1)
traintarget = Variable(traindata_t[:,-1,:]).to(device)
valtarget = Variable(valdata_t[:,-1,:]).to(device)
testtarget = Variable(testdata_t[:,-1,:]).to(device)
traindata_t = Variable(traindata_t[:,:-1,:]).to(device)
valdata_t = Variable(valdata_t[:,:-1,:]).to(device)
testdata_t = Variable(testdata_t[:,:-1,:]).to(device)

edge_index = torch.from_numpy(data_edge_index.values.T).to(device)

def train():
    model.train()
    
    optimizer.zero_grad()
    out = model(traindata_t, edge_index)
    loss = criterion(out, traintarget) 
    loss.backward()
    optimizer.step()
    return loss
    
def test():
    model.eval()
    loss = criterion(model(valdata_t, edge_index), valtarget)
    return loss.item()
    
for epoch in range(200):
    loss = train()
    log = 'Epoch: {:03d}, Train: {:.4f}, Test: {:.4f}'
    print(log.format(epoch+1, loss, test()))