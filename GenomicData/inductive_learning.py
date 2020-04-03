import utils.hops_sampler as hops_sampler
import utils.cancer_data as cancer_data
import torch.nn as nn
from model.gat_conv import GATConv
import torch.nn.functional as F
from torch_geometric.nn import Node2Vec
import torch
import numpy as np
from torch import optim


# Load both data and operate by class pathway
pathwayPATH = './Gene_DATA/sourcePathway.txt'
cancerPATH = './BreastCancer/Data_RNASeq2.mat'

# load the overall pathway and cancer data in object
data = cancer_data.cancer_data(pthwayPATH=pathwayPATH, cancerPATH=cancerPATH)


# sample the protein for the regression problem
hops_samples_obj = hops_sampler.hops_sampler(pathway=data,
                                             batch_size=1,
                                             num_hops=2)



class GATNet(nn.Module):
    def __init__(self, in_channels, out_channels):
        super(GATNet, self).__init__()
        self.conv1 = GATConv(in_channels, 8, heads=8, dropout=0.6, node_dim=1)
        self.conv2 = GATConv(8 * 8, out_channels, heads=1, concat=True,
                             dropout=0.6, node_dim=1)

    def forward(self, x, dataflow):
        block = dataflow[0]
        xt = x[:, block.n_id,:]
        xt = F.elu(
            self.conv1((xt, xt[:, block.res_n_id, :]),
                       block.edge_index,
                       size=block.size))

        #         x = F.dropout(x, p=0.6, training=self.training)
        #         block = data_flow[-2]
        #         x = F.elu(
        #             self.conv2((x, x[:,block.res_n_id,:]),
        #                        block.edge_index,
        #                       size=block.size))
        xt = F.dropout(xt, p=0.6, training=self.training)
        block = dataflow[1]
        xt = self.conv2((xt, xt[:, block.res_n_id, :]),
                        block.edge_index,
                        size=block.size)
        return xt


num_nodes = 15867
pre_embed = Node2Vec(num_nodes, embedding_dim=1, walk_length=20,
                     context_size=10, walks_per_node=10)
pre_embed.load_state_dict(torch.load('./model-dict/node2vec-predembedding/node2vec_1dim.pth'))
pre_embed.eval()
pre_embed_x = pre_embed(torch.arange(num_nodes))

free_x_patient = data.activ_free.T[:, :, np.newaxis]
cancer_x_patient = data.activ_cancer.T[:, :, np.newaxis]


genome_idxs = data.pthway_NameList[data.pthway_NameList['GenomeType'] == 'protein'].index


def scatter(num_nodes, genome, genome_idxs, embed):
    num_patients, num_genome, num_feat = genome.shape
    activ_x = np.zeros((num_patients, num_nodes, num_feat))
    activ_x[:, genome_idxs, :] = genome
    embed = embed.data.numpy()
    embed = np.hstack([embed] * num_patients).T[:, :, np.newaxis]
    activ_x = np.concatenate([activ_x, embed], axis=2)
    return activ_x


free_x_patient_all = scatter(num_nodes, free_x_patient, genome_idxs, pre_embed_x)
cancer_x_patient_all = scatter(num_nodes, cancer_x_patient, genome_idxs, pre_embed_x)


device = 'cpu' if torch.cuda.is_available() else 'cpu'
model = GATNet(2, 1).to(device)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)



def data_splits(samples, idx, train_ratio, test_ratio):
    flow = samples[idx]
    target_index = flow.target
    total_x = np.concatenate([free_x_patient_all, cancer_x_patient_all], axis=0)
    total_target = total_x[:, target_index, 0]
    flow.x = torch.from_numpy(total_x)
    flow.y = torch.from_numpy(total_target).unsqueeze(-1)

    num = flow.y.shape[0]
    row_index = torch.randperm(num)

    x = flow.x[row_index]
    y = flow.y[row_index]
    flow.train_x = x[:int(np.floor(num * train_ratio)), :, :].to(device)
    flow.train_y = y[:int(np.floor(num * train_ratio))].to(device)

    flow.test_x = x[int(np.floor(num * train_ratio)):, :, :].to(device)
    flow.test_y = y[int(np.floor(num * train_ratio)):].to(device)

    return flow

# acquire the data
dataflow = data_splits(hops_samples_obj.samples, 0, 0.8, 0.2).to(device)


def test():
    model.eval()
    pred = model(dataflow.test_x.float().to(device), dataflow.dataflow)
    loss = criterion(pred, dataflow.test_y.float().to(device))
    return loss.item()


# start to training
# Here, we only do test example on `ABL1`

for epoch in range(1000):
    model.train()
    optimizer.zero_grad()
    pred = model(dataflow.train_x.float().to(device), dataflow.dataflow)
    loss = criterion(pred, dataflow.train_y.float().to(device))
    loss.backward()
    optimizer.step()
    test_loss = test()
    print("Epoch {} : Train Loss {} Test Loss {}".format(epoch + 1, loss.item(), test_loss))
