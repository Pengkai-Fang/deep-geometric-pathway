import utils.hops_sampler as hops_sampler
import utils.cancer_data as cancer_data
import torch.nn as nn
from model.gat_conv import GATConv
import torch.nn.functional as F
from torch_geometric.nn import Node2Vec
import torch
import numpy as np
from torch import optim
from model.scatter import scatter
from sklearn.metrics import r2_score

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
        self.conv1 = GATConv(in_channels, 128, heads=12, node_dim=1)
        # self.conv2 = GATConv(16 * 6, 128, heads=6, node_dim=1)
        self.conv3 = GATConv(128*12 , 1, heads=1, node_dim=1, concat=True)

    def forward(self, x, dataflow):
        block = dataflow[0]
        xt = x[:, block.n_id,:]
        xt = F.relu(
            self.conv1((xt, xt[:, block.res_n_id, :]),
                       block.edge_index,
                       size=block.size))

        block = dataflow[1]
        # xt = F.relu(self.conv2((xt, xt[:, block.res_n_id, :]),
        #                 block.edge_index,
        #                 size=block.size))
        #
        #
        # block = dataflow[2]
        xt = self.conv3((xt, xt[:, block.res_n_id, :]),
                        block.edge_index,
                        size=block.size)
        return xt


num_nodes = data.pthway_NameList.shape[0]
# pre_embed = Node2Vec(num_nodes, embedding_dim=16, walk_length=20,
#                      context_size=10, walks_per_node=10)
# pre_embed.load_state_dict(torch.load('./model-dict/node2vec-predembedding/node2vec_16dim.pth'))
# pre_embed.eval()
# pre_embed_x = pre_embed(torch.arange(num_nodes))
pre_embed_x = None
free_x_patient = data.activ_free.T[:, :, np.newaxis]
cancer_x_patient = data.activ_cancer.T[:, :, np.newaxis]


genome_idxs = data.pthway_NameList[data.pthway_NameList['GenomeType'] == 'protein'].index

free_x_patient_all = scatter(num_nodes, free_x_patient, genome_idxs, pre_embed_x)
cancer_x_patient_all = scatter(num_nodes, cancer_x_patient, genome_idxs, pre_embed_x)


device = 'cuda' if torch.cuda.is_available() else 'cpu'
model = GATNet(1, 1).to(device)
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.04)

# define the permutation here
all_patients  = np.concatenate([free_x_patient_all, cancer_x_patient_all], axis=0)
class_patients = np.hstack([np.zeros(free_x_patient.shape[0]), np.ones(cancer_x_patient.shape[0])])
permutation_idx = np.random.permutation(all_patients.shape[0])
all_patients, class_patients = all_patients[permutation_idx], class_patients[permutation_idx]
# generate cv
cv = 5
splits_id = np.hstack([np.hstack([np.arange(cv)]*int(np.floor(all_patients.shape[0] / cv))),
                       np.arange(int(np.floor(all_patients.shape[0] % cv)))])


def data_splits(samples, idx, test_idx, splits_id, all_patients, class_patients):
    #extract the specific flow
    assert test_idx <= np.max(splits_id)

    flow = samples[idx]
    target_index = flow.target
    true_target = all_patients[:, target_index, 0]
    flow.true_target = true_target
    flow.x = torch.from_numpy(all_patients)
    flow.y = torch.from_numpy(true_target).unsqueeze(-1)

    flow.train_x = flow.x[splits_id != test_idx].to(device)
    flow.train_y = flow.y[splits_id != test_idx].to(device)
    flow.train_class = class_patients[splits_id != test_idx]

    flow.test_x = flow.x[splits_id == test_idx].to(device)
    flow.test_y = flow.y[splits_id == test_idx].to(device)
    flow.test_class = class_patients[splits_id == test_idx]

    return flow

def train_R_score(flow, train_pred):
    # compute the R_score in train set
    free_pred = train_pred[flow.train_class == 0].cpu().data.numpy().reshape(-1)
    cancer_pred = train_pred[flow.train_class == 1].cpu().data.numpy().reshape(-1)
    free_true = flow.train_y[flow.train_class == 0].cpu().data.numpy().reshape(-1)
    cancer_true = flow.train_y[flow.train_class == 1].cpu().data.numpy().reshape(-1)
    return r2_score(free_true, free_pred), r2_score(cancer_true, cancer_pred)

def test_R_score(test_pred, target):
    # compute the R_score in train set
    target = target
    free_pred = test_pred[class_patients == 0].reshape(-1)
    cancer_pred = test_pred[class_patients== 1].reshape(-1)
    free_true = target[class_patients == 0].reshape(-1)
    cancer_true = target[class_patients == 1].reshape(-1)
    return r2_score(free_true, free_pred), r2_score(cancer_true, cancer_pred)


def test(dataflow):
    model.eval()
    pred = model(dataflow.test_x.float().to(device), dataflow.dataflow)
    loss = criterion(pred, dataflow.test_y.float().to(device))
    return pred, loss.item()

# start to training
# Here, we only do test example on `ABL1`
def train(dataflow):
    for epoch in range(5000):
        model.train()
        optimizer.zero_grad()
        pred = model(dataflow.train_x.float().to(device), dataflow.dataflow)
        loss = criterion(pred, dataflow.train_y.float().to(device))
        R_score = "Free R^2 score {} Cancer R^2 score {}".format(*train_R_score(dataflow, pred))
        loss.backward()
        optimizer.step()
        test_pred, test_loss = test(dataflow)
        print("Epoch {} : Train Loss {} Test Loss {} {}".format(epoch + 1, loss.item(), test_loss, R_score))
    return test_pred


test_pred_all = np.zeros(all_patients.shape[0])
for idx in range(cv):
# acquire the data
    flow = data_splits(hops_samples_obj.samples, 0,
                       idx, splits_id, all_patients, class_patients).to(device)

    test_pred = train(flow).cpu().data.numpy()
    test_pred_all[splits_id == idx] = test_pred.reshape(-1)

print("Across all {} folds, the overall R^2 score is \n\tFree {} Cancer {}".format(cv, *test_R_score(test_pred_all, flow.true_target)))


import matplotlib.pyplot as plt

plt.figure(figsize=(10,10))
plt.scatter(flow.true_target[class_patients == 1].reshape(-1), test_pred_all[class_patients == 1].reshape(-1), s=5, c='r')
plt.scatter(flow.true_target[class_patients == 0].reshape(-1), test_pred_all[class_patients == 0].reshape(-1), s=5, c='b')
plt.legend(fontsize=13)
# plot the criterion line
temp = flow.true_target
xx = np.linspace(temp.min(), temp.max(), 10000)
yy = xx
plt.plot(xx, yy, 'g')
plt.title(data.pthway_NameList.iloc[flow.target,:]['GenomeName'].values)
plt.xlabel('The truth ground')
plt.ylabel('The LASSO prediction')
plt.show()



