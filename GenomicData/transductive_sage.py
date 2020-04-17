import utils.hops_sampler as hops_sampler
import utils.cancer_data as cancer_data
import torch.nn as nn
from model.sage_conv import SAGEConv
import torch.nn.functional as F
import torch
import numpy as np
from torch import optim
from model.scatter import scatter
from utils.score_w_splits import data_splits, test_R_score, train_R_score
import matplotlib.pyplot as plt

# Load both data and operate by class pathway
pathwayPATH = './Gene_DATA/sourcePathway.txt'
cancerPATH = './BreastCancer/Data_RNASeq2.mat'

# load the overall pathway and cancer data in object
data = cancer_data.cancer_data(pthwayPATH=pathwayPATH, cancerPATH=cancerPATH)


# sample the protein for the regression problem
hops_samples_obj = hops_sampler.hops_sampler(pathway=data,
                                             batch_size=1,
                                             num_hops=2)


class SageNet(nn.Module):
    def __init__(self, in_channels, out_channels, concat=True):
        super(SageNet, self).__init__()

        self.conv1 = SAGEConv(in_channels, 64, normalize=True, concat=concat, node_dim=1)
        self.conv2 = SAGEConv(64, 256, normalize=True, concat=concat, node_dim=1)
        self.conv3 = SAGEConv(256, 512, normalize=True, concat=concat, node_dim=1)

        self.lin1 = nn.Sequential(
            nn.Linear(512*77, 1024),
            nn.ReLU(True),
            nn.Linear(1024, 218),
            nn.ReLU(True),
            nn.Linear(218, out_channels)
        )

        # self.conv3 = SAGEConv(36, 216, normalize=False, concat=False, node_dim=1)
        # self.conv4 = SAGEConv(216, out_channels, normalize=False, concat=False, node_dim=1)

    def forward(self, x, dataflow):
        block = dataflow[0]
        xt = x[:,block.n_id,:]
        xt = F.relu(
            self.conv1(xt, block.edge_index)
        )
        # xt = F.dropout(xt, p=0.5, training=self.training)
        xt = F.relu(
            self.conv2(xt, block.edge_index)
        )
        xt = F.relu(
            self.conv3(xt, block.edge_index)
        )
        xt = xt.view(xt.shape[0],-1)
        xt = self.lin1(xt)

        # xt = F.dropout(xt, p=0.5, training=self.training)
        # block = dataflow[1]
        # xt = F.relu(
        #     self.conv3((xt, None),
        #                block.edge_index,
        #                size=block.size,
        #                res_n_id = block.res_n_id))
        # # xt = F.dropout(xt, p=0.5, training=self.training)
        # block = dataflow[2]
        # xt = self.conv4((xt, None),
        #                 block.edge_index,
        #                 size=block.size,
        #                 res_n_id = block.res_n_id)

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


device = 'cpu' if torch.cuda.is_available() else 'cpu'
model = SageNet(1, 1).to(device)
# criterion = nn.MSELoss()
criterion = nn.SmoothL1Loss()
optimizer = optim.Adam(model.parameters(), lr=0.01)

# define the permutation here
all_patients  = np.concatenate([free_x_patient_all, cancer_x_patient_all], axis=0)
class_patients = np.hstack([np.zeros(free_x_patient.shape[0]), np.ones(cancer_x_patient.shape[0])])
permutation_idx = np.random.permutation(all_patients.shape[0])
all_patients, class_patients = all_patients[permutation_idx], class_patients[permutation_idx]

# overall edge_index
edge_index = torch.from_numpy(hops_samples_obj.edge_index).long().to(device)

# generate cv
cv = 10
splits_id = np.hstack([np.hstack([np.arange(cv)]*int(np.floor(all_patients.shape[0] / cv))),
                       np.arange(int(np.floor(all_patients.shape[0] % cv)))])


def test(dataflow):
    model.eval()
    pred = model(dataflow.test_x.float().to(device), dataflow.dataflow).unsqueeze(-1)
    loss = criterion(pred, dataflow.test_y.float().to(device))
    score = test_R_score(pred, dataflow.test_y.float().to(device), True, class_patients)
    return pred, loss.item(), score

# start to training
# Here, we only do test example on `ABL1`
def train(dataflow):
    x = dataflow.train_x.float().to(device)
    y = dataflow.train_y.float().to(device)
    for epoch in range(500):
        model.train()
        optimizer.zero_grad()
        pred = model(x, dataflow.dataflow).unsqueeze(-1)
        loss = criterion(pred, y)
        R_score = "Free R^2 score {:.2f} Cancer R^2 score {:.2f} ".format(*train_R_score(dataflow, pred))
        loss.backward()
        optimizer.step()
        test_pred, test_loss, score = test(dataflow)
        print("Epoch {} : Train Loss {:.6f} Test Loss {:.6f} {} Test_score {:.2f}".format(epoch + 1, loss.item(),
                                                                        test_loss, R_score, score))
        del pred, loss, test_loss, score
        torch.cuda.empty_cache()
    return test_pred


test_pred_all = np.zeros(all_patients.shape[0])
for idx in range(cv):
# acquire the data
    model = SageNet(1, 1).to(device)
    criterion = nn.SmoothL1Loss()
    optimizer = optim.Adam(model.parameters(), lr=0.002)
    flow = data_splits(hops_samples_obj.samples, 0,
                       idx, splits_id, all_patients, class_patients).to(device)

    test_pred = train(flow).cpu().data.numpy()
    test_pred_all[splits_id == idx] = test_pred.reshape(-1)
    del model, criterion, optimizer
    torch.cuda.empty_cache()
print("Across all {} folds, the overall R^2 score is \n\tFree {} Cancer {}".format(cv, *test_R_score(test_pred_all, flow.true_target, False, class_patients)))

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

# step last - invert the permutation to match on the original input
test_pred_all = test_pred_all[np.argsort(permutation_idx)]
