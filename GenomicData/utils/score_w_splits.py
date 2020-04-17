import numpy as np
import torch
from sklearn.metrics import r2_score
import torch.nn as nn

device = 'cuda' if torch.cuda.is_available() else 'cpu'
criterion = nn.SmoothL1Loss()

def data_splits(samples, idx, test_idx, splits_id, all_patients, class_patients):
    #extract the specific flow
    assert test_idx <= np.max(splits_id)

    flow = samples[idx]
    target_index = flow.target
    true_target = all_patients[:, target_index, 0]
    flow.org = all_patients
    flow.true_target = true_target
    flow.x = torch.from_numpy(all_patients)
    flow.y = torch.from_numpy(true_target).unsqueeze(-1)

    flow.train_x = flow.x[splits_id != test_idx]
    #radomize
    # num_pataients, num_node, dim=flow.train_x.shape

    # prevent overfit and cheat on predication
    flow.train_x[:,target_index,0] = 0
    flow.train_x = flow.train_x.to(device)
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

def test_R_score(test_pred, target, train, class_ary):
    # compute the R_score in train set
    if train:
        return r2_score(target.cpu().data.numpy().reshape(-1), test_pred.cpu().data.numpy().reshape(-1))
    else:
        target = target
        free_pred = test_pred[class_ary == 0].reshape(-1)
        cancer_pred = test_pred[class_ary == 1].reshape(-1)
        free_true = target[class_ary == 0].reshape(-1)
        cancer_true = target[class_ary == 1].reshape(-1)
        return r2_score(free_true, free_pred), r2_score(cancer_true, cancer_pred),\
               criterion(torch.from_numpy(free_pred), torch.from_numpy(free_true)).item(), \
               criterion(torch.from_numpy(cancer_pred), torch.from_numpy(cancer_true)).item()
