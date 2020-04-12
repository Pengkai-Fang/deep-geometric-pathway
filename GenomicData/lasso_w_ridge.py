#%%

from sklearn.linear_model import LassoCV
import pandas as pd
import utils.hops_sampler as hops_sampler
import numpy as np
from sklearn.preprocessing import LabelEncoder
from torch_geometric.data import Data
import utils.hops_sampler as hops_sampler
import utils.cancer_data as pathway
import utils.regression_hops_sampler as regression_sampler_hop
import torch
import matplotlib.pyplot as plt
from sklearn.linear_model import RidgeCV

#%%

# Load both data and operate by class pathway
pathwayPATH = './Gene_DATA/sourcePathway.txt'
cancerPATH = './BreastCancer/Data_RNASeq2.mat'

# load the overall pathway and cancer data in object
data = pathway.cancer_data(pthwayPATH=pathwayPATH, cancerPATH=cancerPATH)

# sample the protein for the regression problem
regression_samplesLoader = regression_sampler_hop.regression_hops_sampler(pathway=data,
                                                 num_hops=2)

#%% md

#### L1_regularization AKA Lasso

#%%

data_loader = regression_samplesLoader.regression_samples
for idx, data_elem in enumerate(data_loader):
    x_1 = data_elem.activ_free.T
    x_2 = data_elem.activ_cancer.T
    X = np.concatenate([x_1, x_2], axis=0)
    y = np.hstack([data_elem.activ_free_target.reshape(-1), data_elem.activ_cancer_target.reshape(-1)])
    reg = LassoCV(cv=5, random_state=0, fit_intercept=True).fit(X, y)
    data_loader[idx].predict_free = reg.predict(x_1)
    data_loader[idx].predict_cancer = reg.predict(x_2)
    data_loader[idx].score_free = reg.score(x_1, data_elem.activ_free_target.reshape(-1))
    data_loader[idx].score_cancer = reg.score(x_2, data_elem.activ_cancer_target.reshape(-1))

#%%

# Showup for some predictions and real value
for idx, obj in enumerate(data_loader):
    plt.figure(figsize=(10,10))
    caner_loss = np.mean((obj.activ_cancer_target.reshape(-1) - obj.predict_cancer) ** 2)
    free_loss = np.mean((obj.activ_free_target.reshape(-1) - obj.predict_free) ** 2)
    plt.scatter(obj.activ_cancer_target.reshape(-1), obj.predict_cancer, s=5, c='r',
                label='cancer R^2 = {:.4f} Loss = {:.4f}'.format(obj.score_cancer, caner_loss))
    plt.scatter(obj.activ_free_target.reshape(-1), obj.predict_free, s=5, c='b',
                label='free R^2 = {:.4f} Loss = {:.4f}'.format(obj.score_free, free_loss))
    plt.legend(fontsize=13)
    # plot the criterion line
    temp = np.hstack([obj.activ_free_target.reshape(-1),
                      obj.activ_cancer_target.reshape(-1)])
    xx = np.linspace(temp.min(), temp.max(), 10000)
    yy = xx
    plt.plot(xx, yy, 'g')
    plt.title(data.pthway_NameList.iloc[obj.test_id,:]['GenomeName'].values)
    plt.xlabel('The truth ground')
    plt.ylabel('The LASSO prediction')
    plt.savefig('./results/Lasso/{}.pdf'.format(data.pthway_NameList.iloc[obj.test_id,:]['GenomeName'].values))
    plt.show()

#%% md

#### L2_regularization AKA Ridge

#%%

data_loader = regression_samplesLoader.regression_samples
for idx, data_elem in enumerate(data_loader):
    x_1 = data_elem.activ_free.T
    x_2 = data_elem.activ_cancer.T
    X = np.concatenate([x_1, x_2], axis=0)
    y = np.hstack([data_elem.activ_free_target.reshape(-1), data_elem.activ_cancer_target.reshape(-1)])
    comb = np.concatenate([X,y.reshape(-1,1)],axis=1)
    np.random.shuffle(comb)
    X = comb[:,:-1]
    y = comb[:,-1].reshape(-1)
    reg = RidgeCV(cv=5, fit_intercept=True, alphas=[1e-3, 1e-2, 1e-1, 1]).fit(X, y)
    data_loader[idx].predict_free = reg.predict(x_1)
    data_loader[idx].predict_cancer = reg.predict(x_2)
    data_loader[idx].score_free = reg.score(x_1, data_elem.activ_free_target.reshape(-1))
    data_loader[idx].score_cancer = reg.score(x_2, data_elem.activ_cancer_target.reshape(-1))

#%%

# Showup for some predictions and real value
for idx, obj in enumerate(data_loader):
    plt.figure(figsize=(10,10))
    caner_loss = np.mean((obj.activ_cancer_target.reshape(-1) - obj.predict_cancer)**2)
    free_loss = np.mean((obj.activ_free_target.reshape(-1) - obj.predict_free) ** 2)
    plt.scatter(obj.activ_cancer_target.reshape(-1), obj.predict_cancer, s=5, c='r',
                label='cancer R^2 = {:.4f} Loss = {:.4f}'.format(obj.score_cancer, caner_loss))
    plt.scatter(obj.activ_free_target.reshape(-1), obj.predict_free, s=5, c='b',
                label='free R^2 = {:.4f} Loss = {:.4f}'.format(obj.score_free, free_loss))
    plt.legend(fontsize=13)
    # plot the criterion line
    temp = np.hstack([obj.activ_free_target.reshape(-1),
                      obj.activ_cancer_target.reshape(-1)])
    xx = np.linspace(temp.min(), temp.max(), 10000)
    yy = xx
    plt.plot(xx, yy, 'g')
    plt.title(data.pthway_NameList.iloc[obj.test_id,:]['GenomeName'].values)
    plt.xlabel('The ground truth')
    plt.ylabel('The Ridge prediction')
    plt.savefig('./results/Ridge/{}.pdf'.format(data.pthway_NameList.iloc[obj.test_id,:]['GenomeName'].values))
    plt.show()