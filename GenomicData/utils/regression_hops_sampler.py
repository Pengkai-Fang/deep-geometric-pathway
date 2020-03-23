import pandas as pd
import utils.hops_sampler as hops_sampler
import numpy as np
from sklearn.preprocessing import LabelEncoder
from torch_geometric.data import Data
from .hops_sampler import hops_sampler
import torch

class regression_hops_sampler(hops_sampler):
    """
    This class in designed for the regression purpose of pathway dataset.
    The general funciton of `regression_sampler_hop` is almost the same as `hops_sampler`.
    Instead of considering the pathway structure, here we still sample the data based on the num_hops, but discard the non-protein genomic
    The `batch_size` here is not allowed to input manually and always set to $1$.

    Inputs: `pathway`: a pathway object
            `num_hops`: the number of hops wanted to sample
            `least_size(optional)`: the option that narrow down some nodes that do not have links
            
    Return: `regression_samples`: the loader of sub-graph without the consideration of graph relation
            `regression_samples.activ_free/cancer`: the training data
            `regression_samples.activ_free/cancer_target`:  the target value
    
    """
    def __init__(self, pathway, num_hops, least_size=5):
        self.data = pathway
        self.num_hops = num_hops
        self.batch_size = 1
        self.least_size = least_size
        
        self._setup()
        self._splits_protein()
        
        # From here, extract the remained graph based on batached_node_id-list
        le = LabelEncoder()
        le.fit(self.data.remained_protein)
        self.regression_samples = []
        
        for test_id, batch in zip(self.batch_splits, self.batched_node_list):
            self.regression_samples.append(self.regeression_sampler(test_id.numpy(), batch, le))
        
    def _splits_protein(self):
        # Start from here is the test of sub-sampling linkage based on num of hops
        # First, split the started node into desgined batch_size
        self.batch_splits = torch.split(torch.tensor(self.save_ids), self.batch_size)
        # function that find all parent hops
        self.batched_node_list = []
        for batch in self.batch_splits:
            # the ids i want to keep in further sub-sampling part
            this_batch_ids = batch.numpy()
            for num in range(self.num_hops):
                this_batch_ids = np.hstack([this_batch_ids, self.node_i[np.in1d(self.node_j, this_batch_ids)]])
            self.batched_node_list.append(np.intersect1d(self.data.remained_protein_id, np.unique(this_batch_ids)))
        
    def regeression_sampler(self, test_id, batch, le):
        pthway = Data()
        pthway.test_id = test_id
        pthway.operate_id = np.setdiff1d(batch, test_id)
        pthway.pthway_NameList = self.data.pthway_NameList.iloc[pthway.operate_id,:]
        if (self.data.pthway_NameList.iloc[pthway.test_id,:]['GenomeType'].values != 'protein'):
            raise IndexError('The test element is not the protein type. \t')
        
        if (np.all(pthway.pthway_NameList['GenomeType'].values == 'protein') is False):
            raise IndexError('Some of Related Genome ar not protein type.\t')
        
        # The pathway name list that exclude the target protein
        pthway.Genome_NameList = pthway.pthway_NameList['GenomeName'].values
        # Only the name of target protein
        pthway.Test_NameList = self.data.pthway_NameList.iloc[pthway.test_id,:]['GenomeName'].values
        
        # The overall trained dataset
        activ_id = le.transform(pthway.Genome_NameList)
        pthway.activ_free = self.data.activ_free[activ_id]
        pthway.activ_cancer = self.data.activ_cancer[activ_id]
        
        # the target value
        activ_id_test = le.transform(pthway.Test_NameList)
        pthway.activ_free_target = self.data.activ_free[activ_id_test]
        pthway.activ_cancer_target = self.data.activ_cancer[activ_id_test]
        
        return pthway