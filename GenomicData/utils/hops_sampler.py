import torch
import numpy as np
import pandas as pd
from torch_geometric.data import Data
from sklearn.preprocessing import LabelEncoder

class hops_sampler(object):
    """
    This function does the sampling method for the dataset based on the number of hops defined
    
    Inputs: `pathway`: a pathway object
            `num_hops`: the number of hops wanted to sample
            `batch_size`: the batch_size of returned dataloader
            `least_size(optional)`: the option that narrow down some nodes that do not have links
            
    Return: `hops_samples`: the loader of sub-graph
    
    """
    def __init__(self, pathway, num_hops, batch_size, least_size=5):
        self.data = pathway
        self.num_hops = num_hops
        self.batch_size = batch_size
        self.least_size = least_size
        self.hops_samples = []
        
        self._setup()
        self._splits()
        
        # From here, extract the remained graph based on batached_node_id-list
        le = LabelEncoder()
        le.fit(self.data.remained_protein)
        
        for sub_batch in self.batched_node_list:
            self.hops_samples.append(self.sampler_generater(sub_batch, le))
        
    def _setup(self):
        # Create deep_geometric data object
        edge_index = self.data.edge_index.T
        
        # the remained id that does have at least one hop
        # Based on computation, we only have 1435 protein nodes satisfy the requirement
        self.save_ids = []
        self.edge_index = edge_index
        self.node_i, self.node_j = self.edge_index
        for elem in self.data.remained_protein_id:
            level_node = self.node_i[np.in1d(self.node_j, elem)]
            level_node = np.setdiff1d(level_node, elem)
            self.save_ids.append(elem) if level_node.size > self.least_size  else self.save_ids
            
        
    def _splits(self):
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
            self.batched_node_list.append(np.unique(this_batch_ids))
        
    def sampler_generater(self, batch, le):
        """
        This function passes batch index number to obtained trained object
        """
        deep_pthway = Data()
        newpthway_Namelist = self.data.pthway_NameList.iloc[batch,:].reset_index(drop=True)
        deep_pthway.genome_Namelist = newpthway_Namelist[newpthway_Namelist['GenomeType'] == 'protein']['GenomeName'].values
        activ_id = le.transform(deep_pthway.genome_Namelist)
        deep_pthway.activ_free = self.data.activ_free[activ_id]
        deep_pthway.activ_cancer = self.data.activ_cancer[activ_id]

        deep_pthway.pth_Namelist = newpthway_Namelist
        Edgelist = self.data.Edgelist
        Namelist_l = list(newpthway_Namelist['GenomeName'].values)
        Edgelist_l = list(Edgelist.iloc[:,0].values)
        Edgelist_ll = list(Edgelist.iloc[:,1].values)
        exclude_list = []
        for idx, (elem, elem2) in enumerate(zip(Edgelist_l, Edgelist_ll)):
            if ((elem not in Namelist_l) or (elem2 not in Namelist_l)):
                exclude_list.append(idx)

        newpthway_Edgelist = Edgelist.drop(exclude_list).reset_index(drop=True)
        deep_pthway.Edgelist = newpthway_Edgelist

        le2 = LabelEncoder()
        le2.fit(deep_pthway.pth_Namelist['GenomeName'].values)
        deep_pthway.edge_index = le2.transform(deep_pthway.Edgelist.iloc[:,:2].values.reshape(-1)).reshape(-1,2)
        deep_pthway.all_elem_className = list(le2.classes_)

        # Label edge_class
        le2 = LabelEncoder()
        le2.fit(deep_pthway.Edgelist['edgeType'])
        deep_pthway.edge_class = le2.transform(deep_pthway.Edgelist['edgeType'])
        deep_pthway.edge_className = list(le2.classes_)

        # Label node class
        le2 = LabelEncoder()
        le2.fit(deep_pthway.pth_Namelist['GenomeType'])
        deep_pthway.node_class = le2.transform(deep_pthway.pth_Namelist['GenomeType'])
        deep_pthway.node_className = list(le2.classes_) 
        
        return deep_pthway