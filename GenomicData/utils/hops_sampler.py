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
    def __init__(self, pathway, num_hops, batch_size, least_size=5, label_all=False):
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
        
        if label_all:
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
            
    def relabel_edge(self, edge_index_ori, nodes_index):
        # The function relabel the original indexxed edge to be simple 0~num-1 index
        if np.all(np.isin(edge_index_ori, nodes_index)) is False:
            raise IndexError("The edge_index in not in the nodes_index, it is not able to relabel the edge.\n")
        le = LabelEncoder()
        le.fit(nodes_index)
        return le.transform(edge_index_ori.reshape(-1)).reshape(2,-1)
    
    def _splits(self):
        # Start from here is the test of sub-sampling linkage based on num of hops
        # First, split the started node into desgined batch_size
        self.batch_splits = torch.split(torch.tensor(self.save_ids), self.batch_size)
        # function that find all parent hops
        self.batched_node_list = []
        self.samples = []
        for batch in self.batch_splits:
            # the ids i want to keep in further sub-sampling part
            this_batch_ids = batch.numpy()
            temp = batch.numpy()
            subset = Data()
            subset.dataflow = []
            for num in range(self.num_hops):
                # Create dataflow structure
                block = Data()
                block.res_n_id = temp
                link_indice = np.in1d(self.node_j, temp)
                block.edge_index_ori = self.edge_index[:,link_indice]               
                temp = self.node_i[link_indice]
                block.n_id = temp
                block.edge_index = self.relabel_edge(block.edge_index_ori, np.hstack([block.n_id, block.res_n_id]))
                subset.dataflow.append(block)
                
                # overall variable
                this_batch_ids = np.hstack([this_batch_ids, self.node_i[np.in1d(self.node_j, this_batch_ids)]])
            self.batched_node_list.append(np.unique(this_batch_ids))
            subset.size_list = [len(obj.res_n_id) for obj in subset.dataflow]
            subset.size_list.append(len(subset.dataflow[-1].n_id))
            subset.property = ("DataFlow({}".format("{} <- "*(self.num_hops)) + "{})").format(*subset.size_list)
            self.samples.append(subset)
        
    def sampler_generater(self, batch, le):
        """
        This function passes batch index number to obtained trained object
        """
        deep_pthway = Data()
        deep_pthway.batch = batch
        
        # find out the edge_index of original pathway index
        row, col = self.edge_index
        row = list(row)
        col = list(col)
        exclude_list = []
        for idx, (a, b) in enumerate(zip(row,col)):
            if ((a not in list(batch)) or (b not in list(batch))):
                exclude_list.append(idx)
        deep_pthway.edge_index_oriIndexxed = np.delete(self.edge_index, exclude_list, 1)
        
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