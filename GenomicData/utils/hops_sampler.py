import torch
import numpy as np
import pandas as pd
from torch_geometric.data import Data
from sklearn.preprocessing import LabelEncoder
import torch

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
        self.device = 'cuda' if torch.cuda.is_available() else 'cpu'
        
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
        self.edge_index = np.unique(self._add_selfLoop(edge_index), axis=1)
        self.node_i, self.node_j = self.edge_index
        for elem in self.data.remained_protein_id:
            level_node = self.node_i[np.in1d(self.node_j, elem)]
            level_node = np.setdiff1d(level_node, elem)
            self.save_ids.append(elem) if level_node.size > self.least_size  else self.save_ids

    def _add_selfLoop(self, edge_index):
        num_node = np.max(edge_index) + 1
        temp = np.vstack([np.arange(num_node)]*2)
        edge_index = np.concatenate([edge_index, temp],axis=1)
        return edge_index

    def relabel_edge(self, edge_index_ori, nodes_index):
        # The function relabel the original indexxed edge to be simple 0~num-1 index

        #--------------------------------------------------------------------------------------------
        # NOTE: REMODIFIED THIS FUNCTION TO MAKE THE edge_index_j starts always from 0 to the num_j-1
        #--------------------------------------------------------------------------------------------
        if np.all(np.isin(edge_index_ori, nodes_index)) is False:
            raise IndexError("The edge_index in not in the nodes_index, it is not able to relabel the edge.\n")

        # create a dictionary
        edge_i, edge_j = edge_index_ori
        n_ids = np.unique(edge_i)
        res_n = np.unique(edge_j)
        table_n_ids = {elem: idx for idx, elem in enumerate(n_ids)}
        table_res_n = {elem:idx for idx, elem in enumerate(res_n)}

        new_edge_i = np.array([table_n_ids.get(elem) for elem in edge_i])[np.newaxis,:]
        new_edge_j = np.array([table_res_n.get(elem) for elem in edge_j])[np.newaxis,:]

        return np.concatenate([new_edge_i,new_edge_j], axis=0)


    def _splits(self):
        #---------------------------------------------------------------------------------------
        ## Ver 3.0 Modify the sampling method from big sample to smaller sample step by step,
        # in the end we get the middle one
        #---------------------------------------------------------------------------------------

        # Start from here is the test of sub-sampling linkage based on num of hops
        # First, split the started node into desired batch_size
        self.batch_splits = torch.split(torch.tensor(self.save_ids), self.batch_size)
        # function that find all parent hops
        self.batched_node_list = []
        self.samples = []
        for batch in self.batch_splits:
            # the ids i want to keep in further sub-sampling part
            subset = Data()
            # get numpy final samples
            subset.target = batch.numpy()
            subset.this_batch_ids = batch.numpy()
            temp = batch.numpy()
            subset.dataflow = []
            subset.edge_index_t = []
            # Get the edge_index layer by layer
            # Test on 1 more hops everytimes
            for num in range(self.num_hops):
                subset.this_batch_ids = np.hstack([subset.this_batch_ids, self.node_i[np.in1d(self.node_j, subset.this_batch_ids)]])
                edge_indice = np.in1d(self.node_j, temp)
                temp = self.node_i[edge_indice]
                if num == 0:
                    this_edge_index = self.edge_index[:,edge_indice]
                else:
                    this_edge_index = np.unique(np.concatenate([this_edge_index, self.edge_index[:,edge_indice]],
                                                               axis=1), axis=1)
                subset.edge_index_t.append(this_edge_index)

            subset.edge_index_t = subset.edge_index_t[::-1]

            #-------------------------------------------------
            # overall propagation
            #-------------------------------------------------
            glb_edge_idx = subset.edge_index_t[0]
            block = Data()
            block.n_id = np.unique(glb_edge_idx.reshape(-1))
            block.edge_index_ori = np.unique(np.concatenate([glb_edge_idx, np.vstack([block.n_id]*2)], axis=-1), axis=1)
            alllist_table = {elem: idx for idx, elem in enumerate(block.n_id)}
            node_i = np.array([alllist_table.get(elem) for elem in block.edge_index_ori[0]])
            node_j = np.array([alllist_table.get(elem) for elem in block.edge_index_ori[1]])
            block.edge_index = torch.from_numpy(np.vstack([node_i,node_j])).long().to(self.device)
            subset.dataflow.append(block)
            #-------------------DONE overall propagation-------------

            for idx, elem in enumerate(subset.edge_index_t):
                block = Data()
                # Assign original edge_index
                block.edge_index_ori = elem[:,elem[1].argsort()]

                # Assign the n_id, and res_n_id
                block.n_id, block.res_n_id = elem
                block.n_id = np.unique(block.n_id)
                block.res_n_id = np.unique(block.res_n_id)

                block.edge_index = torch.from_numpy(self.relabel_edge(block.edge_index_ori, np.hstack([block.n_id, block.res_n_id]))).long().to(self.device)
                all_index_list = np.unique(block.edge_index_ori.reshape(-1))
                alllist_table = {elem:idx for idx,elem in enumerate(all_index_list)}
                block.res_n_id = [alllist_table.get(elem) for elem in block.res_n_id]
                block.size = (len(block.n_id), len(block.res_n_id))
                subset.dataflow.append(block)

            subset.this_batch_ids = np.unique(subset.this_batch_ids)
            self.batched_node_list.append(subset.this_batch_ids)
            subset.size_list = [len(obj.n_id) for obj in subset.dataflow]
            subset.size_list.append(len(subset.dataflow[-1].res_n_id))
            subset.property = ("DataFlow({}".format("{} -> "*(self.num_hops)) + "{})").format(*subset.size_list)
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