import pandas as pd
import numpy as np
from sklearn import preprocessing

class data_fetch():
    """
    This is the method that operates the source data file and stores processed format 
    into desired genometric learning format.
    
    **input**: `argv_list`(optional), argv_list is the list input that only take the possible node
    type that wanted to be exclude from the dataset.
    
    `obj.Namelist`: the entire processed namelist dataset with `'GenomeType'` 
    and `'GenomeName'` 2 columns.
    
    `obj.Edgelist`: the entire process edgelist with the labeled name. The dataframe has `start`, 'end' and
    `edge_type` information.
    
    Other stored variables: `edge_index`, `edge_class`, `node_class`, `node_all`, `edge_className`, `node_className`

    **Note**: the `node_all` is the numpy format node name with ligits labeled node class 
    
            For now, if the `linkskip` is `True`, the following deselect works: `abstract`, `rna`
    """
    def __init__(self, filepath, argv_list = [], linkskip=False):
        self.argv_list = argv_list
        self.linkskip = linkskip
        self.loadpath = filepath
        self.Namelist, self.Edgelist = self._load()
        if len(self.argv_list):
            self._remove(self.Namelist)
        
        self._edge_removal(self.Edgelist)
        
        self.edge_index, self.node_class, self.node_all = self._encoder()
        
    def _load(self):
        # load the data from text
        data_total = pd.read_csv(self.loadpath)

        # extract the info based on regex
        Namelist = data_total['TXT'].str.extract(r'^(?P<GenomeType>[a-z]+)\s(?P<GenomeName>[\w\/\-()+]+)$').dropna(axis=0).reset_index(drop=True)
        Edgelist = data_total['TXT'].str.extract(r'^(?P<edgeStart>[\w\/\-()+]+)\s(?P<edgeEnd>[\w\/\-()+]+)\s(?P<edgeType>[\w\>\|-]+)$').dropna(axis=0).reset_index(drop=True)
        
        return (Namelist, Edgelist)
    
    def _remove(self, Namelist):
        # Here we started to choose to discard some features/node types 
        le = preprocessing.LabelEncoder()
        le.fit(Namelist['GenomeType'])
        all_node_class = list(le.classes_)

        # Input error check
        for elem in self.argv_list:
            if elem not in all_node_class:
                raise NameError('The input argument {} is not defined in node class.'.format(elem))
        
        
        if self.linkskip:
            self._edge_removal(self.Edgelist)
            self.link_skipper(self.Namelist)
        else:
            self.Namelist = Namelist[~Namelist['GenomeType'].isin(self.argv_list)].reset_index(drop=True)
        
    
    def link_skipper(self, namelist):
        """
        The input should be self.Namelist object
        The function will skip linkage and re-link the remained object of the strategy
        """
        
        # step 1,2,3 - relabel everything by index
        le = preprocessing.LabelEncoder()
        namelist = namelist.drop_duplicates()
        namelist = namelist.sort_values(by='GenomeName', ignore_index=True)
        self.Namelist = namelist
        le.fit(namelist['GenomeName'].values)
        self.edge_index = le.transform(self.Edgelist.iloc[:,:2].values.reshape(-1)).reshape(-1,2)
        self.row, self.col = self.edge_index.T
        
        # step 4 - find unsatisfiled node
        del_index_list = namelist[namelist['GenomeType'].isin(self.argv_list)].index.values
        
        # step 5 - loop tracing forward and back
        for k in del_index_list:
            source_list = np.setdiff1d(self.row[self.col == k], k)
            target_list = np.setdiff1d(self.col[self.row == k], k)
            
            if ((len(source_list) == 0) or (len(target_list) == 0)):
                continue
            else:
                self._func3(source_list, target_list)
        
        # step 6 - deselect the namelist
        self.Namelist = self.Namelist[~self.Namelist['GenomeType'].isin(self.argv_list)].reset_index(drop=True)
        self.Edgelist = pd.DataFrame(le.inverse_transform(self.edge_index.reshape(-1)).reshape(-1,2), columns={'edgeStart', 'edgeEnd'})
        
    def _func1(self, m, target_list):
        for x in target_list:
            self._func2(m,x)
            
    def _func2(self, m, x):
        if self.Namelist.iloc[x,:]['GenomeType'] not in self.argv_list:
            temp = np.array([m, x]).reshape(-1,2)
            self.edge_index = np.concatenate([self.edge_index, temp], axis=0)
            return
        else:
            sub_target_list = np.setdiff1d(self.col[self.row == x], x)
            if len(sub_target_list) == 0:
                return
            else:
                self._func1(m, sub_target_list)
                
    def _func3(self, source_list, target_list):
        for m in source_list:
            if self.Namelist.iloc[m,:]['GenomeType'] not in self.argv_list:
                self._func1(m, target_list)
            else:
                sub_source_list = np.setdiff1d(self.row[self.col == m], m)
                if len(sub_source_list) == 0:
                    continue
                else:
                    self._func3(sub_source_list, target_list)
    
    
    def _edge_removal(self, Edgelist):
        # The beforehand check
        # Here we check the edge that does not showup in node list and drop these edges
        Namelist_l = list(self.Namelist['GenomeName'])
        Edgelist_l = list(Edgelist.iloc[:,0].values)
        Edgelist_ll = list(Edgelist.iloc[:,1].values)
        exclude_list = []
        for idx, (elem, elem2) in enumerate(zip(Edgelist_l, Edgelist_ll)):
            if ((elem not in Namelist_l) or (elem2 not in Namelist_l)):
                exclude_list.append(idx)

        self.Edgelist = Edgelist.drop(exclude_list).reset_index(drop=True)
    
    def _encoder(self):
        # Label the index of each edge
        # Make Edgelist --> edge_index

        # Label edge_index
        le = preprocessing.LabelEncoder()
        le.fit(self.Namelist['GenomeName'])
        edge_index = le.transform(self.Edgelist.iloc[:,:2].values.reshape(-1)).reshape(-1,2)
        
        # if we skip some links, that will make the the edgeType messed up. So, in this case, we dont consider the edge type anymore
        if self.linkskip is False:
            # Label edge_class
            le = preprocessing.LabelEncoder()
            le.fit(self.Edgelist['edgeType'])
            edge_class = le.transform(self.Edgelist['edgeType'])
            self.edge_className = list(le.classes_)
            self.edge_class = edge_class
        
        # Label node class
        le = preprocessing.LabelEncoder()
        le.fit(self.Namelist['GenomeType'])
        node_class = le.transform(self.Namelist['GenomeType'])
        self.node_className = list(le.classes_)
        
        # The combination of node info
        node_all = np.concatenate([self.Namelist['GenomeName'].values.reshape(-1,1), node_class.reshape(-1,1)], axis=1)
        
        return (edge_index, node_class, node_all)


