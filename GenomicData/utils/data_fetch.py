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
    """
    def __init__(self, filepath, argv_list = []):
        self.argv_list = argv_list
        self.loadpath = filepath
        self.Namelist, self.Edgelist = self._load()
        if len(self.argv_list):
            self._remove(self.Namelist)
        
        self._edge_removal(self.Edgelist)
        self.edge_index, self.edge_class, self.node_class, self.node_all = self._encoder()
        
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

        # Start to exclude the node type
        self.Namelist = Namelist[~Namelist['GenomeType'].isin(self.argv_list)].reset_index(drop=True)
        
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

        # Label edge_class
        le = preprocessing.LabelEncoder()
        le.fit(self.Edgelist['edgeType'])
        edge_class = le.transform(self.Edgelist['edgeType'])
        self.edge_className = list(le.classes_)
        
        # Label node class
        le = preprocessing.LabelEncoder()
        le.fit(self.Namelist['GenomeType'])
        node_class = le.transform(self.Namelist['GenomeType'])
        self.node_className = list(le.classes_)
        
        # The combination of node info
        node_all = np.concatenate([self.Namelist['GenomeName'].values.reshape(-1,1), node_class.reshape(-1,1)], axis=1)
        
        return (edge_index, edge_class, node_class, node_all)


