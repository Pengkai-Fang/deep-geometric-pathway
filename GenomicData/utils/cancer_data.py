import scipy.io as sio
import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
from .data_fetch import data_fetch


class cancer_data():
    """
    This cancer_data class handle the both pathway data and cancer data at the same time
    
    Inputs : `cancerPATH` : the cancerdata path in `./BreastCancer`
                 `pthwayPATH` : the pathway data path in `./Gene_DATA`
                 
                 `exclude_type`: the specific pathway node type does not want to include
    
    Parameters:
             `Edgelist`, `activ_free`, `activ_cancer`, `edge_index`, `pthway_NameList`
             `all_elem_className`, `edge_class`, `edge_className`, `node_class`, `node_className`
    
    """
    def __init__(self, cancerPATH, pthwayPATH, exlcude_type=[]):
        self.cancerPATH = cancerPATH
        self.pthwayPATH = pthwayPATH
        self.argv_list = exlcude_type
        
        self._load()
        self._removal()
        self._shuffle_w_construct()
        self._Label_encoder()
        
    def _load(self):
        raw_data = sio.loadmat(self.cancerPATH)
        self.pthway = data_fetch(filepath=self.pthwayPATH,
                                 argv_list=self.argv_list)
        
        self.namelist_free = raw_data['SampleList0']
        self.namelist_cancer = raw_data['SampleList1']
        self.geneNamelist = raw_data['gNamePureListRNASeq2'][0]
        self.activ_free = raw_data['geneRNASeq2RawMatrix0']
        self.activ_cancer = raw_data['geneRNASeq2RawMatrix1']
        
        self.protein_list = self.pthway.Namelist
        
        if 'protein' not in self.pthway.node_className:
            raise NameError('The pathway data does not include the protein info')
        
        
    def _removal(self):
        # Obtained the data that only inside of protein_list    
        self.protein_list = self.protein_list[self.protein_list['GenomeType'] == 'protein']['GenomeName'].values
        self.geneNamelist = [np.array2string(self.geneNamelist[i]).replace("['",'').replace("']",'') for i in range(self.geneNamelist.size)]
        self.geneNamelist = np.array(self.geneNamelist).astype(np.object)

        included_gene = np.in1d(self.geneNamelist, self.protein_list)

        # Get two input data in the end, 
        # TODO reshuffle required 
        self.activ_free = self.activ_free[included_gene]
        self.activ_cancer = self.activ_cancer[included_gene]
        
        # Here we get the new genenNamelist that only include the contained protein name
        self.geneNamelist = self.geneNamelist[included_gene]

        included_protein = np.in1d(self.protein_list, self.geneNamelist)

        self.remained_protein = self.protein_list[included_protein]
        otherpth_list = self.pthway.Namelist[self.pthway.Namelist['GenomeType'] != 'protein']['GenomeName'].values
        remained_protein_all = np.concatenate([self.remained_protein, otherpth_list])
        
        # save the temporary edgelist
        Edgelist = self.pthway.Edgelist
        Namelist_l = list(remained_protein_all)
        Edgelist_l = list(Edgelist.iloc[:,0].values)
        Edgelist_ll = list(Edgelist.iloc[:,1].values)
        exclude_list = []
        for idx, (elem, elem2) in enumerate(zip(Edgelist_l, Edgelist_ll)):
            if ((elem not in Namelist_l) or (elem2 not in Namelist_l)):
                exclude_list.append(idx)

        self.Edgelist = Edgelist.drop(exclude_list).reset_index(drop=True)
        
        
    def _shuffle_w_construct(self):
        # shuffle the data to be the same order of both protein and genome
        tracked_index = np.argsort(self.geneNamelist)
        self.geneNamelist = self.geneNamelist[tracked_index]
        self.activ_free = self.activ_free[tracked_index]
        self.activ_cancer = self.activ_cancer[tracked_index]
        
        # reconstruct nodeset
        self.remained_protein = np.sort(self.remained_protein)
        
        protein_df = pd.DataFrame({'GenomeType': 'protein', 
                                   'GenomeName': self.remained_protein})
        self.pthway_NameList = pd.concat([protein_df, self.pthway.Namelist[self.pthway.Namelist['GenomeType'] != 'protein']], 
                               axis=0, 
                               ignore_index=True)
        
    def _Label_encoder(self):
        # Start to index label
        self.pthway_NameList = self.pthway_NameList.sort_values(by='GenomeName').reset_index(drop=True)
        # There are some fucking duplicated data !!!!!! FUCK, TAKE ME SO MUCH TIME TO FIND OUT
        self.pthway_NameList = self.pthway_NameList[~self.pthway_NameList.duplicated()]
        
        le = LabelEncoder()
        le.fit(self.pthway_NameList['GenomeName'].values)
        self.edge_index = le.transform(self.Edgelist.iloc[:,:2].values.reshape(-1)).reshape(-1,2)
        self.all_elem_className = list(le.classes_)
        
        # Make the id and matrix matched
        self.remained_protein_id = le.transform(self.remained_protein)
        
        # Label edge_class
        le = LabelEncoder()
        le.fit(self.Edgelist['edgeType'])
        self.edge_class = le.transform(self.Edgelist['edgeType'])
        self.edge_className = list(le.classes_)

        # Label node class
        le = LabelEncoder()
        le.fit(self.pthway_NameList['GenomeType'])
        self.node_class = le.transform(self.pthway_NameList['GenomeType'])
        self.node_className = list(le.classes_) 
