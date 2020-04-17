import numpy as np
import scipy.io as sio

class mutation_Loader(object):
    """
    This class is designed for loading mutation object as same as the cancer_data
    The function is similar to utils.cancerdata class
    """

    def __init__(self, cancerData_obj, mutation_PATH):
        self.cancerData_obj = cancerData_obj
        self.path = mutation_PATH
        self._load()
        self._slice_w_reshuffle()

    def _load(self):
        raw_data = sio.loadmat(self.path)
        self.mut_matrix = raw_data['MutMatrix'][:,:1079]
        self.mut_namelist = raw_data['geneNamelist'][0]
        self.mut_namelist = [np.array2string(self.mut_namelist[i]).replace("['", '').replace("']", '') for i in
                             range(self.mut_namelist.size)]
        self.mut_namelist = np.array(self.mut_namelist).astype(np.object)

    def _slice_w_reshuffle(self):
        kept_index = np.in1d(self.mut_namelist, self.cancerData_obj.geneNamelist)
        self.mut_matrix = self.mut_matrix[kept_index]
        self.mut_namelist = self.mut_namelist[kept_index]
        tracked_index = np.argsort(self.mut_namelist )
        self.mut_namelist = self.mut_namelist [tracked_index]
        self.mut_matrix  = self.mut_matrix [tracked_index]
        assert np.all(self.mut_namelist == self.cancerData_obj.geneNamelist)

        # TODO make the mut_matrix have the all types genomes
        genome_idxs = self.cancerData_obj.pthway_NameList[self.cancerData_obj.pthway_NameList['GenomeType'] == 'protein'].index
        total_mut = np.empty((self.mut_matrix.shape[1], self.cancerData_obj.pthway_NameList.shape[0]))
        total_mut[:] = np.nan
        total_mut[:, genome_idxs] = self.mut_matrix.T
        self.mut_matrix = total_mut