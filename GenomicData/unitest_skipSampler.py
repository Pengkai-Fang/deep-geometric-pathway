import utils.hops_sampler as hops_sampler
import utils.cancer_data as cancer_data

# Load both data and operate by class pathway
pathwayPATH = './Gene_DATA/sourcePathway.txt'
cancerPATH = './BreastCancer/Data_RNASeq2.mat'

# load the overall pathway and cancer data in object
data = cancer_data.cancer_data(pthwayPATH=pathwayPATH, cancerPATH=cancerPATH)


# sample the protein for the regression problem
hops_samples_obj = hops_sampler.hops_sampler(pathway=data,
                                             batch_size=1,
                                             num_hops=3,
                                             keep_type=['protein'])


hops_samples_obj.samples
