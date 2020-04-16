import utils.cancer_data as cancer_data
import mutation_analysis.mutation_loader as mutation_loader


# Load both data and operate by class pathway
pathwayPATH = './Gene_DATA/sourcePathway.txt'
cancerPATH = './BreastCancer/Data_RNASeq2.mat'
mutPATH = './BreastCancer/Data_MUT.mat'

# load the overall pathway and cancer data in object
data = cancer_data.cancer_data(pthwayPATH=pathwayPATH, cancerPATH=cancerPATH)

mut_data = mutation_loader.mutation_Loader(cancerData_obj=data, mutation_PATH=mutPATH)
