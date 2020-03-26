# deep-geometric-pathway

## Genomic graph-learning LOG

#### March 25th

1. Implemented the `Node2vec` Embedding method to pre-train the feature of the pathway-graph
2. The saved embedding manifold representation is [node2vec_embed.pdf](GenomicData/results/node2vec-embedding/March-25th/node2vec_embed.pdf) . For now, the embedding is not sparse.
3. The `node2vec` method is in file `.node2vec_feature_represenation.iphb`
4. Added folder `./model` to save the `.pth` model parameters



#### March 24th

1. Implemented the `RidgeCV` with previous `pathway_Lasso` . The results are shown in `./results/Ridge`

2. Learned the `Node2Vec` embedding strategy

   

#### March 23rd

1. Implement `regression_hops_sampler()` function. The function can return an always $1$ batch size object loader. The `return` serves for the regression training without considering the geometric structure, but each training set of loader does from the `num_hops`. So, it is kind of fair for further geometric pathway training.  (**DONE**)
2. Normalized the data, fixed the possible $NAN$ issue
3. Implemented the `pathway_LassoCV.py` to compute the regression, all the plot results are shown as GROUND TRUTH vs PREDICTION. They are in `./results/Lasso` folder 

Example:

See detail is `./regression_hop_sampler.ipynb` ,`./pathway_LassoCV.ipynb`. Also some results figures: [['ABL1'].pdf](GenomicData/results/Lasso/['ABL1'].pdf)  [['APAF1'].pdf](GenomicData/results/Lasso/['APAF1'].pdf) 



#### March 22nd

1. Implement `cancer_data()` function, which operates the `Breast Caner` data in `./BreastCancer/Data_RNASeq2.mat` (**DONE**)
2. Implement `hops_sampler()` function that sample the overall *genomeic/pathway* info and return the `dataloader` in specific `batch_size` required. (**DONE**)

Example:

The detailed example is available in `./test_notebook/` and `./sampler_test.ipynb`



#### March 21st

 1. The pure text raw data `./GenomicData/Gene_DATA/sourcePathway.txt` processed method (**DONE**)
 2. The utility method `data_fetch()` class has been **DONE** 
 3. `data_fetch` is able to exclude the node type that is no necessary wanted

 Example:

 `$(LIST)` is optional input list type. e.g., `['protein', 'abstract']`

 ```python
 from utils.data_fetch as dffter
 dffter = data_fetch(argv_list=$(LIST))
 ```

 

 #### March 12th

 1. Have more clear idea of dataset
 2. Get the node information and pattern information
 3. Learned Regular expression(Regex)
 4. Convert rough data from `.mat` to `.txt` for future purpose
 5. The link and node info are in `./GenomicData/Gene_DATA/sourcePathway.txt`
 6. The further *genomic activity feature* is stored in `./GenomicData/Data_RNASeq2.mat` **NEED TO PROCESS LATER**

 The following is the *regex pattern* that can make extract information from source file really quick and efficient.

 ```regex
 // General node pattern
 ^([a-z]+)\s([\w\/\-()+]+)$
 
 // Link pattern extraction
 ^([\w\/\-()+]+)\s([\w\/\-()+]+)\s([\w\>\|-]+)$ 
 ```

