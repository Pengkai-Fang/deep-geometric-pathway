# deep-geometric-pathway

## Genomic graph-learning LOG


#### **TODO**

1. ~~Figure out the same input for LASSO, check out the result~~
2. ~~Modify the model~~
3. Try modify the GraphSAGE, and try the results


#### April 1th-5th

1. Add graph-convoluational related files in folder `./model`
2. Implemented *GATGraph* based pathway inductive learning
3. Get the first graph learning result stored in `./results/Graph-conv/`
4. Redo the inductive sampling method. Now, the sampling is still based on the number of hops. However, the dataflow is from the big-global to the wanted predicted node step by step.
5. The main code is stored `./inductive-learning.py`

The graph-learning is better than LASSO prediction. Also, the R^2 score is 0.68, 0.65 for cancer free and cancer obtained patients. The both R^2 score is higher than LASSO prediction.

![graph-learning](https://github.com/Zaoyee/deep-geometric-pathway/blob/master/GenomicData/results/Graph-conv/%5B'ABL1'%5D.png)

![LASSO](https://github.com/Zaoyee/deep-geometric-pathway/blob/master/GenomicData/results/Lasso/%5B'ABL1'%5D.pdf)


#### March 29th

1. Add `utils.visualize.show_pth()` function, that is able to print out the specific sampled pathway

   The detail of how to use this function is in `./test_notebook/visualize.ipynb`. Example visualization [here](GenomicData/readme-figs/ALAS1_hops2.pdf).

#### March 28th

1. Add `skip_link` option for `data_fetch()` function in order to realize the skip link in some situations instead of deleting every links that conncets to unwanted node type. For now, the `skip_link` is only work for `abstract` and `rna`. 
2. Add $R^2$ score computation of regression results. Updated the code in `pathway_Lasso_w_Ridge.ipynb`. 

$R^2$ is defined as: $R^2 = 1 - {\sum(y_{pred} - y_{true})^2 \over \sum(y_{true} - \bar{y_{true}})^2}$, the best situation is that $R^2 = 1$ 



The following is the idea of the option `skip_link = True`, when we want to exclude `abstract` and `rna`

![link_skip](https://github.com/Zaoyee/deep-geometric-pathway/blob/master/GenomicData/readme-figs/link_skip.png)

#### March 27th

**TODO**

1. visualize the sampled pathway compared to the .m code
2. ~~show average consistency for regrssion~~ ( **DONE** )
3. verify the pathway input with different parameter setting
4. ~~when exlcude the linkage type, do we get rid of the rest remained linkage type that connected~~ ( **DONE** )



#### March 26th

1. Update the `utils.hop_sampler()` to be able to handle the `DataFlow` structure.
2. The sample method will return a flow way to train the data or train the sample once for all
3. Be able to **visualize** the dataflow of each sample, and handle the *batch_size* 
4. `hops_sampler` added `label_all` option, default is `False` , which allows the class compute the overall relabeled every node and edge.

Example: The following is some `dataflow` visualization in case of sample hops is `5`. 

Note: The results below miss one flow for each data flow. **FIXED**

```reStructuredText
 Data(dataflow=[5], property=DataFlow(3 <- 20 <- 172 <- 727 <- 2720), size_list=[6]),
 Data(dataflow=[5], property=DataFlow(3 <- 61 <- 431 <- 1946 <- 4938), size_list=[6]),
 Data(dataflow=[5], property=DataFlow(3 <- 29 <- 218 <- 1210 <- 4184), size_list=[6]),
 Data(dataflow=[5], property=DataFlow(3 <- 21 <- 200 <- 1151 <- 3873), size_list=[6]),
 Data(dataflow=[5], property=DataFlow(3 <- 46 <- 316 <- 1620 <- 4866), size_list=[6]),
 Data(dataflow=[5], property=DataFlow(3 <- 24 <- 187 <- 1276 <- 4455), size_list=[6]),
 Data(dataflow=[5], property=DataFlow(3 <- 33 <- 185 <- 794 <- 2711), size_list=[6]),
```



####  March 25th

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

