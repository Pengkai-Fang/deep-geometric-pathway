# deep-geometric-pathway
 
 ## Genomic graph-learning LOG
 
 
 
 #### March 22nd
 
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
 
