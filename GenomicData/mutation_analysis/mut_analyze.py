import numpy as np
import sympy as sp

def mut_analyze(flow, activ_pred, true_target, mut_loader, cancer_data_obj):
    """
    Start the analysis of every flow related to the mutation matrix
    Inputs:
            `flow`: The extracted flow after training. The detail of flow is listed in function `utils.score_w_splits.data_splits()`
            `activ_pred`: The prediction of the activity level of every cancer obtained patient
            `mut_loader`: The mutation_loader object. Source to `mutation_analysis.mutation_loader`
            `cancer_data_obj`: The cancer data object that is called by function `utils.cancer_data`

    # TODO NOTE THIS FUNCTION DID NOT CONSIDER THE TARGET GENOME SELF
    """
    data = cancer_data_obj.pthway_NameList
    dataflow = flow.dataflow[::-1]
    mut_mat = mut_loader.mut_matrix

    # step 0 - preprocess the mut matrix to have column index of each non-nan value
    temp_mut_mat = np.vstack([np.arange(mut_mat.shape[1])] * mut_mat[0])
    mut_mat[mut_mat != np.nan] = temp_mut_mat[mut_mat != np.nan]

    # ~TODO~ The flow index does not match on genome_matrix (DONE)
    # step 1 - compute the inconsistency score and only take the cancer patient samples
    score = inconsistent_score(activ_pred, true_target)[111:]
    assert score.shape[0] == 1079

    # step 2 - get current flow available index
    target_id = flow.target
    sampled_idx = np.unique(flow.this_batch_ids) # include the target id
    num_layer = len(dataflow) # since the first element of the dataflow is the global aggregation layer

    # get the each layer index separately
    temp = target_id
    layer_storage = [temp]
    temp_sum = layer_storage.copy()
    for num in range(num_layer-1):
        temp = np.setdiff1d(np.unique(dataflow[num]), temp_sum)
        temp_sum = np.hstack([temp, temp_sum])
        layer_storage.append(temp)

    # step 3.0
    # # overall splits the inconsistency
    # total_num_mutated = np.sum(np.isin(mut_mat, sampled_idx), axis=1)
    # score = score / total_num_mutated
    alpha = solve_alpha(num_layer)

    # step 3.1 - based on the `layer_storage` compute the most affected genome
    idx_store = []
    score_store = []
    for layer in layer_storage[1:]:
        # For each sample, get the mutated genome number
        tf_table = np.isin(mut_mat, layer)
        this_layer_mutated_num = np.sum(tf_table, axis=1)
        this_layer_mutated_idx = [mut_mat[idx, tf_table[idx]] for idx in range(mut_mat.shape[0])]

        # split the score
        layer_score = score/this_layer_mutated_num * (alpha**layer)
        layer_score_samples = [np.ones(elem.shape[0])*layer_score[idx] for idx, elem in enumerate(this_layer_mutated_idx)]
        this_layer_mutated_idx = np.hstack(this_layer_mutated_idx)
        layer_score_samples = np.hstack(layer_score_samples)
        idx_store.append(this_layer_mutated_idx)
        score_store.append(layer_score_samples)

    idx_store = np.hstack(idx_store)
    score_store = np.hstack(score_store)

    # step 4 - scatter mean and sort the list
    idx_store, score_store = scatter_mean(score_store, idx_store)
    sorted_indice = np.argsort(score_store)
    idx_store = idx_store[sorted_indice]
    score_store = score_store[sorted_indice]

    # step 5 - show the result
    prt_str = "For test genome {}: \n\tThe most affected Genome is {} that has the highest inconsistency score {:.3f}" \
              "\n\t The second affected Genome is {} that has the inconsistency score {:.3f}"
    prt_str = prt_str.format(data.iloc[target_id]['GenomeName'].values,
                             data.iloc[idx_store[0]]['GenomeName'].values,
                             score_store[0], data.iloc[idx_store[1]]['GenomeName'].values,
                             score_store[1])
    print(prt_str)

def scatter_mean(score, idx):
    idx_list = np.unique(idx)
    scattered_score = np.zeros(idx_list.shape)
    for i, elem in enumerate(idx_list):
        scattered_score[i] = score[idx == elem].mean()

    return idx_list, scattered_score

def solve_alpha(num_layer):
    assert num_layer < 6
    x = sp.symbols('x')
    ls = [x**idx for idx in range(1, num_layer+1)]
    solution = sp.solve(sum(ls)-1, x)
    for elem in solution:
        try:
            if float(elem) > 0:
                return float(elem)
        except:
            continue
    return 0.5

def inconsistent_score(pred, true):
    return np.abs(pred - true) / np.mean(np.abs(pred-true))

