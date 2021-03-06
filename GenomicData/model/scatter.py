import numpy as np

def scatter(num_nodes, genome, genome_idxs, embed):
    num_patients, num_genome, num_feat = genome.shape
    activ_x = np.zeros((num_patients, num_nodes, num_feat))
    activ_x[:, genome_idxs, :] = genome
    # randomize
    # other_idxs= np.setdiff1d(np.arange(num_nodes), genome_idxs)
    # activ_x[:, other_idxs, :] = np.random.randn(num_patients,len(other_idxs),1)
    if embed is not None:
        embed = embed.data.numpy()
        embed = np.hstack([embed] * num_patients).reshape(num_patients, num_nodes, embed.shape[1])
        activ_x = np.concatenate([activ_x, embed], axis=2)
    return activ_x