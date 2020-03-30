import networkx as nx
import matplotlib.pyplot as plt


def show_pth(pathway_obj, *argv):
    """
    The simple function that print out the pathway dataset
    The different node types are visualized by different color
    
    **Inputs:**   `pathway_obj`: the pathwayway `object` that called by function `cancer_data()`
                  `batch`: the `batch` of the returned hops_sampler object
                          Example: hops_samples_obj = hops_sampler(pathway = data, 
                                                                   batch_size = 1, 
                                                                   num_hops = 2)
                                   batch = hops_samples_obj.samples[0]
    
    """
    if len(argv) == 1:
        batch = argv[0]
    elif len(argv) == 2:
        dataflow = argv[0]
        batch = dataflow[argv[1]]
    else:
        raise TypeError("The input number is incorrect.")
    
    g = nx.DiGraph()
    # we can put these input some global variable parts
    pathway_info_namelist = pathway_obj.pthway_NameList
    for block in batch.dataflow:
        temp_edge_index = block.edge_index_ori
        temp_edge = pathway_info_namelist.iloc[temp_edge_index.reshape(-1),:]['GenomeName'].values.reshape(2,-1).T
        temp_edge = ([x,y] for x, y in zip(list(temp_edge[:,0]), list(temp_edge[:,1])))
        g.add_edges_from(temp_edge)

    node_color = [pathway_obj.node_class[pathway_obj.pthway_NameList[pathway_obj.pthway_NameList['GenomeName'] == name].index] for name in g.nodes()]    
    plt.figure(figsize=(30,30))
    nx.draw_planar(g, with_labels = True, node_size=1500, node_color=node_color)
    plt.show()