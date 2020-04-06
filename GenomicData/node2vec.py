import torch
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from torch_geometric.nn import Node2Vec
import utils.hops_sampler as hops_sampler
import utils.cancer_data as pathway


# Load both data and operate by class pathway
pathwayPATH = './Gene_DATA/sourcePathway.txt'
cancerPATH = './BreastCancer/Data_RNASeq2.mat'

# load the overall pathway and cancer data in object
data = pathway.cancer_data(pthwayPATH=pathwayPATH, cancerPATH=cancerPATH)

# sample the protein for the regression problem
hops_samples_obj = hops_sampler.hops_sampler(pathway=data,
                                             batch_size=2,  # must>1
                                             num_hops=3,
                                             label_all=True)

hops_sample_loader = hops_samples_obj.hops_samples


device = 'cuda' if torch.cuda.is_available() else 'cpu'
num_nodes = hops_samples_obj.data.pthway_NameList.shape[0]
model = Node2Vec(num_nodes, embedding_dim=16, walk_length=20,
                 context_size=10, walks_per_node=10)
model = model.to(device)
optimizer = torch.optim.Adam(model.parameters(), lr=0.01)


hops_sample_loader[0].edge_index_oriIndexxed



def train():
    model.train()
    total_loss = 0
    for subset in hops_sample_loader:
        optimizer.zero_grad()
        batch = torch.from_numpy(subset.batch).long().to(device)
        edge_index = torch.from_numpy(subset.edge_index_oriIndexxed).long().to(device)
        loss = model.loss(edge_index, batch)
        loss.backward()
        optimizer.step()
        total_loss += loss.item()
    return total_loss / len(hops_sample_loader)



for epoch in range(1, 500):
    loss = train()
    print('Epoch: {:02d}, Loss: {:.4f}'.format(epoch, loss))


model.eval()
with torch.no_grad():
    z = model(torch.arange(num_nodes, device=device))
    z = TSNE(n_components=2).fit_transform(z.cpu().numpy())
    y = data.node_class


# def plot_points(colors):
#     model.eval()
#     with torch.no_grad():
#         z = model(torch.arange(num_nodes, device=device))
#         z = TSNE(n_components=3).fit_transform(z.cpu().numpy())
#         y = data.node_class

colors = [
    '#ffc0cb', '#bada55', '#008080', '#420420', '#7fe5f0', '#065535', '#ffd700'
]
# ax = Axes3D(plt.figure())
plt.figure(figsize=(8, 8))
for i in range(7):
    plt.scatter(z[y == i, 0], z[y == i, 1], s=20, color=colors[i])
plt.axis('off')
plt.savefig('./results/node2vec-embedding/March-25th/node2vec_embed.pdf')
plt.show()



torch.save(model.state_dict(), './model-dict/node2vec-predembedding/node2vec_16dim.pth')