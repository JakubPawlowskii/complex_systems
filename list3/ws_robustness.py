import networkit as nk
import numpy as np
import copy

# generate Watts-Strogatz graph
def generate_ws_graph(n, avg_k, beta):
    G = nk.generators.WattsStrogatzGenerator(n, avg_k//2, beta).generate()
    return G


def remove_random_nodes(G, num_nodes_to_remove):
    nodes_to_remove = nk.graphtools.randomNodes(G, num_nodes_to_remove)
    for node in nodes_to_remove:
        G.removeNode(node)

# function removing given number of nodes with highest degree
def remove_nodes_with_highest_degree(G, num_nodes_to_remove):
    alg = nk.centrality.DegreeCentrality(G)
    alg.run()
    degrees = alg.ranking()
    nodes_to_remove = degrees[:num_nodes_to_remove]
    for node in nodes_to_remove:
        G.removeNode(node[0])
        
# function removing given number of nodes with highest betweenness centrality
def remove_nodes_with_highest_betweenness_centrality(G, num_nodes_to_remove):
    alg = nk.centrality.Betweenness(G)
    alg.run()
    betweenness = alg.ranking()
    nodes_to_remove = betweenness[:num_nodes_to_remove]
    for node in nodes_to_remove:
        G.removeNode(node[0])
        
# function removing given number of nodes with highest closeness centrality
def remove_nodes_with_highest_closeness_centrality(G, num_nodes_to_remove):
    alg = nk.centrality.Closeness(G, True,  nk.centrality.ClosenessVariant.Generalized) 
    alg.run()
    closeness = alg.ranking()
    nodes_to_remove = closeness[:num_nodes_to_remove]
    for node in nodes_to_remove:
        G.removeNode(node[0])


def size_of_giant_component(G):
    alg = nk.components.ConnectedComponents(G)
    alg.run()
    components = list(alg.getComponentSizes().values())
    return np.max(components)


def main():
    num_nodes = 10000
    L = 100
    beta = 0.01
    avg_ks = [2, 4]
    fs = np.linspace(0.0, 0.99, 100)
    Pdegree = np.zeros((len(avg_ks), len(fs)))
    Pbetweenness = np.zeros((len(avg_ks), len(fs)))
    Pcloseness = np.zeros((len(avg_ks), len(fs)))
    Prandom = np.zeros((len(avg_ks), len(fs)))
    
    for i,avg_k in enumerate(avg_ks):
        for j,f in enumerate(fs):
            print("Running avg_k = {}, f = {:.3f}".format(avg_k, f))
            for _ in range(L):
                Gdegree = generate_ws_graph(num_nodes, avg_k, beta)
                Grandom = copy.deepcopy(Gdegree)
                Gbetweenness = copy.deepcopy(Gdegree)
                Gcloseness = copy.deepcopy(Gdegree)
                num_nodes_to_remove = int(f * num_nodes)
                remove_nodes_with_highest_degree(Gdegree, num_nodes_to_remove)
                remove_random_nodes(Grandom, num_nodes_to_remove)
                remove_nodes_with_highest_betweenness_centrality(Gbetweenness, num_nodes_to_remove)
                remove_nodes_with_highest_closeness_centrality(Gcloseness, num_nodes_to_remove)
                
                Pdegree[i,j] += size_of_giant_component(Gdegree)
                Prandom[i,j] += size_of_giant_component(Grandom)
                Pbetweenness[i,j] += size_of_giant_component(Gbetweenness)
                Pcloseness[i,j] += size_of_giant_component(Gcloseness)
                
            Pdegree[i,j] /= L
            Prandom[i,j] /= L
            Pbetweenness[i,j] /= L
            Pcloseness[i,j] /= L    
    
    # for each avg_k, save to file the fraction of nodes removed vs. size of giant component
    for i,avg_k in enumerate(avg_ks):
        np.savetxt("wr_robustness_N_{}_L_{}_beta_{:.3f}_k_{}.dat".format(num_nodes, L, beta, avg_k), np.array([fs, Prandom[i,:], Pdegree[i,:], Pbetweenness[i,:], Pcloseness[i,:]]).T)
          
        
        
    
if __name__ == "__main__":
    main()
