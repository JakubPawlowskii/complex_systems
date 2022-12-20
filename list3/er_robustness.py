import networkit as nk
import numpy as np
import copy

def generate_er_graph(n, avg_k):
    p = float(avg_k) / n
    G = nk.generators.ErdosRenyiGenerator(n, p).generate()
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
    avg_ks = [0.5,1,2,4]
    fs = np.linspace(0.0, 0.99, 100)
    Pdegree = np.zeros((len(avg_ks), len(fs)))
    Pbetweenness = np.zeros((len(avg_ks), len(fs)))
    Pcloseness = np.zeros((len(avg_ks), len(fs)))
    Prandom = np.zeros((len(avg_ks), len(fs)))
    
    # number of nodes to remove at each step
    num_nodes_to_remove = int(num_nodes * (fs[1] - fs[0]))
    node_idx_to_remove = [np.arange(i*num_nodes_to_remove,(i+1)*num_nodes_to_remove) for i in range(num_nodes//num_nodes_to_remove - 1)]
    # print(node_idx_to_remove)
    # print(fs)
    print(len(node_idx_to_remove))
    # print(len(fs))
    for i,avg_k in enumerate(avg_ks):
        for _ in range(L):
            
            Gdegree = generate_er_graph(num_nodes, avg_k)
            alg = nk.centrality.DegreeCentrality(Gdegree)
            alg.run()
            degrees = alg.ranking()
            # print(degrees)
            Grandom = copy.deepcopy(Gdegree)
            
            Gbetweenness = copy.deepcopy(Gdegree)    
            alg = nk.centrality.Betweenness(Gbetweenness)
            alg.run()
            betweenness = alg.ranking()
            
            Gcloseness = copy.deepcopy(Gdegree)
            alg = nk.centrality.Closeness(Gcloseness, True,  nk.centrality.ClosenessVariant.Generalized) 
            alg.run()
            closeness = alg.ranking()

            Pdegree[i,0] += size_of_giant_component(Gdegree)
            Prandom[i,0] += size_of_giant_component(Grandom)
            Pbetweenness[i,0] += size_of_giant_component(Gbetweenness)
            Pcloseness[i,0] += size_of_giant_component(Gcloseness)
            for j,node_idxs in enumerate(node_idx_to_remove):
                print("Running avg_k = {}, f = {:.3f}".format(avg_k, fs[j+1]))
                
                remove_deg = [degrees[idx][0] for idx in node_idxs]
                for node in remove_deg:
                    Gdegree.removeNode(node)
            
                for _ in range(num_nodes_to_remove):
                    node = nk.graphtools.randomNode(Grandom)
                    Grandom.removeNode(node)
                        
                remove_bet = [betweenness[idx][0] for idx in node_idxs]
                for node in remove_bet:
                    Gbetweenness.removeNode(node)
                
                remove_clo = [closeness[idx][0] for idx in node_idxs]
                for node in remove_clo:
                    Gcloseness.removeNode(node)
        
                Pdegree[i,j] += size_of_giant_component(Gdegree)
                Prandom[i,j] += size_of_giant_component(Grandom)
                Pbetweenness[i,j] += size_of_giant_component(Gbetweenness)
                Pcloseness[i,j] += size_of_giant_component(Gcloseness)
            
            for j in range(len(fs)):    
                Pdegree[i,j] /= L
                Prandom[i,j] /= L
                Pbetweenness[i,j] /= L
                Pcloseness[i,j] /= L    
    
    # for each avg_k, save to file the fraction of nodes removed vs. size of giant component
    for i,avg_k in enumerate(avg_ks):
        np.savetxt("v2_er_robustness_N_{}_L_{}_k_{}.dat".format(num_nodes, L, avg_k), np.array([fs, Prandom[i,:], Pdegree[i,:], Pbetweenness[i,:], Pcloseness[i,:]]).T)
          
        
    
if __name__ == "__main__":
    main()
