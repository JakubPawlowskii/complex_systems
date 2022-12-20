import networkit as nk
import numpy as np
import copy

# generate barabasi-albert graph
def generate_ba_graph(n, avg_k):
    G = nk.generators.BarabasiAlbertGenerator(avg_k, n).generate()
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
    avg_ks = [1, 2, 4]
    fs = np.linspace(0.0, 0.99, 100)  
    Pdegree = np.zeros((len(avg_ks), len(fs)))
    Pbetweenness = np.zeros((len(avg_ks), len(fs)))
    Pcloseness = np.zeros((len(avg_ks), len(fs)))
    Prandom = np.zeros((len(avg_ks), len(fs)))
    
    for i,avg_k in enumerate(avg_ks):
        for _ in range(L):

            G = generate_ba_graph(num_nodes, avg_k)
            strategies = []
            
            random = np.random.permutation(num_nodes)
            strategies.append(random)
            
            alg = nk.centrality.DegreeCentrality(G)
            alg.run()
            degrees = alg.ranking()
            strategies.append([x[0] for x in degrees])
            
            alg = nk.centrality.Betweenness(G)
            alg.run()
            betweenness = alg.ranking()
            strategies.append([x[0] for x in betweenness])
            
            alg = nk.centrality.Closeness(G, True,  nk.centrality.ClosenessVariant.Generalized) 
            alg.run()
            closeness = alg.ranking()
            strategies.append([x[0] for x in closeness])

            for j,f in enumerate(fs):
                print("Running avg_k = {}, f = {:.3f}".format(avg_k, f))
                for k,strategy in enumerate(strategies):
                    Gp = copy.deepcopy(G)
                    num_nodes_to_remove = int(num_nodes*f)
                    nodes_to_remove = strategy[:num_nodes_to_remove]
                    for node in nodes_to_remove:
                        Gp.removeNode(node)
                    if k == 0:
                        Prandom[i,j] += size_of_giant_component(Gp)
                    elif k == 1:
                        Pdegree[i,j] += size_of_giant_component(Gp)
                    elif k==2:
                        Pbetweenness[i,j] += size_of_giant_component(Gp)
                    elif k==3:
                        Pcloseness[i,j] += size_of_giant_component(Gp)
                    else:
                        print("Error: unknown strategy")
                        
            for j in range(len(fs)):    
                Pdegree[i,j] /= L
                Prandom[i,j] /= L
                Pbetweenness[i,j] /= L
                Pcloseness[i,j] /= L    
    
    # for each avg_k, save to file the fraction of nodes removed vs. size of giant component
    for i,avg_k in enumerate(avg_ks):
        np.savetxt("v3_ba_robustness_N_{}_L_{}_k_{}.dat".format(num_nodes, L, avg_k), np.array([fs, Prandom[i,:], Pdegree[i,:], Pbetweenness[i,:], Pcloseness[i,:]]).T)
          
        
        
    
if __name__ == "__main__":
    main()
