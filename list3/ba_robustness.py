import networkit as nk
import numpy as np
import copy

# generate barabasi-albert graph
def generate_ba_graph(n, avg_k):
    m = avg_k//2
    n0 = 2*m + 1
    G = nk.generators.BarabasiAlbertGenerator(m, n, n0).generate()
    return G


def size_of_giant_component(G):
    alg = nk.components.ConnectedComponents(G)
    alg.run()
    components = list(alg.getComponentSizes().values())
    return np.max(components)


def main():
    num_nodes = 10000
    L = 10
    avg_ks = [2,4]
    fs = np.linspace(0.0, 0.99, 100)  
    Pdegree = np.zeros((len(avg_ks), len(fs)))
    Pbetweenness = np.zeros((len(avg_ks), len(fs)))
    Pcloseness = np.zeros((len(avg_ks), len(fs)))
    Prandom = np.zeros((len(avg_ks), len(fs)))
    
    for i,avg_k in enumerate(avg_ks):
        for l in range(L):
            print("Running avg_k = {}, iter = {}".format(avg_k, l+1))
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
                for k,strategy in enumerate(strategies):
                    Gp = copy.deepcopy(G)
                    num_nodes_to_remove = int(num_nodes*f)
                    nodes_to_remove = strategy[:num_nodes_to_remove]
                    for node in nodes_to_remove:
                        Gp.removeNode(node)
                    if k == 0:
                        Prandom[i,j] += size_of_giant_component(Gp)/float(L)
                    elif k == 1:
                        Pdegree[i,j] += size_of_giant_component(Gp)/float(L)
                    elif k==2:
                        Pbetweenness[i,j] += size_of_giant_component(Gp/float(L)
                    elif k==3:
                        Pcloseness[i,j] += size_of_giant_component(Gp)/float(L)
                    else:
                        print("Error: unknown strategy")
                        
            # for j in range(len(fs)):
                # Pdegree[i,j] /= L
                # Prandom[i,j] /= L
                # Pbetweenness[i,j] /= L
                # Pcloseness[i,j] /= L
    
    # for each avg_k, save to file the fraction of nodes removed vs. size of giant component
    for i,avg_k in enumerate(avg_ks):
        np.savetxt("v5_ba_robustness_N_{}_L_{}_k_{}.dat".format(num_nodes, L, avg_k), np.array([fs, Prandom[i,:], Pdegree[i,:], Pbetweenness[i,:], Pcloseness[i,:]]).T)
          
        
        
    
if __name__ == "__main__":
    main()
