import networkit as nk
import sys
# cmd args
# 1. network type
#   - 'cg' : Complete graph
#   - 'er' : Erdos-Renyi graph
#   - 'ws' : Watts-Strogatz graph
# 2. number of nodes
# 3. average degree of node
# 4. probability of rewriting (only for ws)
# 5. number of generated graphs

NUMBER_OF_NODES = 'number of nodes'
NETWORK_TYPE = 'network type'
EDGE_PROBABILITY = 'edge probability'
AVERAGE_DEGREE = 'average degree'
REWIRE_PROBABILITY = 'rewire probability'


def generator(params):
    if params[NETWORK_TYPE] == 'cg':
        return nk.generators.ErdosRenyiGenerator(params[NUMBER_OF_NODES], 1.0)
    if params[NETWORK_TYPE] == 'er':
        return nk.generators.ErdosRenyiGenerator(params[NUMBER_OF_NODES], params[EDGE_PROBABILITY])
    if params[NETWORK_TYPE] == 'ws':
        return nk.generators.WattsStrogatzGenerator(params[NUMBER_OF_NODES], params[AVERAGE_DEGREE] / 2,
                                                    params[REWIRE_PROBABILITY])
    raise AttributeError("Unsupported network type.")


arg_list = sys.argv
graph_properties = dict()
network_type = arg_list[1]
name = ''
metis_name = ''
if network_type == 'cg':
    name = 'complete graph'
    metis_name = 'complete'
    graph_properties[NETWORK_TYPE] = 'cg'
    graph_properties[NUMBER_OF_NODES] = int(arg_list[2])
if network_type == 'er':
    name = 'erdos-renyi graph'
    metis_name = 'erdos_renyi'
    graph_properties[NETWORK_TYPE] = 'er'
    graph_properties[NUMBER_OF_NODES] = int(arg_list[2])
    graph_properties[AVERAGE_DEGREE] = int(arg_list[3])
    graph_properties[EDGE_PROBABILITY] = float(arg_list[3]) / (float(arg_list[2]) - 1)
if network_type == 'ws':
    name = 'watts-strogatz graph'
    metis_name = 'watts_strogatz'
    graph_properties[NETWORK_TYPE] = 'ws'
    graph_properties[NUMBER_OF_NODES] = int(arg_list[2])
    graph_properties[AVERAGE_DEGREE] = int(arg_list[3])
    graph_properties[REWIRE_PROBABILITY] = float(arg_list[4])

number_of_graphs = int(arg_list[-1])

try:
    gen = generator(graph_properties)
except AttributeError:
    sys.exit('No such network. Aborting...')

for i in range(number_of_graphs):
    filename = metis_name + str(i + 1) + '.metis'
    print(filename)
    number_of_components = 2
    while number_of_components > 1:    
        G = nk.Graph(gen.generate())
        components = nk.components.ConnectedComponents(G)
        components.run()
        number_of_components = components.numberOfComponents()
    
    nk.writeGraph(G, filename, nk.Format.METIS)
    with open(filename, "r+") as f:
        s = f.read()
        f.seek(0)
        f.write("#"+name+"\n" + s)

print(graph_properties)
print('Number of generated graphs = ', number_of_graphs)
