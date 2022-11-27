#include <iostream>

#include "../utils/graph.hpp"

int main()
{
    auto g = graph{};
    


    std::cout << "Number of nodes: " << g.get_number_of_nodes() << std::endl;
    std::cout << "Number of edges: " << g.get_number_of_links() << std::endl;
    std::cout << "Average degree: " << g.get_average_degree() << std::endl;
    std::cout << "Average clustering: " << g.get_average_clustering() << std::endl;

    g.from_adj_list_to_nodes_and_edges();


    return 0;
}