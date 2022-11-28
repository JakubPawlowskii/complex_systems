#include <iostream>

#include "../utils/graph.hpp"

int main()
{
    auto g = graph{};
    g.watts_strogatz(1000, 150, 0.000001);
    // g.k_regular_graph(1000, 10);

    std::cout << "Number of nodes: " << g.get_number_of_nodes() << std::endl;
    std::cout << "Number of edges: " << g.get_number_of_links() << std::endl;
    std::cout << "Average degree: " << g.get_average_degree() << std::endl;
    std::cout << "Average clustering: " << g.get_average_clustering() << std::endl;

    g.from_adj_list_to_nodes_and_edges();

    dijkstra searcher{g, 0};

    auto res1 = searcher.search_path(20);
    auto res2 = searcher.search_path(500);

    if(res1.has_value())
    {
        auto [path, distance] = res1.value();
        std::cout << "Path from 0 to 20: ";
        for(auto node : path)
        {
            std::cout << node << " ";
        }
        std::cout << std::endl;
        std::cout << "Distance: " << distance << std::endl;
    }

    if(res2.has_value())
    {
        auto [path, distance] = res2.value();
        std::cout << "Path from 0 to 500: ";
        for(auto node : path)
        {
            std::cout << node << " ";
        }
        std::cout << std::endl;
        std::cout << "Distance: " << distance << std::endl;
    }

    return 0;
}