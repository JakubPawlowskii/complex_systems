#include <iostream>

#include "../utils/graph.hpp"


int main()
{
    auto g = graph{};
    g.read_from_metis("/home/jakubp/Code/Studies/semester_2/complex_systems/list2/networks/facebook.metis");
    
    std::cout << "Number of nodes: " << g.get_number_of_nodes() << std::endl;
    std::cout << "Number of edges: " << g.get_number_of_links() << std::endl;
    std::cout << "Average degree: " << g.get_average_degree() << std::endl;
    std::cout << "Average clustering: " << g.get_average_clustering() << std::endl;

    g.from_adj_list_to_nodes_and_edges();

    std::size_t shortest_path = std::numeric_limits<std::size_t>::max();

    dijkstra searcher(g, g.vertex_data(10));
    auto res = searcher.search_path(g.vertex_data(1000));
    if(res.has_value())
    {
        auto [path,len] = res.value();
        for(auto v : path ) std::cout << v << " --- ";
        std::cout << std::endl;
        std::cout << "Path length is " << len << std::endl;

    }

    // for(std::size_t i = 0; i < g.get_number_of_nodes(); i++)
    // {
    //     std::cout << "Node " << i << " has " << g.get_node_neighbors(i).size() << " neighbors" << std::endl;
    //     std::cout << std::flush;
    //     auto src = g.vertex_data(i);
    //     dijkstra searcher{g, src};
    //     for(std::size_t j = i+1; j < g.get_number_of_nodes(); j++)
    //     {
    //         auto dest = g.vertex_data(j);
    //         auto res = searcher.search_path(dest);

    //         if(res.has_value())
    //         {
    //             auto path_len = res.value().second;

    //             if (path_len < shortest_path)
    //             {
    //                 shortest_path = path_len;
    //             }
    //         }
    //     }
    // }

    // std::cout << "Shortest path: " << shortest_path << std::endl;

   return 0;
}