#include <iostream>

// #include "../extern/graph_base.hpp"

#include "../utils/graph.hpp"

int main()
{
    auto g = graph<int,int>{};
    g.add_vertex(0);
    g.add_vertex(1);
    g.add_vertex(2);
    g.add_vertex(3);
    g.add_vertex(4);

    g.add_undirected_edge(0, 1,1);
    g.add_undirected_edge(0, 2,2);
    g.add_undirected_edge(0, 3,3);
    g.add_undirected_edge(1, 2,4);
    g.add_undirected_edge(1, 3,5);
    g.add_undirected_edge(2, 3,6);
    g.add_undirected_edge(2, 4,7);
    g.add_undirected_edge(3, 4,8);

    // auto res = graph.vertex_ids();
    // for (auto i : res)
    //     std::cout << i << std::endl;

    // auto res2 = graph.edge_ids();
    // for (auto i : res2)
    //     std::cout << i << std::endl;

    auto res3 = g.vertex_data(0);
    std::cout << res3 << std::endl;

    g.clear();

   return 0;
}