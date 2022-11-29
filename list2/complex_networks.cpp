#include <iostream>

#include "../utils/graph.hpp"

int main()
{

    std::size_t num_nodes = 5000;
    // std::vector<double> probs{0.0001, 0.0005, 0.01, 0.003, 0.001, 0.01, 0.05, 0.4};
    std::vector<double> probs{0.005, 0.01, 0.1, 0.4};

    // vector of degrees for each node given by N*p
    std::vector<std::size_t> degrees(probs.size());
    for (std::size_t i = 0; i < degrees.size(); ++i)
    {
        degrees[i] = static_cast<std::size_t>(num_nodes * probs[i]/2);
    }
    std::vector<double> betas = {0.001, 0.01, 0.1, 0.5};

    // sort the probabilities and make them unique
    std::sort(probs.begin(), probs.end());
    probs.erase(std::unique(probs.begin(), probs.end()), probs.end());
    for (auto beta : betas)
    {
        for (auto k : degrees)
        {
            std::cout << "k = " << k << ", beta = " << beta << std::endl;
            // std::cout << "p = " << p << std::endl;
            auto g = graph{};
            g.watts_strogatz(num_nodes, k, beta);
            // g.erdos_renyi_simple(num_nodes, p);
            std::string graph_name = "watts_strogatz_k_" + std::to_string(k) + "_beta_" + std::to_string(beta) + ".txt";
            // std::string graph_name = "erdos_renyi_p_" + std::to_string(p) + ".txt";
            std::cout << "Number of nodes: " << g.get_number_of_nodes() << std::endl;
            std::cout << "Number of edges: " << g.get_number_of_links() << std::endl;
            std::cout << "Average degree: " << g.get_average_degree() << std::endl;
            // std::cout << "Average clustering: " << g.get_average_clustering() << std::endl;
            auto path_to_file = "/home/jakubp/Code/Studies/semester_2/complex_systems/list2/networks/";
            g.print_adjacency_list(path_to_file + graph_name);
        }
    }

    // g.from_adj_list_to_nodes_and_edges();

    // dijkstra searcher{g, 0};

    // auto res1 = searcher.search_path(20);
    // auto res2 = searcher.search_path(500);

    // if(res1.has_value())
    // {
    //     auto [path, distance] = res1.value();
    //     std::cout << "Path from 0 to 20: ";
    //     for(auto node : path)
    //     {
    //         std::cout << node << " ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "Distance: " << distance << std::endl;
    // }

    // if(res2.has_value())
    // {
    //     auto [path, distance] = res2.value();
    //     std::cout << "Path from 0 to 500: ";
    //     for(auto node : path)
    //     {
    //         std::cout << node << " ";
    //     }
    //     std::cout << std::endl;
    //     std::cout << "Distance: " << distance << std::endl;
    // }

    return 0;
}