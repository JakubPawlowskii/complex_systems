/*
 * @file graph.h
 * @author Jakub Pawłowski
 * @version 0.3 25.10.2022
 * */

#ifndef GRAPH_H
#define GRAPH_H

#include <random>
#include <fstream>
#include <vector>
#include <cassert>
#include <type_traits>

#include "../extern/graph_base.hpp"
#include "../extern/dijkstra.hpp"
#include "prng.hpp"

enum class directed
{
    yes,
    no
};

/**
 * Stores network in form of adjacency list.
 * Generates graph from .METIS file.
 */

template <directed D = directed::no, class V = uint, class E = uint>
class graph : public graph_base<V, E>
{
private:
    std::string network_name_;
    std::size_t number_of_nodes_;
    std::size_t number_of_links_;
    double average_degree_ = -1.0;     // average degree of a node, calculated from node_degree array
    double average_clustering_ = -1.0; // average clustering coefficient of a node, calculated from node_degree array
    std::size_t max_degree_ = 0;
    /**
     * 2-D adjacency list
     */
    std::vector<std::vector<V>> adjacency_list_{};

public:
    graph() = default;
    ~graph() = default;

    /*
     * Initializes network from provided .metis file/
     */
    template <class U = V, class T = E>
    void read_from_metis(const std::string &graph_filename, std::enable_if_t<std::is_integral_v<U> && std::is_integral_v<T>, void> * = nullptr)
    {
        this->clear_edges();
        this->clear();

        std::string line;
        std::ifstream file;
        std::string delimiter = " ";
        std::size_t pos;

        file.open(graph_filename);
        if (!file.is_open())
        {
            perror("Error opening file with graph");
            exit(EXIT_FAILURE);
        }

        /*
         * Read number of nodes and number of links from first line of .metis file.
         * Initialize class variables.
         */
        getline(file, line);
        line.erase(line.find('%'), 1);
        network_name_ = line;

        getline(file, line);
        std::string network_parameters[2];
        int count = 0;
        while ((pos = line.find(delimiter)) != std::string::npos)
        {
            network_parameters[count++] = line.substr(0, pos);
            line.erase(0, pos + delimiter.length());
        }

        number_of_nodes_ = std::stoi(network_parameters[0]);
        number_of_links_ = std::stoi(network_parameters[1]);

        if constexpr (D == directed::yes)
            number_of_links_ *= 2;
        adjacency_list_.resize(number_of_nodes_);

        /*
         * Read adjacency list from rest of .metis file
         */

        // std::size_t sum_of_degrees = 0;
        std::size_t node_num = 0;
        max_degree_ = 0;
        std::size_t current_degree = 0;
        while (getline(file, line))
        {
            // std::cout << line << std::endl;
            adjacency_list_[node_num] = std::vector<V>{};
            while ((pos = line.find(delimiter)) != std::string::npos)
            {
                std::string token = line.substr(0, pos);
                adjacency_list_[node_num].push_back(std::stoi(token) - 1);
                line.erase(0, pos + delimiter.length());
                current_degree++;
            }

            std::string token = line.substr(0, pos);
            try
            {
                adjacency_list_[node_num].push_back(std::stoi(token) - 1);
                line.erase(0, pos + delimiter.length());
                current_degree++;
            }
            catch (const std::exception &e)
            {
            }

            node_num++;
            if (current_degree > max_degree_)
                max_degree_ = current_degree;
            current_degree = 0;
        }
    }

    // G(N, p) - Erdős-Rényi random graph
    template <class U = V, class T = E>
    void erdos_renyi_simple(const std::size_t &number_of_nodes, const double &p_edge, std::enable_if_t<std::is_integral_v<U> && std::is_integral_v<T>, void> * = nullptr)
    {
        this->clear_edges();
        this->clear();

        number_of_nodes_ = number_of_nodes;
        number_of_links_ = 0;
        adjacency_list_.resize(number_of_nodes_);

        prng rng;
        if constexpr (D == directed::yes)
        {
            for (std::size_t i = 0; i < number_of_nodes_; i++)
            {
                for (std::size_t j = 0; j < number_of_nodes_; j++)
                {
                    if (i == j)
                        continue;
                    if (rng.random_double() < p_edge)
                    {
                        adjacency_list_[i].push_back(j);
                        number_of_links_++;
                    }
                }
            }
        }
        else
        {
            for (std::size_t i = 0; i < number_of_nodes_; i++)
            {
                for (std::size_t j = i + 1; j < number_of_nodes_; j++)
                {
                    if (rng.random_double() < p_edge)
                    {
                        adjacency_list_[i].push_back(j);
                        adjacency_list_[j].push_back(i);
                        number_of_links_ += 2;
                    }
                }
            }
        }
    }

    // k-regular graph (ring lattice)
    template <class U = V, class T = E>
    void k_regular_graph(std::size_t num_nodes, std::size_t k, std::enable_if_t<std::is_integral_v<U> && std::is_integral_v<T>, void> * = nullptr)
    {
        this->clear_edges();
        this->clear();

        if (k > num_nodes - 1)
        {
            throw std::invalid_argument("k must be <= to n-1");
        }

        number_of_nodes_ = num_nodes;
        number_of_links_ = 0;
        adjacency_list_.resize(number_of_nodes_);

        if constexpr (D == directed::yes)
        {
            if (k % 2)
            {
                for (std::size_t i = 0; i < number_of_nodes_; i++)
                {
                    for (std::size_t j = 0; j < k / 2; j++)
                    {
                        adjacency_list_[i].push_back((i + j + 1) % number_of_nodes_);
                        number_of_links_++;
                    }
                    if (i < number_of_nodes_ / 2)
                    {
                        adjacency_list_[i].push_back(i + number_of_nodes_/2);
                        number_of_links_++;
                    }
                }
            }
            else
            {
                for (std::size_t i = 0; i < number_of_nodes_; i++)
                {
                    for (std::size_t j = 0; j < k / 2; j++)
                    {
                        adjacency_list_[i].push_back((i + j + 1) % number_of_nodes_);

                        number_of_links_++;
                    }
                }
            }
        }
        else
        {
            if (k % 2)
            {
                for (std::size_t i = 0; i < number_of_nodes_; i++)
                {
                    for (std::size_t j = 0; j < k / 2; j++)
                    {
                        adjacency_list_[i].push_back((i + j+1) % number_of_nodes_);
                        adjacency_list_[(i + j + 1) % number_of_nodes_].push_back(i);
                        number_of_links_++;
                    }
                    if (i < number_of_nodes_ / 2)
                    {
                        adjacency_list_[i].push_back(i + number_of_nodes_/2);
                        adjacency_list_[i + number_of_nodes_/2].push_back(i);
                        number_of_links_++;
                    }
                }
            }
            else
            {
                for (std::size_t i = 0; i < number_of_nodes_; i++)
                {
                    for (std::size_t j = 0; j < k / 2; j++)
                    {
                        adjacency_list_[i].push_back((i + j+1) % number_of_nodes_);
                        adjacency_list_[(i + j+1) % number_of_nodes_].push_back(i);

                        number_of_links_ += 2;
                    }
                }
            }
        }
    }

    // Watts-Strogatz small-world graph
    template <class U = V, class T = E>
    void watts_strogatz(std::size_t num_edges, std::size_t k, double beta, std::enable_if_t<std::is_integral_v<U> && std::is_integral_v<T>, void> * = nullptr)
    {
        k_regular_graph(num_edges, k);

        // rewire edges with probability beta
        prng rng;
        for (std::size_t i = 0; i < number_of_nodes_; i++)
        {
            for (std::size_t j = 0; j < adjacency_list_[i].size(); j++)
            {
                if (rng.random_double() < beta)
                {
                    std::size_t new_target = rng.random_uint_64() % number_of_nodes_;
                    while (new_target == i || new_target == adjacency_list_[i][j])
                    {
                        new_target = rng.random_uint_64() % number_of_nodes_;
                    }
                    adjacency_list_[i][j] = new_target;
                }
            }
        }


    }
    /*
     * Returns the number of nodes
     */
    [[nodiscard]] std::size_t get_number_of_nodes() const
    {
        return number_of_nodes_;
    }

    /*
     * Returns the number of links
     */
    [[nodiscard]] std::size_t get_number_of_links()
    {
        return number_of_links_;
    }

    /*
     * Returns average degree of nodes in graph
     */
    double get_average_degree()
    {
        if (average_degree_ < 0)
            calculate_average_degree();
        return average_degree_;
    }
    /*
     * Returns average clustering coefficient of nodes in graph
     */
    double get_average_clustering()
    {
        if (average_clustering_ < 0.0)
            calculate_average_clustering();
        return average_clustering_;
    }

    /*
     * Return distribution of node degree
     */
    std::vector<uint> get_degree_distribution()
    {
        std::vector<uint> degree_distribution;
        degree_distribution.resize(max_degree_);
        for (uint i = 0; i < number_of_nodes_; i++)
        {
            degree_distribution[adjacency_list_[i].size()]++;
        }
        return degree_distribution;
    }

    /*
     * Returns vector of neighbors of node
     */
    std::vector<V> get_node_neighbors(std::size_t node_idx)
    {
        if (node_idx >= number_of_nodes_)
        {
            throw std::invalid_argument("graph::get_node_neighbors: Node number is out of range");
        }
        return adjacency_list_[node_idx];
    }

    /*
     * Prints adjacency list to standard output.
     */
    void print_adjacency_list()
    {
        for (std::size_t i = 0; i < number_of_nodes_; i++)
        {
            for (auto v : adjacency_list_[i])
            {
                std::cout << v << " ";
            }
            std::cout << std::endl;
        }
    }

    /*
     * Prints adjacency list to file.
     */
    void print_adjacency_list(const std::string &filename)
    {
        std::ofstream file;
        file.open(filename);
        if (!file.is_open())
        {
            perror("Error writing adjacency list to file.");
            exit(EXIT_FAILURE);
        }

        for (std::size_t i = 0; i < number_of_nodes_; i++)
        {
            file << i << " ";
            for (auto v : adjacency_list_[i])
            {
                file << v << " ";
            }
            file << std::endl;
        }
        file.close();
    }

    bool is_linked(const V src, const V dest)
    {
        if (src >= number_of_nodes_ || dest >= number_of_nodes_)
        {
            throw std::invalid_argument("graph::is_linked: Node number is out of range");
        }

        for (auto v : adjacency_list_[src])
        {
            if (v == dest)
                return true;
        }
        return false;
    }
    std::string get_name() const
    {
        return network_name_;
    }

    template <class U = V, class T = E>
    void from_adj_list_to_nodes_and_edges(std::enable_if_t<std::is_integral_v<U> && std::is_integral_v<T>, void> * = nullptr)
    {
        for (std::size_t i = 0; i < number_of_nodes_; i++)
        {
            this->add_vertex(i);
        }
        for (std::size_t i = 0; i < number_of_nodes_; i++)
        {
            auto vert = this->vertex_data(i);
            if constexpr (D == directed::yes)
            {
                for (auto v : adjacency_list_[i])
                {
                    this->add_edge(vert, v, 1);
                }
            }
            else
            {
                for (auto v : adjacency_list_[i])
                {
                    this->add_undirected_edge(vert, v, 1);
                }
            }
        }

        adjacency_list_.clear();
    }

    std::size_t get_node_degree(std::size_t node_idx)
    {
        return adjacency_list_[node_idx].size();
    }

    double get_node_clustering_coeff(std::size_t node_idx)
    {
        auto neighbors = adjacency_list_[node_idx];
        double link_between_neighbors = 0;
        for (auto src : neighbors)
        {
            for (auto dsc : neighbors)
            {
                if (is_linked(src, dsc))
                    link_between_neighbors += 1.0;
            }
        }
        if constexpr (D == directed::no)
            link_between_neighbors /= 2;

        std::size_t degree = get_node_degree(node_idx);
        if (degree <= 1)
            return 0;
        else
        {
            auto den = get_node_degree(node_idx) * (get_node_degree(node_idx) - 1);
            return 2.0 * link_between_neighbors / den;
        }
    }

private:
    /*
     * Calculates average degree of node
     */
    void calculate_average_degree()
    {
        double degree = 0;
        for (std::size_t i = 0; i < number_of_nodes_; i++)
            degree += get_node_degree(i);
        degree /= number_of_nodes_;
        average_degree_ = degree;
    }

    void calculate_average_clustering()
    {
        double clustering = 0;
        for (std::size_t i = 0; i < number_of_nodes_; i++)
            clustering += get_node_clustering_coeff(i);
        clustering /= number_of_nodes_;
        average_clustering_ = clustering;
    }
};

#endif // GRAPH_H
