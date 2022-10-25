/*
 * @file graph.cpp
 * @author Jakub Pawłowski
 * @version 0.3 25.10.2022
 * */

#include <iostream>
#include "graph.hpp"

/*
 * Basic constructor for StaticNetwork class.
 * Initializes network from provided .metis file/
 */
graph::graph(const std::string &graph_filename)
{
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
    line.erase(line.find('#'), 1);
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
    node_degree_.resize(number_of_nodes_);
    node_index_.resize(number_of_nodes_);
    adjacency_list_.resize(2 * number_of_links_);

    /*
     * Read adjacency list from rest of .metis file
     */

    std::size_t sum_of_degrees = 0;
    std::size_t node_num = 0;
    while (getline(file, line))
    {
        node_index_[node_num] = sum_of_degrees;
        count = 0;
        while ((pos = line.find(delimiter)) != std::string::npos)
        {
            std::string token = line.substr(0, pos);
            adjacency_list_[node_index_[node_num] + count] = std::stoi(token) - 1;
            count++;
            line.erase(0, pos + delimiter.length());
        }
        node_degree_[node_num] = count;
        sum_of_degrees += count;
        node_num++;
    }

    max_degree_ = *std::max(node_degree_.begin(), node_degree_.end());
}

/*
 * Returns the number of nodes
 */
std::size_t graph::get_number_of_nodes() const
{
    return number_of_nodes_;
}

/*
 * Returns the number of links
 */
std::size_t graph::get_number_of_links() const
{
    return number_of_links_;
}

/*
 * Returns vector of neighbors of node
 */
std::vector<uint> graph::get_node_neighbors(uint node)
{
    if (node >= number_of_nodes_)
    {
        throw std::invalid_argument("graph::get_node_neighbors: Node number is out of range");
    }
    std::vector<uint> neighbors;
    uint degree = static_cast<uint>(node_degree_[node]);
    uint index = static_cast<uint>(node_index_[node]);
    neighbors.resize(degree);
    for (uint i = 0; i < degree; i++)
    {
        neighbors[i] = adjacency_list_[index + i];
    }
    return neighbors;
}

/*
 * Prints adjacency list to standard output.
 */
void graph::print_adjacency_list()
{
    for (std::size_t i = 0; i < number_of_nodes_; i++)
    {
        uint degree = node_degree_[i];
        uint index = node_index_[i];
        for (uint j = 0; j < degree; j++)
        {
            std::cout << adjacency_list_[index + j] << " ";
        }
        std::cout << std::endl;
    }
}

/*
 * Prints adjacency list to file.
 */
void graph::print_adjacency_list(const std::string &filename)
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
        uint degree = node_degree_[i];
        uint index = node_index_[i];
        for (uint j = 0; j < degree; j++)
        {
            file << adjacency_list_[index + j] << " ";
        }
        file << std::endl;
    }
}

bool graph::is_linked(const int src, const int dest)
{
    if (src >= number_of_nodes_ || dest >= number_of_nodes_)
    {
        throw std::invalid_argument("graph::is_linked: Node number is out of range");
    }
    uint src_deg = node_degree_[src];
    uint src_index = node_index_[src];
    for (uint i = 0; i < src_deg; i++)
    {
        if (adjacency_list_[src_index + i] == dest)
            return true;
    }
    return false;
}

/*
 * Returns average degree of node
 */
uint graph::get_average_degree()
{
    if (average_degree_ == 0)
        calculate_average_degree();
    return average_degree_;
}

uint graph::get_average_clustering()
{
    if (average_clustering_ == 0)
        calculate_average_clustering();
    return average_clustering_;
}

/*
 * Calculates average degree of node
 */
void graph::calculate_average_degree()
{
    uint degree = 0;
    for (uint i = 0; i < number_of_nodes_; i++)
        degree += node_degree_[i];
    degree /= number_of_nodes_;
    average_degree_ = degree;
}

// void graph::calculate_average_clustering(){

// }

/*
 * Calculates distribution of node degrees
 */
std::vector<uint> graph::get_degree_distribution()
{
    std::vector<uint> degree_distribution;
    degree_distribution.resize(max_degree_);
    for (uint i = 0; i < number_of_nodes_; i++)
    {
        degree_distribution[node_degree_[i]]++;
    }
    return degree_distribution;
}


std::string graph::get_name() const
{
    return network_name_;
}
