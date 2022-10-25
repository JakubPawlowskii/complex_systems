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


/**
 * Stores network in form of adjacency list.
 * Generates graph from .METIS file.
 */
class graph {
private:
    std::string network_name_;
    std::size_t number_of_nodes_ = 0;      // number of nodes in graph
    std::size_t number_of_links_ = 0;      // number of links in graph
    uint average_degree_ = 0;       // average degree of a node, calculated from node_degree array
    uint average_clustering_ = 0;   // average clustering coefficient of a node, calculated from node_degree array
    std::size_t max_degree_ = 0;
    std::vector<uint> node_degree_;      // array of node degrees
    std::vector<uint> node_index_;       // array of indexes of first link of each node in 1-D adjacency_list
    /**
     * 1-D adjacency list
     * unrolled 2-D adjacency list, where all lists of neighbors are written in one row
     * accessed through get_node_neighbors
     */
    std::vector<uint> adjacency_list_;

public:
    graph(const std::string &graph_filename);

    ~graph() = default;

    [[nodiscard]] std::size_t get_number_of_nodes() const;

    [[nodiscard]] std::size_t get_number_of_links() const;

    [[nodiscard]] uint get_average_degree();
    [[nodiscard]] uint get_average_clustering();

    [[nodiscard]] std::vector<uint> get_degree_distribution();

    [[nodiscard]] std::vector<uint> get_node_neighbors(uint node);

    void print_adjacency_list();
    void print_adjacency_list(const std::string &filename);

    bool is_linked(int src, int dest);

    [[nodiscard]] std::string get_name() const;

private:
    void calculate_average_degree();
    void calculate_average_clustering();

};


#endif //GRAPH_H
