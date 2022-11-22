// Copyright (C) 2020 Martijn Saelens
//
// SPDX-License-Identifier: MIT

#ifndef DIJKSTRA_HPP
#define DIJKSTRA_HPP

#include "graph_base.hpp"

#include <cassert>
#include <limits>
#include <optional>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

template <class V = uint, class E = uint>
class dijkstra
{
private:
    constexpr static std::ptrdiff_t no_parent = -1;
    constexpr static E infinity = std::numeric_limits<E>::max();

public:
    // You're only sure you have the shortest path if you processed the whole
    // graph_base. Therefore, the Dijsktra algorithm is executed in the constructor.
    //
    dijkstra(const graph_base<V, E>& original_graph_base, std::ptrdiff_t start_vertex_id)
    : m_original_graph_base{&original_graph_base}, m_start_vertex_id{start_vertex_id}
    {
        assert(m_original_graph_base->is_vertex_present(m_start_vertex_id));

        std::unordered_set<std::ptrdiff_t> determined_vertex_ids;

        for (const auto& vertex_id : m_original_graph_base->vertex_ids())
        {
            m_vertex_parent_id[vertex_id] = no_parent;
            m_vertex_distance[vertex_id] = infinity;
        }
        m_vertex_distance[m_start_vertex_id] = E{0};

        auto cmp = [this](std::ptrdiff_t lhs, std::ptrdiff_t rhs)
        {
            return (m_vertex_distance.at(lhs) > m_vertex_distance.at(rhs));
        };
        std::priority_queue<std::ptrdiff_t, std::vector<std::ptrdiff_t>, decltype(cmp)> vertex_ids_to_process(cmp);
        vertex_ids_to_process.push(m_start_vertex_id);

        while (!vertex_ids_to_process.empty())
        {
            std::ptrdiff_t vertex_id = vertex_ids_to_process.top();
            vertex_ids_to_process.pop();

            if (!determined_vertex_ids.contains(vertex_id))
            {
                determined_vertex_ids.insert(vertex_id);

                for (const auto& edge_id : m_original_graph_base->edge_ids(vertex_id))
                {
                    std::ptrdiff_t neighbor_vertex_id = m_original_graph_base->edge_target_vertex_id(edge_id);

                    E new_distance = (m_vertex_distance[vertex_id] + m_original_graph_base->edge_data(edge_id));
                    if (new_distance < m_vertex_distance[neighbor_vertex_id])
                    {
                        m_vertex_distance[neighbor_vertex_id] = new_distance;
                        m_vertex_parent_id[neighbor_vertex_id] = vertex_id;
                        vertex_ids_to_process.push(neighbor_vertex_id);
                    }
                }
            }
        }
    }

    std::optional<std::pair<std::vector<std::ptrdiff_t>, std::ptrdiff_t>>
    search_path(std::ptrdiff_t target_vertex_id) const
    {
        assert(m_original_graph_base->is_vertex_present(target_vertex_id));

        if (m_start_vertex_id == target_vertex_id)
        {
            return std::optional<std::pair<std::vector<std::ptrdiff_t>, std::ptrdiff_t>>{{{m_start_vertex_id}, 0}};
        }

        if (m_vertex_distance.at(target_vertex_id) == infinity)
        {
            return std::nullopt;
        }

        return std::optional<std::pair<std::vector<std::ptrdiff_t>, std::ptrdiff_t>>{
                {construct_path(target_vertex_id), m_vertex_distance.at(target_vertex_id)}};
    }

private:
    std::vector<std::ptrdiff_t> construct_path(std::ptrdiff_t target_vertex_id) const
    {
        std::vector<std::ptrdiff_t> path_vertex_ids;

        std::ptrdiff_t vertex_id = target_vertex_id;
        while (vertex_id != m_start_vertex_id)
        {
            path_vertex_ids.push_back(vertex_id);
            vertex_id = m_vertex_parent_id.at(vertex_id);
        }

        path_vertex_ids.push_back(m_start_vertex_id);

        std::reverse(path_vertex_ids.begin(), path_vertex_ids.end());

        return path_vertex_ids;
    }

    const graph_base<V, E>* m_original_graph_base;

    std::ptrdiff_t m_start_vertex_id;
    std::unordered_map<std::ptrdiff_t, std::ptrdiff_t> m_vertex_parent_id;
    std::unordered_map<std::ptrdiff_t, E> m_vertex_distance;
};

#endif