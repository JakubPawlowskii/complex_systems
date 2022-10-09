#include <iostream>
#include <fstream>
#include <memory>
#include "../utils/prng.hpp"
#include "../utils/lattices2d.hpp"
#include "../utils/timer.hpp"

using namespace lattices2d;
constexpr uint empty = 0;
constexpr uint tree = 1;

template <neighbors mode>
class site_percolation
{
private:
    double _probability;
    bool _percolates = false;
    uint _shortest_path_length = 0;
    int _max_cluster_size = 0;
    lattice_type _lattice_type;
    std::vector<int> _cluster_sizes;
    std::vector<int> _cluster_sizes_distribution;
    std::unique_ptr<lattice_2d<uint, mode, bc::open>> _lattice;
    std::unique_ptr<lattice_2d<uint, mode, bc::open>> _burned;
    std::unique_ptr<lattice_2d<uint, mode, bc::open>> _clusters;

public:
    site_percolation() = default;
    site_percolation(std::size_t rows, std::size_t cols, lattice_type lattice, double probability) : _probability(probability), _lattice_type(lattice)
    {
        switch (_lattice_type)
        {
        case lattice_type::square:
            _lattice = std::unique_ptr<lattice_2d<uint, mode, bc::open>>(new square_lattice<uint, mode, bc::open>(rows, cols, {empty, tree}, {1 - _probability, _probability}));
            break;
        case lattice_type::triangular:
            _lattice = std::unique_ptr<lattice_2d<uint, mode, bc::open>>(new triangular_lattice<uint, mode, bc::open>(rows, cols, {empty, tree}, {1 - _probability, _probability}));
            break;
        }
        _cluster_sizes = {0, 0};
    }
    ~site_percolation() = default;

    void print_lattice()
    {
        _lattice->print_lattice();
    }
    void print_burned()
    {
        _burned->print_lattice();
    }
    void print_clusters()
    {
        _clusters->print_lattice();
    }

    bool is_percolating() const
    {
        return _percolates;
    }
    std::size_t get_max_cluster_size() const
    {
        return _max_cluster_size;
    }
    std::size_t get_shortest_path_length() const
    {
        return _shortest_path_length;
    }
    std::vector<int> get_cluster_sizes() const
    {
        return _cluster_sizes;
    }
    std::vector<int> get_cluster_sizes_distribution() const
    {
        return _cluster_sizes_distribution;
    }

    void burning_algorithm()
    {
        uint t = 2;
        _burned = _lattice->clone();
        std::vector<coord> currently_burning;
        std::vector<coord> currently_burning_new;

        for (std::size_t i = 0; i < _burned->get_n_cols(); i++)
        {
            if (_burned->at({0, i}) == tree)
            {
                _burned->at({0, i}) = t;
                currently_burning.push_back(std::make_pair(0, i));
            }
        }
        bool is_burning = true;
        while (is_burning)
        {
            is_burning = false;
            for (auto fire : currently_burning)
            {
                auto valid_nn = _burned->get_nn(fire);
                valid_nn.erase(std::remove_if(valid_nn.begin(), valid_nn.end(), [&](coord c)
                                              { return _burned->at(c) != tree; }),
                               valid_nn.end());
                for (auto nn : valid_nn)
                {
                    if (nn.first == _burned->get_n_rows() - 1)
                    {
                        _percolates = true;
                        _shortest_path_length = t;
                        return;
                    }
                    _burned->at(nn) = t + 1;
                    is_burning = true;
                    currently_burning_new.push_back(nn);
                }
            }
            currently_burning = currently_burning_new;
            currently_burning_new.clear();
            t++;
        }
        _shortest_path_length = t - 1;
        return;
    }
    void find_clusters(bool show_true_clusters = false)
    {
        int k = 2;
        _clusters = _lattice->clone();

        std::function<int(int, const std::vector<int> &)> detect_cluster_number =
            [&detect_cluster_number](int k, const std::vector<int> &cluster_sizes)
        {
            if (cluster_sizes[k] >= 0)
                return k;
            else
                return detect_cluster_number(std::abs(cluster_sizes[k]), cluster_sizes);
        };

        for (std::size_t row = 0; row < _clusters->get_n_rows(); row++)
        {
            for (std::size_t col = 0; col < _clusters->get_n_cols(); col++)
            {
                // std::cout << std::endl;
                // _clusters->print_lattice();

                if (_clusters->at({row, col}) == tree)
                {
                    auto valid_nn = _lattice->get_nn({row, col});

                    valid_nn.erase(std::remove_if(valid_nn.begin(), valid_nn.end(), [&row, &col, this](coord c)
                                                  { 
                                                    bool flag = _lattice->at(c) != tree;
                                                    if(c.first > row) return true;
                                                    if(c.second > col) return true;
                                                    // if(this->_lattice_type == lattice_type::triangular && (row % 2 == 0) && c.first < row && c.second == col) return true;
                                                    return flag; }),
                                   valid_nn.end());
                    int k0, k1, k2, k3 = 0;
                    switch (valid_nn.size())
                    {
                    case 0:
                        _cluster_sizes.push_back(1);
                        _clusters->at({row, col}) = k;
                        k++;
                        break;
                    case 1:
                        k0 = detect_cluster_number(_clusters->at(valid_nn[0]), _cluster_sizes);
                        _cluster_sizes[k0]++;
                        _clusters->at({row, col}) = k0;
                        break;
                    case 2:
                        k1 = detect_cluster_number(_clusters->at(valid_nn[0]), _cluster_sizes);
                        k2 = detect_cluster_number(_clusters->at(valid_nn[1]), _cluster_sizes);
                        if (k1 == 0 || k2 == 0)
                        {
                            int k_choice = k1 > 0 ? k1 : k2;
                            _cluster_sizes[k_choice]++;
                            _clusters->at({row, col}) = k_choice;
                        }
                        else if (k1 == k2)
                        {
                            _cluster_sizes[k1]++;
                            _clusters->at({row, col}) = k1;
                        }
                        else if (k1 < k2)
                        {
                            _cluster_sizes[k1] += _cluster_sizes[k2] + 1;
                            _cluster_sizes[k2] = -k1;
                            _clusters->at({row, col}) = k1;
                        }
                        else
                        {
                            _cluster_sizes[k2] += _cluster_sizes[k1] + 1;
                            _cluster_sizes[k1] = -k2;
                            _clusters->at({row, col}) = k2;
                        }
                        break;
                    case 3:
                        k1 = detect_cluster_number(_clusters->at(valid_nn[0]), _cluster_sizes);
                        k2 = detect_cluster_number(_clusters->at(valid_nn[1]), _cluster_sizes);
                        k3 = detect_cluster_number(_clusters->at(valid_nn[2]), _cluster_sizes);
                        if (k1 == 0 || k2 == 0 || k3 == 0)
                        {
                            int k_choice = k1 > 0 ? k1 : k2 > 0 ? k2 : k3;
                            _cluster_sizes[k_choice]++;
                            _clusters->at({row, col}) = k_choice;
                        }
                        else if (k1 == k2 && k2 == k3)
                        {
                            _cluster_sizes[k1]++;
                            _clusters->at({row, col}) = k1;
                        }
                        else if (k1 == k2)
                        {
                            _cluster_sizes[k1] += _cluster_sizes[k3] + 1;
                            _cluster_sizes[k3] = -k1;
                            _clusters->at({row, col}) = k1;
                        }
                        else if (k1 == k3)
                        {
                            _cluster_sizes[k1] += _cluster_sizes[k2] + 1;
                            _cluster_sizes[k2] = -k1;
                            _clusters->at({row, col}) = k1;
                        }
                        else if (k2 == k3)
                        {
                            _cluster_sizes[k2] += _cluster_sizes[k1] + 1;
                            _cluster_sizes[k1] = -k2;
                            _clusters->at({row, col}) = k2;
                        }
                        else if (k1 < k2 && k1 < k3)
                        {
                            _cluster_sizes[k1] += _cluster_sizes[k2] + _cluster_sizes[k3] + 2;
                            _cluster_sizes[k2] = -k1;
                            _cluster_sizes[k3] = -k1;
                            _clusters->at({row, col}) = k1;
                        }
                        else if (k2 < k1 && k2 < k3)
                        {
                            _cluster_sizes[k2] += _cluster_sizes[k1] + _cluster_sizes[k3] + 2;
                            _cluster_sizes[k1] = -k2;
                            _cluster_sizes[k3] = -k2;
                            _clusters->at({row, col}) = k2;
                        }
                        else
                        {
                            _cluster_sizes[k3] += _cluster_sizes[k1] + _cluster_sizes[k2] + 2;
                            _cluster_sizes[k1] = -k3;
                            _cluster_sizes[k2] = -k3;
                            _clusters->at({row, col}) = k3;
                        }
                        break;
                    default:
                        throw std::runtime_error("Invalid number of neighbors.");
                        break;
                    }
                }
                else
                {
                    continue;
                }
            }
        }
        _cluster_sizes_distribution.reserve(_cluster_sizes.size());
        for (int size : _cluster_sizes)
        {
            if (_max_cluster_size < size)
                _max_cluster_size = size;
            if (size > 0)
                _cluster_sizes_distribution[size]++;
        }
        if (show_true_clusters)
        {
            for (std::size_t row = 0; row < _clusters->get_n_rows(); row++)
            {
                for (std::size_t col = 0; col < _clusters->get_n_cols(); col++)
                {
                    int k = detect_cluster_number(_clusters->at({row, col}), _cluster_sizes);
                    _clusters->at({row, col}) = k;
                }
            }
        }
    }
};

int main(int argc, char **argv)
{

    timer t;
    site_percolation<neighbors::calculate_on_the_fly> percolation(50, 50, lattice_type::square, 0.99);
    // percolation.print_burned();
    // percolation.burning_algorithm();
    // std::cout << std::endl;
    // percolation.print_burned();
    // std::cout << "Percolates: " << percolation.is_percolating() << std::endl;
    // std::cout << "Shortest path length: " << percolation.shortest_path_length() << std::endl;

    // percolation.print_lattice();
    percolation.find_clusters(true);
    std::cout << std::endl;
    percolation.print_clusters();
    auto cl = percolation.get_cluster_sizes();
    for (uint i = 0; i < cl.size(); i++)
    {
        std::cout << i << " " << cl[i] << std::endl;
    }

    std::cout << "Time elapsed: " << t.elapsed() << " s" << std::endl;

    return 0;
}