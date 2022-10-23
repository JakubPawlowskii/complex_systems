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
            _burned = _lattice->clone();
            _clusters = _lattice->clone();
            break;
        case lattice_type::triangular:
            _lattice = std::unique_ptr<lattice_2d<uint, mode, bc::open>>(new triangular_lattice<uint, mode, bc::open>(rows, cols, {empty, tree}, {1 - _probability, _probability}));
            _burned = _lattice->clone();
            _clusters = _lattice->clone();
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
                            int k_choice = k1 > 0 ? k1 : k2 > 0 ? k2
                                                                : k3;
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
        _cluster_sizes_distribution.resize(_lattice->get_n_sites(), 0);
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

    void naive_cluster_counting()
    {
        uint cluster = 10;

        std::function<std::size_t(std::size_t,std::size_t)> area = [&area, &cluster, this](std::size_t row, std::size_t col) -> std::size_t
        {
            if (_clusters->at({row, col}) == tree)
            {
                _clusters->at({row, col}) = cluster;
                auto nn = _clusters->get_nn({row, col});
                std::size_t area_count = 1;
                for (auto &n : nn)
                {
                    area_count += area(n.first, n.second);
                }
                return area_count;
            }
            else
                return 0;
        };

        for (std::size_t row = 0; row < _clusters->get_n_rows(); row++)
        {
            for (std::size_t col = 0; col < _clusters->get_n_cols(); col++)
            {
                std::size_t current_cluster_size = area(row, col);
                if (current_cluster_size > 0)
                {
                    _cluster_sizes.push_back(current_cluster_size);
                    cluster++;
                }
            }
        }
        _cluster_sizes_distribution.resize(_lattice->get_n_sites(), 0);
        for (int size : _cluster_sizes)
        {
            if (_max_cluster_size < size)
                _max_cluster_size = size;
            if (size > 0)
                _cluster_sizes_distribution[size]++;
        }
    }
};

int main(int argc, char **argv)
{

    if (argc != 2)
    {
        throw std::invalid_argument("Invalid number of cmd arguments. Must be 1.");
    }

    size_t L = 4;
    uint mcs = 10;
    double p0 = 0.01;
    double pk = 1.0;
    double dp = 0.1;

    std::string init_filename(argv[1]);
    std::ifstream file(init_filename);

    for (std::string line; std::getline(file, line);)
    {
        if (line.substr(0, 2) == "#L")
        {
            L = std::stoi(line.substr(2));
        }
        else if (line.substr(0, 2) == "#T")
        {
            mcs = std::stoi(line.substr(2));
        }
        else if (line.substr(0, 3) == "#p0")
        {
            p0 = std::stod(line.substr(3));
        }
        else if (line.substr(0, 3) == "#pk")
        {
            pk = std::stod(line.substr(3));
        }
        else if (line.substr(0, 3) == "#dp")
        {
            dp = std::stod(line.substr(3));
        }
        else
        {
            std::cout << line.substr(0, 2) << std::endl;
            throw std::invalid_argument("Invalid values in initialization file. Allowed values are #L, #T, #p0, #pk, #dp.");
        }
    }

    std::vector<double> p_values;
    for (double p = p0; p <= pk; p += dp)
    {
        p_values.push_back(p);
    }

    std::vector<double> cluster_sizes_distribution;
    cluster_sizes_distribution.resize(L * L, 0);

    std::vector<double> max_cluster_size(p_values.size(), 0);
    // std::vector<double> p_flow(p_values.size(),0.0);

    std::cout << "# ----------------------------------------" << std::endl;
    std::cout << "# Monte Carlo simulation of percolation on square lattice" << std::endl;
    std::cout << "# ----------------------------------------" << std::endl;
    std::cout << "# Parameters from " << init_filename << " : " << std::endl;
    std::cout << "# L = " << L << std::endl;
    std::cout << "# MCS = " << mcs << std::endl;
    std::cout << "# dp = " << dp << std::endl;
    std::cout << "# p0 = " << p0 << std::endl;
    std::cout << "# pk = " << pk << std::endl;
    std::cout << "# Starting simulation..." << std::endl;

    timer time;

    std::string filename;
    for (uint i = 0; i < p_values.size(); i++)
    {
        double p = p_values[i];
        std::cout << "# Simulating p = " << p << std::endl;
        for (uint j = 0; j < mcs; j++)
        {
            site_percolation<neighbors::calculate_on_the_fly> percolation(L, L, lattice_type::triangular, p);
            // percolation.burning_algorithm();
            // percolation.find_clusters();
            percolation.naive_cluster_counting();

            // p_flow[i] += static_cast<double>(percolation.is_percolating());
            max_cluster_size[i] += static_cast<double>(percolation.get_max_cluster_size());

            auto cluster_distribution = percolation.get_cluster_sizes_distribution();

            for (uint s = 0; s < cluster_distribution.size(); s++)
            {
                cluster_sizes_distribution[s] += static_cast<double>(cluster_distribution[s]);
            }
        }

        // filename = "naive_cluster_sizes_distribution_triangle_L" + std::to_string(L) + "_mcs" + std::to_string(mcs) + "_p" + std::to_string(p_values[i]) + ".dat";
        filename = "naive_cluster_sizes_distribution_square_L" + std::to_string(L) + "_mcs" + std::to_string(mcs) + "_p" + std::to_string(p_values[i]) + ".dat";
        std::ofstream file(filename);
        for (uint s = 0; s < cluster_sizes_distribution.size(); s++)
        {
            if (cluster_sizes_distribution[s] > 0)
                file << s << " " << cluster_sizes_distribution[s] / mcs << std::endl;
        }
        file.close();
        cluster_sizes_distribution.clear();
        cluster_sizes_distribution.resize(L * L, 0);

        max_cluster_size[i] /= mcs;
        // p_flow[i] /= mcs;
    }

    // save p_flow to file together with p
    // name the file with simulation parameters
    // filename = "p_flow_triangle_L" + std::to_string(L) + "_mcs" + std::to_string(mcs) + ".dat";
    // // filename = "p_flow_square_L" + std::to_string(L) + "_mcs" + std::to_string(mcs) + ".dat";
    // std::ofstream p_flow_file(filename);
    // for (uint i = 0; i < p_values.size(); i++)
    // {
    //     p_flow_file << p_values[i] << " " << p_flow[i] << std::endl;
    // }
    // p_flow_file.close();

    // save max cluster size to file together with p
    // name the file with simulation parameters
    // filename = "naive_max_cluster_size_triangle_L" + std::to_string(L) + "_mcs" + std::to_string(mcs) + ".dat";
    filename = "naive_max_cluster_size_square_L" + std::to_string(L) + "_mcs" + std::to_string(mcs) + ".dat";
    std::ofstream max_cluster_size_file(filename);
    for (uint i = 0; i < p_values.size(); i++)
    {
        max_cluster_size_file << p_values[i] << " " << max_cluster_size[i] << std::endl;
    }
    max_cluster_size_file.close();

    std::cout << "# Time elapsed: " << time.elapsed() << " s" << std::endl;

    return 0;
}