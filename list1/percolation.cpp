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
    std::size_t _max_cluster_size = 0;
    std::unique_ptr<lattice_2d<uint, mode, bc::open>> _burned;
    std::unique_ptr<lattice_2d<uint, mode, bc::open>> _clusters;

public:
    site_percolation() = default;
    site_percolation(std::size_t rows, std::size_t cols, lattice_type lattice, double probability) : _probability(probability)
    {
        switch (lattice)
        {
        case lattice_type::square:
            _burned =  std::unique_ptr<lattice_2d<uint, mode, bc::open>>(new square_lattice<uint, mode, bc::open>(rows, cols, {empty, tree},{1-_probability,_probability}));
            _clusters = _burned->clone();
            break;
        case lattice_type::triangular:
            _burned =  std::unique_ptr<lattice_2d<uint, mode, bc::open>>(new triangular_lattice<uint, mode, bc::open>(rows, cols, {empty, tree},{1-_probability,_probability}));
            _clusters = _burned->clone();
            break;
        }
    }
    ~site_percolation() = default;

    void print_burned(){
        _burned->print_lattice();
    }

    bool is_percolating() const
    {
        return _percolates;
    }
    std::size_t max_cluster_size() const
    {
        return _max_cluster_size;
    }
    std::size_t shortest_path_length() const
    {
        return _shortest_path_length;
    }
    void burning_algorithm(){
        uint t = 2;
        std::vector<coord> currently_burning;
        std::vector<coord> currently_burning_new;

        for (std::size_t i = 0; i < _burned->get_n_cols() ; i++)
        {
            if (_burned->at({0, i}) == tree)
            {
                _burned->at({0, i}) = t;
                currently_burning.push_back(std::make_pair(0,i));
            }
        }
        bool is_burning = true;
        while(is_burning)
        {
            is_burning = false;
            for(auto fire : currently_burning)
            {
                auto valid_nn = _burned->get_nn(fire);
                valid_nn.erase(std::remove_if(valid_nn.begin(), valid_nn.end(), [&](coord c) { return _burned->at(c) != tree; }), valid_nn.end());
                for(auto nn : valid_nn)
                {
                    if(nn.first == _burned->get_n_rows() - 1)
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


};

int main(int argc, char **argv)
{
    timer t;
    site_percolation<neighbors::calculate_on_the_fly> percolation(10,10, lattice_type::square, 0.59275);
    // percolation.print_burned();
    percolation.burning_algorithm();
    std::cout << std::endl;
    // percolation.print_burned();
    std::cout << "Percolates: " << percolation.is_percolating() << std::endl;
    std::cout << "Shortest path length: " << percolation.shortest_path_length() << std::endl;
    std::cout << "Time elapsed: " << t.elapsed() << " s" << std::endl;

    return 0;
}