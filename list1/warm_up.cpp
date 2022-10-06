#include <iostream>
#include <fstream>
#include "../utils/prng.hpp"
#include "../utils/lattices.hpp"

enum class site
{
    empty,
    dog,
    visited_by_flea,
    current_flea
};

std::ostream &operator<<(std::ostream &out, const site &site)
{
    switch (site)
    {
    case site::empty:
        // out << "\033[1;37me\033[0m";
        out << "\033[38:5:15me\033[0m";
        break;
    case site::dog:
        // out << "\033[1;31md\033[0m";
        out << "\033[38:5:160md\033[0m";
        break;
    case site::visited_by_flea:
        // out << "\033[1;32mf\033[0m";
        out << "\033[38:5:21mf\033[0m";
        break;
    case site::current_flea:
        // out << "\033[1;32mf\033[0m";
        out << "\033[38:5:82mf\033[0m";
        break;
    }

    return out;
}

int main(int argc, char **argv)
{

    if (argc != 2)
    {
        throw std::invalid_argument("Invalid number of cmd arguments. Must be 1.");
    }

    uint L = 4;
    double p = 0.1;

    std::string init_filename(argv[1]);
    std::ifstream file(init_filename);

    for (std::string line; std::getline(file, line);)
    {
        if (line.substr(0, 2) == "#L")
        {
            L = std::stoi(line.substr(2));
        }
        else if (line.substr(0, 2) == "#p")
        {
            p = std::stod(line.substr(2));
        }
        else
        {
            throw std::invalid_argument("Invalid values in initialization file. Must be only L and p.");
        }
    }

    // Task 1
    // square_lattice<site> lattice(L, {p,1-p}, {site::dog,site::empty});
    // lattice.print_lattice_params();
    // lattice.print_lattice();

    // Task 2

    square_lattice<site> lattice(L, {p, 1 - p}, {site::dog, site::empty});
    auto find_first = [&](site type) -> coord
    {
        for (uint i = 0; i < L; i++)
        {
            for (uint j = 0; j < L; j++)
            {
                if (lattice.lattice_sites[i][j] == type)
                {
                    return {i, j};
                }
            }
        }
        return {-1, -1};
    };

    lattice.print_lattice();
    std::cout << std::endl;
    coord current_pos = find_first(site::dog);
    lattice.set_lattice_site(current_pos, site::current_flea);
    lattice.print_lattice();
    std::cout << std::endl;

    auto possible_moves = lattice.possible_moves_nn(current_pos, {site::dog, site::visited_by_flea});
    if (possible_moves.empty())
    {
        std::cout << "Cannot move." << std::endl;
    }
    else
    {
        std::vector<double> probabilities;
        probabilities.resize(possible_moves.size());
        for (uint i = 0; i < possible_moves.size(); i++)
        {
            probabilities[i] = 1.0 / possible_moves.size();
        }

        coord next_pos = lattice.random_choice<coord>(possible_moves, probabilities);
        lattice.set_lattice_site(current_pos, site::visited_by_flea);
        lattice.set_lattice_site(next_pos, site::current_flea);
        current_pos = next_pos;
        lattice.print_lattice();
    }

    return 0;
}