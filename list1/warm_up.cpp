#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include "../utils/prng.hpp"
#include "../utils/lattices2d.hpp"

enum class site
{
    empty,
    dog,
    visited_by_flea,
    flea
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
        out << "\033[38:5:39mf\033[0m";
        break;
    case site::flea:
        // out << "\033[1;32mf\033[0m";
        out << "\033[1m\033[38:5:82mf\033[0m";
        break;
    }

    return out;
}

int main(int argc, char **argv)
{
    using namespace std::this_thread;
    using namespace std::chrono_literals;
    using std::chrono::system_clock;

    using namespace lattices2d;

    if (argc != 2)
    {
        throw std::invalid_argument("Invalid number of cmd arguments. Must be 1.");
    }

    size_t L = 4;
    double p = 0.1;
    uint max_t = 100;

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
        else if(line.substr(0,2) == "#t")
        {
            max_t = std::stoi(line.substr(2));
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

    square_lattice<site, neighbors_calculation::calculate_on_the_fly, bc::open> lattice(L, L, {site::dog, site::empty}, {p, 1-p});
    for (auto neighbor : lattice.get_nn({0,0}))
    {
        std::cout << neighbor << std::endl;
    }

    auto find_first = [&](site type) -> coord
    {
        for (uint i = 0; i < L; i++)
        {
            for (uint j = 0; j < L; j++)
            {
                if (lattice({i,j}) == type)
                {
                    return {i, j};
                }
            }
        }
        return {-1, -1};
    };

    
    std::cout << "Probability of dog: " << p << std::endl;
    lattice.print_lattice();
    std::cout << "Press any key to start simulation." << std::endl;
    std::cin.get();

    coord current_flea_pos = find_first(site::dog);
    lattice(current_flea_pos) = site::flea;
    bool flag = true;
    for(uint t = 0; t <= max_t; t++)
    {
        std::cout << "\033[2J\033[;H"; // clearing the screen
        std::cout << "t = " << t << std::endl;
        lattice.print_lattice();
        sleep_for(0.05s);
        auto possible_moves = lattice.get_nn(current_flea_pos);
        std::remove_if(possible_moves.begin(), possible_moves.end(), [&](coord c) { return lattice(c) == site::empty; });

        if (possible_moves.empty())
        {
            std::cout << "Flea cannot move. Ending simulation." << std::endl;
            flag = false;
            break;
        }
        else
        {
            std::vector<double> probabilities;
            probabilities.resize(possible_moves.size());
            for (uint i = 0; i < possible_moves.size(); i++)
            {
                probabilities[i] = 1.0 / possible_moves.size();
            }

            lattice(current_flea_pos) = site::visited_by_flea;
            current_flea_pos = lattice.random_choice<coord>(possible_moves, probabilities);
            lattice(current_flea_pos) =  site::flea;
        }

    }
    if(flag) std::cout << "Max time reached. Ending simulation." << std::endl;
    return 0;
}