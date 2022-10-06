#ifndef LATTICES_H
#define LATTICES_H

#include "prng.hpp"
#include <vector>
#include <algorithm>
using coord = std::pair<uint, uint>;

template<typename T>
struct square_lattice
{
    uint size = 0;
    std::vector<double> initial_probabilities = 0;
    std::vector<T> initial_states;
    T** lattice_sites;
    prng gen;

     
    square_lattice(uint L, std::vector<double> probabilities, std::vector<T> choices):
    size(L), initial_probabilities(probabilities), initial_states(choices)
    {

        if (initial_probabilities.size() != initial_states.size())
        {
            throw std::invalid_argument("Invalid number of initial_probabilities. Must be equal to number of choices.");
        }

        double sum = 0;
        for (auto p : initial_probabilities)
        {
            if (p < 0 || p > 1)
            {
                throw std::invalid_argument("Invalid probabilities. Must be in [0,1].");
            }

            sum += p;
        }
        if (std::abs(sum - 1.0) > 1e-14)
        {
            throw std::invalid_argument("Invalid probabilities. Must sum to 1.");
        }
        
        lattice_sites = new T*[size];
        for(uint i = 0; i < size; i++) lattice_sites[i] = new T[size];
        
        for(uint i = 0; i < size; i++)
        {
            for(uint j = 0; j < size; j++)
            {
                lattice_sites[i][j] = random_choice(initial_states, initial_probabilities);
            }
        }

    }

    ~square_lattice() {
        for(uint i = 0; i < size; i++) delete[] lattice_sites[i];        
        delete[] lattice_sites;
    }
    
    template<typename U>
    U random_choice(std::vector<U> states, std::vector<double> probabilities)
    {
        double r = gen.random_double();
        double sum = 0;
        for (uint i = 0; i < probabilities.size(); i++)
        {
            sum += probabilities[i];
            if (r < sum)
            {
                return states[i];
            }
        }
        return states[states.size() - 1];
    }
    
    void set_lattice_site(coord coord, T val)
    {
        lattice_sites[coord.first][coord.second] = val;
    }
    void print_lattice(){
        for(uint i = 0; i < size; i++)
        {
            for(uint j = 0; j < size; j++)
            {
                std::cout << lattice_sites[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    void save_lattice_to_file(std::string filename){
        std::ofstream output(filename);
        for(uint i = 0; i < size; i++)
        {
            for(uint j = 0; j < size; j++)
            {
                output << lattice_sites[i][j] << " ";
            }
            output << std::endl;
        }
        output.close();
    }
    void print_lattice_params(){
        std::cout << "# Lattice size is " << size << " x " << size << std::endl;
        for(uint i = 0; i < initial_probabilities.size(); i++)
        {
            std::cout << "# Probability of site being occupied by " << initial_states[i] << " is " << initial_probabilities[i] << std::endl;
        }
    }



    std::vector<coord> possible_moves_nn(coord current_position, std::vector<T> allowed_sites)
    {
        std::vector<coord> possible_moves;
        
        if (current_position.first > 0)
        {
            if (std::find(allowed_sites.begin(), allowed_sites.end(), lattice_sites[current_position.first - 1][current_position.second]) != allowed_sites.end())
            {
                possible_moves.push_back(std::make_pair(current_position.first - 1, current_position.second));
            }
        }
        if (current_position.first < size - 1)
        {
            if (std::find(allowed_sites.begin(), allowed_sites.end(), lattice_sites[current_position.first + 1][current_position.second]) != allowed_sites.end())
            {
                possible_moves.push_back(std::make_pair(current_position.first + 1, current_position.second));
            }
        }
        if (current_position.second > 0)
        {
            if (std::find(allowed_sites.begin(), allowed_sites.end(), lattice_sites[current_position.first][current_position.second - 1]) != allowed_sites.end())
            {
                possible_moves.push_back(std::make_pair(current_position.first, current_position.second - 1));
            }
        }
        if (current_position.second < size - 1)
        {
            if (std::find(allowed_sites.begin(), allowed_sites.end(), lattice_sites[current_position.first][current_position.second + 1]) != allowed_sites.end())
            {
                possible_moves.push_back(std::make_pair(current_position.first, current_position.second + 1));
            }
        }
        return possible_moves;
    }

};


#endif //LATTICES_H
