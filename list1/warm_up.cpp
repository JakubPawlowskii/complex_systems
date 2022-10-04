#include <iostream>
#include <fstream>
#include "../utils/prng.hpp"

enum class site{
    empty, dog, flea
};

std::ostream& operator<< (std::ostream& out, const site& site)
{
    switch (site)
    {
    case site::empty:
        // out << "e";
        out << "\033[1;37me\033[0m";
        break;
    case site::dog:
        out << "\033[1;31md\033[0m";
        break;
    default:
        out << "\033[1;32mf\033[0m";
    }

    return out;
}

struct lattice
{
    uint size = 0;
    double filling_factor = 0;
    site** lattice_sites;
    prng gen;

    lattice(uint L, double p):
    size(L), filling_factor(p)
    {
        if(filling_factor > 1 || filling_factor < 0)
        {
            throw std::invalid_argument("Probability must be in [0,1].");
        }
        lattice_sites = new site*[size];
        for(uint i = 0; i < size; i++) lattice_sites[i] = new site[size];
        
        for(uint i = 0; i < size; i++)
        {
            for(uint j = 0; j < size; j++)
            {
                lattice_sites[i][j] = (gen.random_double() < filling_factor) ? site::dog : site::empty;
            }
        }

    }

    ~lattice() {
        for(uint i = 0; i < size; i++) delete[] lattice_sites[i];        
        delete[] lattice_sites;
    }

    void set_lattice_site(std::pair<uint,uint> coord, site val)
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
        std::cout << "# Probability of site being occupied by dog is " << filling_factor << std::endl;
    }

    std::pair<uint,uint> random_possible_jump(std::pair<uint,uint> current_coord)
    {
        
    }

};


int main(int argc, char** argv)
{

    if(argc != 2)
    {
        throw std::invalid_argument("Invalid number of cmd arguments. Must be 1.");
    }

    uint L = 4;
    double p = 0.1;

    std::string init_filename(argv[1]);
    std::ifstream file(init_filename);

    for (std::string line; std::getline(file, line); ) 
    {
        if( line.substr(0,2) == "#L" ){
            L = std::stoi(line.substr(2));
        }
        else if (line.substr(0,2) == "#p")
        {
            p = std::stod(line.substr(2));
        }
        else{
            throw std::invalid_argument("Invalid values in initialization file. Must be only L and p.");
        }
    }

// Task 1
    // lattice l(L,p);
    // l.print_lattice_params();
    // l.print_lattice();

// Task 2
    lattice l(L,p);
    
    std::pair<uint, uint> coord = {0,0};
    l.set_lattice_site(coord, site::flea);



    

    



    return 0;
}