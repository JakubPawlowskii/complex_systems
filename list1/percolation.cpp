#include <iostream>
#include <fstream>
#include "prng.hpp"
#include "lattices2d.hpp"
#include "timer.hpp"

enum class site
{
    empty,
    tree,
    burned_tree,
    fire
};

std::ostream &operator<<(std::ostream &out, const site &site)
{
    switch (site)
    {
    case site::empty:
        out << 0;
        break;
    case site::tree:
        out << 1;
        break;
    case site::burned_tree:
        out << 2;
        break;
    case site::fire:
        out << 3;
        break;
    }
    return out;
}

int main(int argc, char** argv)
{


    return 0;
}