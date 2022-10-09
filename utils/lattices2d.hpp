#ifndef LATTICES_H
#define LATTICES_H

#include "prng.hpp"
#include <vector>
#include <algorithm>

using coord = std::pair<std::size_t, std::size_t>;
std::ostream &operator<<(std::ostream &out, const coord &c)
{
    out << "(" << c.first << ", " << c.second << ")";
    return out;
}

namespace lattices2d
{

    enum class lattice_type
    {
        square,
        triangular,
        // hexagonal
    };

    enum class bc
    {
        periodic,
        open
    };

    enum class neighbors
    {
        precalculate,
        calculate_on_the_fly
    };

    template <typename T, neighbors mode, bc boundary_conditions>
    class lattice_2d
    {
    protected:
        std::size_t _rows;
        std::size_t _cols;
        std::size_t _n_sites;
        std::vector<T> _lattice_sites;
        lattice_type _lattice_type;
        std::vector<std::vector<coord>> _nn;
        std::vector<std::vector<coord>> _nnn;

    public:
        prng gen;

        virtual ~lattice_2d() = default;
        virtual std::unique_ptr<lattice_2d<T, mode, boundary_conditions>> clone() const = 0;
        virtual std::unique_ptr<lattice_2d<T, mode, boundary_conditions>> move_clone() = 0;

        T &operator()(coord coordinate)
        {
            return _lattice_sites[coordinate.first * _cols + coordinate.second];
        }
        T &at(coord coordinate)
        {
            return _lattice_sites[coordinate.first * _cols + coordinate.second];
        }

        std::size_t get_n_rows()
        {
            return _rows;
        }
        std::size_t get_n_cols()
        {
            return _cols;
        }
        std::size_t get_n_sites()
        {
            return _n_sites;
        }
        std::vector<T> get_lattice_sites()
        {
            return _lattice_sites;
        }
        lattice_type get_lattice_type()
        {
            return _lattice_type;
        }
        std::string get_lattice_type_str()
        {
            switch (_lattice_type)
            {
            case lattice_type::square:
                return "square";
            case lattice_type::triangular:
                return "triangular";
            // case lattices::hexagonal:
            //     return "hexagonal";
            default:
                return "unknown";
            }
        }

        std::vector<coord> get_nn(coord coordinate)
        {
            if constexpr (mode == neighbors::precalculate)
            {
                return _nn[coordinate.first * _cols + coordinate.second];
            }
            else
            {
                return calculate_nn(coordinate);
            }
        }
        std::vector<coord> get_nnn(coord coordinate)
        {
            if constexpr (mode == neighbors::precalculate)
            {
                return _nnn[coordinate.first * _cols + coordinate.second];
            }
            else
            {
                return calculate_nnn(coordinate);
            }
        }

        virtual std::vector<coord> calculate_nn(coord coordinate) = 0;
        virtual std::vector<coord> calculate_nnn(coord coordinate) = 0;

        std::size_t get_nn_number(coord coordinate)
        {
            return _nn[coordinate.first * _cols + coordinate.second].size();
        }
        std::size_t get_nnn_number(coord coordinate)
        {
            return _nnn[coordinate.first * _cols + coordinate.second].size();
        }

        template <typename U>
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

        virtual void print_lattice() = 0;

        void save_lattice_to_file(std::string filename)
        {
            std::ofstream file;
            file.open(filename);
            for (std::size_t i = 0; i < _rows; i++)
            {
                for (std::size_t j = 0; j < _cols; j++)
                {
                    file << _lattice_sites[i * _cols + j] << ",";
                }
                file << std::endl;
            }
            file.close();
        }

        void fill_lattice(std::vector<T> states, std::vector<double> probabilities)
        {
            if (states.size() != probabilities.size())
            {
                throw std::invalid_argument("Initial values and probabilities must have the same size.");
            }

            for (auto p : probabilities)
            {
                if (p < 0 || p > 1)
                {
                    throw std::invalid_argument("Probabilities must be between 0 and 1.");
                }
            }

            if (std::accumulate(probabilities.begin(), probabilities.end(), 0.0) != 1.0)
            {
                throw std::invalid_argument("Probabilities must sum to 1.");
            }

            for (std::size_t i = 0; i < _rows; i++)
            {
                for (std::size_t j = 0; j < _cols; j++)
                {
                    _lattice_sites[i * _cols + j] = random_choice(states, probabilities);
                }
            }
        }

        void calculate_neighbors()
        {
            _nn.resize(_n_sites);
            _nnn.resize(_n_sites);
            for (std::size_t i = 0; i < _rows; i++)
            {
                for (std::size_t j = 0; j < _cols; j++)
                {
                    _nn[i * _cols + j] = calculate_nn({i, j});
                    _nnn[i * _cols + j] = calculate_nnn({i, j});
                }
            }
        }
    };

    template <typename T, neighbors mode, bc boundary_conditions>
    class square_lattice final : public virtual lattice_2d<T, mode, boundary_conditions>
    {
    public:
        square_lattice()
        {
            this->_lattice_type = lattice_type::square;
            this->_rows = 0;
            this->_cols = 0;
            this->_n_sites = 0;
        }
        square_lattice(std::size_t rows, std::size_t cols, std::vector<T> initial_values, std::vector<double> probabilities)
        {
            this->_rows = rows;
            this->_cols = cols;
            this->_n_sites = rows * cols;
            this->_lattice_sites.resize(this->_n_sites);
            this->_lattice_type = lattice_type::square;

            if constexpr (mode == neighbors::precalculate)
            {
                this->calculate_neighbors();
            }
            this->fill_lattice(initial_values, probabilities);
        }

        // copy constructor
        square_lattice(const square_lattice &other)
        {
            this->_lattice_type = other._lattice_type;
            this->_rows = other._rows;
            this->_cols = other._cols;
            this->_n_sites = other._n_sites;
            this->_lattice_sites = other._lattice_sites;
            if constexpr (mode == neighbors::precalculate)
            {
                this->_nn = other._nn;
                this->_nnn = other._nnn;
            }
        }
        // move constructor
        square_lattice(square_lattice &&other) noexcept
        {
            this->_lattice_type = other._lattice_type;
            this->_rows = other._rows;
            this->_cols = other._cols;
            this->_n_sites = other._n_sites;
            this->_lattice_sites = std::move(other._lattice_sites);
            if constexpr (mode == neighbors::precalculate)
            {
                this->_nn = std::move(other._nn);
                this->_nnn = std::move(other._nnn);
            }
        }
        // copy assignment operator
        square_lattice &operator=(const square_lattice &other)
        {
            this->_lattice_type = other._lattice_type;
            this->_rows = other._rows;
            this->_cols = other._cols;
            this->_n_sites = other._n_sites;
            this->_lattice_sites = other._lattice_sites;
            if constexpr (mode == neighbors::precalculate)
            {
                this->_nn = other._nn;
                this->_nnn = other._nnn;
            }
            return *this;
        }
        // move assignment operator
        square_lattice &operator=(square_lattice &&other)
        {
            this->_lattice_type = other._lattice_type;
            this->_rows = other._rows;
            this->_cols = other._cols;
            this->_n_sites = other._n_sites;
            this->_lattice_sites = std::move(other._lattice_sites);
            if constexpr (mode == neighbors::precalculate)
            {
                this->_nn = std::move(other._nn);
                this->_nnn = std::move(other._nnn);
            }
            return *this;
        }

        ~square_lattice() = default;
        virtual std::unique_ptr<lattice_2d<T, mode, boundary_conditions>> clone() const override
        {
            return std::unique_ptr<lattice_2d<T, mode, boundary_conditions>>(new square_lattice<T,mode,boundary_conditions>(*this));
        }
        virtual std::unique_ptr<lattice_2d<T, mode, boundary_conditions>> move_clone() override
        {
            return std::unique_ptr<lattice_2d<T, mode, boundary_conditions>>(new square_lattice<T,mode,boundary_conditions>(std::move(*this)));
        }

        void print_lattice() override
        {
            for (std::size_t i = 0; i < this->_rows; i++)
            {
                for (std::size_t j = 0; j < this->_cols; j++)
                {
                    std::cout << this->_lattice_sites[i * this->_cols + j] << " ";
                }
                std::cout << std::endl;
            }
        }

        std::vector<coord> calculate_nn(coord coordinate) override
        {
            std::vector<coord> nn;
            nn.reserve(4);
            if constexpr (boundary_conditions == bc::open)
            {
                if (coordinate.first > 0)
                {
                    coord neighbor = std::make_pair(coordinate.first - 1, coordinate.second);
                    nn.push_back(neighbor);
                }
                if (coordinate.first < this->_rows - 1)
                {
                    coord neighbor = std::make_pair(coordinate.first + 1, coordinate.second);
                    nn.push_back(neighbor);
                }
                if (coordinate.second > 0)
                {
                    coord neighbor = std::make_pair(coordinate.first, coordinate.second - 1);
                    nn.push_back(neighbor);
                }
                if (coordinate.second < this->_cols - 1)
                {
                    coord neighbor = std::make_pair(coordinate.first, coordinate.second + 1);
                    nn.push_back(neighbor);
                }
            }
            else if constexpr (boundary_conditions == bc::periodic)
            {
                coord neighbor = std::make_pair((coordinate.first + this->_rows - 1) % this->_rows, coordinate.second);
                nn.push_back(neighbor);
                neighbor = std::make_pair((coordinate.first + 1) % this->_rows, coordinate.second);
                nn.push_back(neighbor);
                neighbor = std::make_pair(coordinate.first, (coordinate.second + this->_cols - 1) % this->_cols);
                nn.push_back(neighbor);
                neighbor = std::make_pair(coordinate.first, (coordinate.second + 1) % this->_cols);
                nn.push_back(neighbor);
            }
            return nn;
        }

        std::vector<coord> calculate_nnn(coord coordinate) override
        {
            std::vector<coord> nnn;
            nnn.reserve(4);
            if constexpr (boundary_conditions == bc::open)
            {
                if (coordinate.first > 0 && coordinate.second > 0)
                {
                    coord neighbor = std::make_pair(coordinate.first - 1, coordinate.second - 1);
                    nnn.push_back(neighbor);
                }
                if (coordinate.first > 0 && coordinate.second < this->_cols - 1)
                {
                    coord neighbor = std::make_pair(coordinate.first - 1, coordinate.second + 1);
                    nnn.push_back(neighbor);
                }
                if (coordinate.first < this->_rows - 1 && coordinate.second > 0)
                {
                    coord neighbor = std::make_pair(coordinate.first + 1, coordinate.second - 1);
                    nnn.push_back(neighbor);
                }
                if (coordinate.first < this->_rows - 1 && coordinate.second < this->_cols - 1)
                {
                    coord neighbor = std::make_pair(coordinate.first + 1, coordinate.second + 1);
                    nnn.push_back(neighbor);
                }
            }
            else if (boundary_conditions == bc::periodic)
            {
                coord neighbor = std::make_pair((coordinate.first + this->_rows - 1) % this->_rows, (coordinate.second + this->_cols - 1) % this->_cols);
                nnn.push_back(neighbor);
                neighbor = std::make_pair((coordinate.first + this->_rows - 1) % this->_rows, (coordinate.second + 1) % this->_cols);
                nnn.push_back(neighbor);
                neighbor = std::make_pair((coordinate.first + 1) % this->_rows, (coordinate.second + this->_cols - 1) % this->_cols);
                nnn.push_back(neighbor);
                neighbor = std::make_pair((coordinate.first + 1) % this->_rows, (coordinate.second + 1) % this->_cols);
                nnn.push_back(neighbor);
            }
            return nnn;
        }
    };

    template <typename T, neighbors mode, bc boundary_conditions>
    class triangular_lattice : public lattice_2d<T, mode, boundary_conditions>
    {
    public:
        triangular_lattice()
        {
            this->_lattice_type = lattice_type::triangular;
            this->_rows = 0;
            this->_cols = 0;
            this->_n_sites = 0;
        }
        triangular_lattice(std::size_t rows, std::size_t cols, std::vector<T> initial_values, std::vector<double> probabilities)
        {
            this->_rows = rows;
            this->_cols = cols;
            this->_n_sites = rows * cols;
            this->_lattice_sites.resize(this->_n_sites);
            this->_lattice_type = lattice_type::square;

            if constexpr (mode == neighbors::precalculate)
            {
                this->calculate_neighbors();
            }
            this->fill_lattice(initial_values, probabilities);
        }

        // copy constructor
        triangular_lattice(const triangular_lattice &other)
        {
            this->_rows = other._rows;
            this->_cols = other._cols;
            this->_n_sites = other._n_sites;
            this->_lattice_sites = other._lattice_sites;
            this->_lattice_type = other._lattice_type;
            this->_nn = other._nn;
            this->_nnn = other._nnn;
        }
        // move constructor
        triangular_lattice(triangular_lattice &&other) noexcept
        {
            this->_rows = other._rows;
            this->_cols = other._cols;
            this->_n_sites = other._n_sites;
            this->_lattice_sites = std::move(other._lattice_sites);
            this->_lattice_type = other._lattice_type;
            this->_nn = std::move(other._nn);
            this->_nnn = std::move(other._nnn);
        }

        // copy assignment operator
        triangular_lattice &operator=(const triangular_lattice &other)
        {
            if (this == &other)
                return *this;
            this->_rows = other._rows;
            this->_cols = other._cols;
            this->_n_sites = other._n_sites;
            this->_lattice_sites = other._lattice_sites;
            this->_lattice_type = other._lattice_type;
            this->_nn = other._nn;
            this->_nnn = other._nnn;
            return *this;
        }
        // move assignment operator
        triangular_lattice &operator=(triangular_lattice &&other) noexcept
        {
            if (this == &other)
                return *this;
            this->_rows = other._rows;
            this->_cols = other._cols;
            this->_n_sites = other._n_sites;
            this->_lattice_sites = std::move(other._lattice_sites);
            this->_lattice_type = other._lattice_type;
            this->_nn = std::move(other._nn);
            this->_nnn = std::move(other._nnn);
            return *this;
        }
        ~triangular_lattice() = default;
        virtual std::unique_ptr<lattice_2d<T, mode, boundary_conditions>> clone() const override
        {
            return std::unique_ptr<lattice_2d<T, mode, boundary_conditions>>(new triangular_lattice<T,mode,boundary_conditions>(*this));
        }
        virtual std::unique_ptr<lattice_2d<T, mode, boundary_conditions>> move_clone() override
        {
            return std::unique_ptr<lattice_2d<T, mode, boundary_conditions>>(new triangular_lattice<T,mode,boundary_conditions>(std::move(*this)));
        }

        void print_lattice() override
        {
            std::string offset;
            for (std::size_t i = 0; i < this->_rows; i++)
            {
                std::cout << offset;
                for (std::size_t j = 0; j < this->_cols; j++)
                {
                    std::cout << this->_lattice_sites[i * this->_cols + j] << " ";
                }
                if (i % 2 == 0)
                {
                    offset = " ";
                }
                else
                {
                    offset = "";
                }
                std::cout << std::endl;
            }
        }

        std::vector<coord> calculate_nn(coord coordinate) override
        {
            std::vector<coord> nn;
            nn.reserve(6);
            coord neighbor;
            if constexpr (boundary_conditions == bc::open)
            {

                if (coordinate.first % 2 == 0)
                {
                    if (coordinate.first > 0)
                    {
                        neighbor = std::make_pair(coordinate.first - 1, coordinate.second);
                        nn.push_back(neighbor);
                    }
                    if (coordinate.first > 0 && coordinate.second > 0)
                    {
                        neighbor = std::make_pair(coordinate.first - 1, coordinate.second - 1);
                        nn.push_back(neighbor);
                    }

                    if (coordinate.first < this->_rows - 1)
                    {
                        neighbor = std::make_pair(coordinate.first + 1, coordinate.second);
                        nn.push_back(neighbor);
                    }
                    if (coordinate.first < this->_rows - 1 && coordinate.second > 0)
                    {
                        neighbor = std::make_pair(coordinate.first + 1, coordinate.second - 1);
                        nn.push_back(neighbor);
                    }

                    if (coordinate.second > 0)
                    {
                        neighbor = std::make_pair(coordinate.first, coordinate.second - 1);
                        nn.push_back(neighbor);
                    }
                    if (coordinate.second < this->_cols - 1)
                    {
                        neighbor = std::make_pair(coordinate.first, coordinate.second + 1);
                        nn.push_back(neighbor);
                    }
                }
                else
                {
                    if (coordinate.first > 0)
                    {
                        neighbor = std::make_pair(coordinate.first - 1, coordinate.second);
                        nn.push_back(neighbor);
                    }
                    if (coordinate.first > 0 && coordinate.second < this->_cols - 1)
                    {
                        neighbor = std::make_pair(coordinate.first - 1, coordinate.second + 1);
                        nn.push_back(neighbor);
                    }

                    if (coordinate.first < this->_rows - 1)
                    {
                        neighbor = std::make_pair(coordinate.first + 1, coordinate.second);
                        nn.push_back(neighbor);
                    }
                    if (coordinate.first < this->_rows - 1 && coordinate.second < this->_cols)
                    {
                        neighbor = std::make_pair(coordinate.first + 1, coordinate.second + 1);
                        nn.push_back(neighbor);
                    }

                    if (coordinate.second < this->_cols - 1)
                    {
                        neighbor = std::make_pair(coordinate.first, coordinate.second + 1);
                        nn.push_back(neighbor);
                    }
                    if (coordinate.second > 0)
                    {
                        neighbor = std::make_pair(coordinate.first, coordinate.second - 1);
                        nn.push_back(neighbor);
                    }
                }
            }
            else if constexpr (boundary_conditions == bc::periodic)
            {
                if (coordinate.first % 2 == 0)
                {
                    neighbor = std::make_pair((coordinate.first + this->_rows - 1) % this->_rows, coordinate.second);
                    nn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + this->_rows - 1) % this->_rows, (coordinate.second + this->_cols - 1) % this->_cols);
                    nn.push_back(neighbor);

                    neighbor = std::make_pair((coordinate.first + 1) % this->_rows, coordinate.second);
                    nn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + 1) % this->_rows, (coordinate.second + this->_cols - 1) % this->_cols);
                    nn.push_back(neighbor);

                    neighbor = std::make_pair(coordinate.first, (coordinate.second + 1) % this->_cols);
                    nn.push_back(neighbor);
                    neighbor = std::make_pair(coordinate.first, (coordinate.second + this->_cols - 1) % this->_cols);
                    nn.push_back(neighbor);
                }
                else
                {
                    neighbor = std::make_pair((coordinate.first + this->_rows - 1) % this->_rows, coordinate.second);
                    nn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + this->_rows - 1) % this->_rows, (coordinate.second + 1) % this->_cols);
                    nn.push_back(neighbor);

                    neighbor = std::make_pair((coordinate.first + 1) % this->_rows, coordinate.second);
                    nn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + 1) % this->_rows, (coordinate.second + 1) % this->_cols);
                    nn.push_back(neighbor);

                    neighbor = std::make_pair(coordinate.first, (coordinate.second + 1) % this->_cols);
                    nn.push_back(neighbor);
                    neighbor = std::make_pair(coordinate.first, (coordinate.second + this->_cols - 1) % this->_cols);
                    nn.push_back(neighbor);
                }
            }
            return nn;
        }

        std::vector<coord> calculate_nnn(coord coordinate) override
        {
            std::vector<coord> nnn;
            nnn.reserve(6);
            coord neighbor;
            if constexpr (boundary_conditions == bc::periodic)
            {
                if (coordinate.first % 2 == 0)
                {
                    neighbor = std::make_pair((coordinate.first + this->_rows - 2) % this->_rows, coordinate.second);
                    nnn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + this->_rows - 1) % this->_rows, (coordinate.second + this->_cols - 2) % this->_cols);
                    nnn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + this->_rows - 1) % this->_rows, (coordinate.second + 1) % this->_cols);
                    nnn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + this->_rows + 1) % this->_rows, (coordinate.second + 1) % this->_cols);
                    nnn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + 2) % this->_rows, coordinate.second);
                    nnn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + 1) % this->_rows, (coordinate.second + this->_cols - 2) % this->_cols);
                    nnn.push_back(neighbor);
                }
                else
                {
                    neighbor = std::make_pair((coordinate.first + this->_rows - 2) % this->_rows, coordinate.second);
                    nnn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + this->_rows - 1) % this->_rows, (coordinate.second + this->_cols + 2) % this->_cols);
                    nnn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + this->_rows - 1) % this->_rows, (coordinate.second - 1) % this->_cols);
                    nnn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + 1) % this->_rows, (coordinate.second + 2) % this->_cols);
                    nnn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + 2) % this->_rows, coordinate.second);
                    nnn.push_back(neighbor);
                    neighbor = std::make_pair((coordinate.first + 1) % this->_rows, (coordinate.second + this->_cols - 1) % this->_cols);
                    nnn.push_back(neighbor);
                }
            }
            else if constexpr (boundary_conditions == bc::open)
            {
                if (coordinate.first % 2 == 0)
                {
                    if (coordinate.first > 1)
                    {
                        neighbor = std::make_pair(coordinate.first - 2, coordinate.second);
                        nnn.push_back(neighbor);
                    }
                    if (coordinate.first > 0 && coordinate.second > 1)
                    {
                        neighbor = std::make_pair(coordinate.first - 1, coordinate.second - 2);
                        nnn.push_back(neighbor);
                    }
                    if (coordinate.first > 0 && coordinate.second < this->_cols - 1)
                    {
                        neighbor = std::make_pair(coordinate.first - 1, coordinate.second + 1);
                        nnn.push_back(neighbor);
                    }
                    if (coordinate.first < this->_rows - 1 && coordinate.second < this->_cols - 1)
                    {
                        neighbor = std::make_pair(coordinate.first + 1, coordinate.second + 1);
                        nnn.push_back(neighbor);
                    }
                    if (coordinate.first < this->_rows - 2)
                    {
                        neighbor = std::make_pair(coordinate.first + 2, coordinate.second);
                        nnn.push_back(neighbor);
                    }
                    if (coordinate.first < this->_rows - 1 && coordinate.second > 1)
                    {
                        neighbor = std::make_pair(coordinate.first + 1, coordinate.second - 2);
                        nnn.push_back(neighbor);
                    }
                }
                else
                {
                    if (coordinate.first > 1)
                    {
                        neighbor = std::make_pair(coordinate.first - 2, coordinate.second);
                        nnn.push_back(neighbor);
                    }
                    if (coordinate.first > 0 && coordinate.second < this->_cols - 2)
                    {
                        neighbor = std::make_pair(coordinate.first - 1, coordinate.second + 2);
                        nnn.push_back(neighbor);
                    }
                    if (coordinate.first > 0 && coordinate.second > 0)
                    {
                        neighbor = std::make_pair(coordinate.first - 1, coordinate.second - 1);
                        nnn.push_back(neighbor);
                    }
                    if (coordinate.first < this->_rows - 1 && coordinate.second < this->_cols - 2)
                    {
                        neighbor = std::make_pair(coordinate.first + 1, coordinate.second + 2);
                        nnn.push_back(neighbor);
                    }
                    if (coordinate.first < this->_rows - 2)
                    {
                        neighbor = std::make_pair(coordinate.first + 2, coordinate.second);
                        nnn.push_back(neighbor);
                    }
                    if (coordinate.first < this->_rows - 1 && coordinate.second > 0)
                    {
                        neighbor = std::make_pair(coordinate.first + 1, coordinate.second - 1);
                        nnn.push_back(neighbor);
                    }
                }
            }
            return nnn;
        }
    };

};

#endif // LATTICES_H
