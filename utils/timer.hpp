/*
 * @file timer.hpp
 * @author
 * @version 1.0 30.12.2020
 * */

#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>

namespace lqs {

//! Simple class for measuring runtime
/*!
 Source: https://www.learncpp.com/cpp-tutorial/timing-your-code/
 */
class timer {
private:
    //! Type aliases to make accessing nested type easier
    using clock_t = std::chrono::high_resolution_clock; //!< Type alias for high_resolution_cloc from std::chrono
    using second_t = std::chrono::duration<double, std::ratio<1> >; //!< Type alias for one second in std::chrono
    std::chrono::time_point<clock_t> _m_beg; //!< Private members holding the starting point of time measurement. Initialized during object creation.
public:
    timer() { _m_beg = clock_t::now(); } //!< Constructor

    void reset(); //!< Function resetting current time

    [[nodiscard]] double elapsed() const; //!< Function returning the time elapsed between moment of call and _m_beg

    static double absolute_time(); //!< Function returning absolute time, without subtracting _m_beg
};

}
#endif //TIMER_HPP
