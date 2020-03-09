#pragma once

/*
 *
 * Some template functions, that work as usual math functions, but type specific
 * To work properly, these types have to implement:
 *  - Constructor(0)
 *  - >, <, =
 *  - multiplication with floating point number
 */

namespace TYPESPECIFIC{
    template<class T>
        T abs(const T& number){
            return number >= T(0) ? number : (-1.0) * number;
        }
}
