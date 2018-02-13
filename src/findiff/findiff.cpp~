#include "typedef.h"

//Central FD coefficients for 1st derivative 4th accuracy order. Antisymmetric coefficients.
const Real d1_4_2C = -1.0/12; //Coeff for +-2h, where h is the stepsize.
const Real d1_4_1C = 2.0/3; //+-1h

//Central FD coefficients for 2nd derivative 4th accuracty order.
//symmetric coefficients.
const Real d2_4_2C = -1.0/12;
const Real d2_4_1C = 4.0/3;
const Real d2_4_0C = -5.0/2;

//First derivative stencil.
//Might want to enable different accuracy orders later. 
Real der1(Real m2h, Real m1h, Real p1h, Real p2h)
{
        /*
        // debug prints:
        std::cout << m2h << std::endl;
        std::cout << m1h << std::endl;
        std::cout << p1h << std::endl;
        std::cout << p2h << std::endl;
        std::cout << d1_4_2C*m2h << std::endl;
        std::cout << d1_4_1C*m1h << std::endl;
        std::cout << d1_4_1C*p1h << std::endl;
        std::cout << d1_4_2C*p2h << std::endl;
        std::cout << d1_4_2C*(p2h-m2h) << std::endl;
        std::cout << d1_4_1C*(p1h-m1h) << std::endl;
        */
        
        return d1_4_2C*(p2h-m2h)+d1_4_1C*(p1h-m1h);
}

Real der2(Real m2h, Real m1h, Real mid, Real p1h, Real p2h)
{

        return d2_4_2C*(m2h+p2h)+d2_4_1C*(m1h+p1h)+d2_4_0C*mid;
}
