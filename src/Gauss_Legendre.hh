#ifndef Gauss_Legendre_hh
#define Gauss_Legendre_hh

#include <vector>
#include <iostream>

using namespace std;

class Gauss_Legendre
{
public:

    Gauss_Legendre(unsigned number_of_ordinates);

    void get_number_of_ordinates(unsigned &number_of_ordinates_temp)
    {
        number_of_ordinates_temp = number_of_ordinates;
    }
    
    void get_ordinates(vector<double> &ordinates_temp,
                       vector<double> &weights_temp)
    {
        weights_temp = weights;
        ordinates_temp = ordinates;
    }
    
private:
 
    void assign_ordinates(vector<double> &ordinates_temp,
                     vector<double> &weights_temp);
    
    unsigned number_of_ordinates;
    vector<double> weights;
    vector<double> ordinates;
};

#endif


