#include <vector>
#include <string>
#include <iostream>

using namespace std;

unsigned check(unsigned size,
               vector<double> &check_vector,
               string vector_name)
{
    if (check_vector.size() != size)
    {
        cout << vector_name << " size incorrect: " << check_vector.size() << " instead of " << size << endl; 
        return 1;
    }
    else
    {
        return 0;
    }
}

unsigned check(unsigned size,
               vector<string> &check_vector,
               string vector_name)
{
    if (check_vector.size() != size)
    {
        cout << vector_name << " size incorrect: " << check_vector.size() << " instead of " << size << endl; 
        return 1;
    }
    else
    {
        return 0;
    }
}
