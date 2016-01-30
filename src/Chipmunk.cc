#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <fstream>

#include "Gauss_Legendre.hh"
#include "Homogeneous_Problem.hh"
#include "Parser.hh"
#include "Sn_Stochastic.hh"
#include "Sn_Transport.hh"

int main(int argc, char *argv[])
{
    if (argc != 10)
    {
        cout << "usage: chipmunk [length, sigma_t_1, sigma_t_2, c_1, c_2, chord_length_1, chord_length_2, internal_source, boundary_source]" << endl;
        return 1;
    }
    
    vector<double> args(argc, 0.0);
    for (unsigned i = 1; i < argc; ++i)
    {
        args[i] = atof(argv[i]);
    }
    
    double length = args[1];
    vector<double> sigma_t(&args[2], &args[4]);
    vector<double> scattering_percentage(&args[4], &args[6]);
    vector<double> chord_length(&args[6], &args[8]);
    double internal_source = args[8];
    double boundary_source = args[9];
    
    unsigned
        number_of_cells = 1000,
        number_of_groups = 1,
        number_of_ordinates = 16,
        number_of_materials = 2,
        number_of_benchmarks = 1000,
        max_iterations = 50000;
    
    double tolerance = 1e-6;
    
    string description_path = "chip_desc.txt";
    string output_path = "chip_out.txt";
    
    // srand (time(NULL));
    // string rand_time = to_string(rand() % 100000000);
    // string temp_desc = "io/temp_desc" + rand_time;
    // string temp_out = "io/temp_out" + rand_time;
    
    Homogeneous_Problem homo_prob(number_of_cells,
                                  number_of_groups,
                                  number_of_ordinates,
                                  number_of_materials,
                                  number_of_benchmarks,
                                  max_iterations,
                                  tolerance,
                                  description_path,
                                  output_path);
    
    homo_prob.set_length(length);
    homo_prob.set_sigma_t(sigma_t);
    homo_prob.set_scattering_percentage(scattering_percentage);
    homo_prob.set_chord_length(chord_length);
    homo_prob.set_source(internal_source,
                         boundary_source);
    
    vector<string> benchmark_cases = {"chord", "percentage", "exponential"};
    vector<string> standard_cases = {"atomic", "lp", "lp_mesh"};
    vector<string> mesh_cases = {"mesh", "dist"};
    
    for (unsigned i = 0; i < benchmark_cases.size(); ++i)
    {
        homo_prob.run_benchmark(benchmark_cases[i]);
    }
    
    for (unsigned i = 0; i < standard_cases.size(); ++i)
    {
        homo_prob.run_standard(standard_cases[i]);
    }
    
    for (unsigned i = 0; i < mesh_cases.size(); ++i)
    {
        homo_prob.run_mesh(mesh_cases[i]);
    }

    homo_prob.calculate_statistics();
}
