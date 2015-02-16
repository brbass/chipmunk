#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>

#include "Parser.hh"
#include "Sn_Stochastic.hh"
#include "Sn_Transport.hh"
#include "Gauss_Legendre.hh"

int main(int argc, char** argv)
{
    using namespace std;

    if (argc != 2)
    {
        cout << "usage: chipmunk input_folder" << endl;
        return 1;
    }

    string input_folder = argv[1]; 
    
    unsigned
        number_of_cells,
        number_of_groups,
        number_of_ordinates,
        number_of_materials,
        number_of_benchmarks,
        max_iterations,
        iterations;
    double
        tolerance;
    vector<double>
        psi,
        psi_total,
        psi_squared,
        psi_squared_total,
        internal_source,
        boundary_sources,
        cell_length,
        sigma_t,
        sigma_s,
        nu_sigma_f,
        chi,
        ordinates,
        weights,
        chord_length,
        leakage;
    vector<string>
        boundary_conditions,
        methods;
    string
        problem_type;
    unsigned number_of_nodes = 2;
    
    parse_data(input_folder + "/problem_type", problem_type);

    if (problem_type == "stochastic_material")
    {
        parse_folder(input_folder,
                     number_of_cells,
                     number_of_groups,
                     number_of_ordinates,
                     number_of_materials,
                     number_of_benchmarks,
                     max_iterations,
                     tolerance,
                     psi,
                     psi_total,
                     internal_source,
                     boundary_sources,
                     boundary_conditions,
                     cell_length,
                     sigma_t,
                     sigma_s,
                     nu_sigma_f,
                     chi,
                     chord_length,
                     methods);

        Gauss_Legendre ordinate_set(number_of_ordinates);
        ordinate_set.get_ordinates(ordinates,weights);
        
        Sn_Stochastic SS(internal_source,
                         boundary_sources,
                         boundary_conditions,
                         cell_length,
                         sigma_t,
                         sigma_s,
                         nu_sigma_f,
                         chi,
                         ordinates,
                         weights,
                         chord_length,
                         number_of_cells,
                         number_of_groups,
                         number_of_ordinates,
                         number_of_materials,
                         number_of_benchmarks,
                         max_iterations,
                         tolerance);

        vector<double> transition_probability_dist(number_of_cells * number_of_materials, 0);
        
        SS.run_constant_mesh(psi,
                             psi_total,
                             leakage,
                             transition_probability_dist,
                             iterations);

        output_data(input_folder + "/transition_probability_dist", transition_probability_dist);
        output_data(input_folder + "/psi", psi);
        output_data(input_folder + "/leakage", leakage);
        output_data(input_folder + "/psi_total", psi_total);
        output_data(input_folder + "/iterations", iterations);
        output_data(input_folder + "/ordinates", ordinates);
        output_data(input_folder + "/weights", weights);
    }
    else if (problem_type == "stochastic_benchmark")
    {
        parse_folder(input_folder,
                     number_of_cells,
                     number_of_groups,
                     number_of_ordinates,
                     number_of_materials,
                     number_of_benchmarks,
                     max_iterations,
                     tolerance,
                     psi,
                     psi_total,
                     internal_source,
                     boundary_sources,
                     boundary_conditions,
                     cell_length,
                     sigma_t,
                     sigma_s,
                     nu_sigma_f,
                     chi,
                     chord_length,
                     methods);

        
        Gauss_Legendre ordinate_set(number_of_ordinates);
        ordinate_set.get_ordinates(ordinates,weights);
        
        Sn_Stochastic SS(internal_source,
                         boundary_sources,
                         boundary_conditions,
                         cell_length,
                         sigma_t,
                         sigma_s,
                         nu_sigma_f,
                         chi,
                         ordinates,
                         weights,
                         chord_length,
                         number_of_cells,
                         number_of_groups,
                         number_of_ordinates,
                         number_of_materials,
                         number_of_benchmarks,
                         max_iterations,
                         tolerance);

        vector<double> phi_benchmark(number_of_nodes * number_of_cells * number_of_groups * number_of_materials, 0);
        vector<double> phi_benchmark_total(number_of_nodes * number_of_cells * number_of_groups, 0);
        vector<double> leakage_benchmark(2, 0);
        
        for (unsigned d = 0; d < methods.size(); ++d)
        {
            vector<double> psi_temp(psi);
            vector<double> psi_temp_total(psi_total);
            vector<double> leakage_temp(2, 0.0);
            vector<double> phi_temp(number_of_nodes * number_of_cells * number_of_groups * number_of_materials, 0);
            vector<double> phi_temp_total(number_of_nodes * number_of_cells * number_of_groups, 0);
            vector<double> error_phi(number_of_materials, 0);
            double error_phi_total = 0;
            vector<double> error_leakage(number_of_materials, 0);
            
            if (methods[d] == "chord" || methods[d] == "percentage" || methods[d] == "exponential")
            {
                vector<double> psi_temp_squared(psi);
                vector<double> psi_temp_squared_total(psi_total);
                vector<double> leakage_squared_temp(2, 0.0);
                vector<unsigned> iterations_temp(number_of_benchmarks, 0);
                vector<double> phi_temp_squared(phi_temp);
                vector<double> phi_temp_squared_total(phi_temp_total);
                
                if (methods[d] == "chord")
                {
                    SS.run_chord_benchmark(psi_temp,
                                           psi_temp_total,
                                           psi_temp_squared,
                                           psi_temp_squared_total,
                                           phi_temp_squared,
                                           phi_temp_squared_total,
                                           leakage_temp,
                                           leakage_squared_temp,
                                           iterations_temp);
                }
                else if (methods[d] == "percentage")
                {
                    SS.run_percentage_benchmark(psi_temp,
                                                psi_temp_total,
                                                psi_temp_squared,
                                                psi_temp_squared_total,
                                                phi_temp_squared,
                                                phi_temp_squared_total,
                                                leakage_temp,
                                                leakage_squared_temp,
                                                iterations_temp);
                }
                else if (methods[d] == "exponential")
                {
                    SS.run_exponential_benchmark(psi_temp,
                                                 psi_temp_total,
                                                 psi_temp_squared,
                                                 psi_temp_squared_total,
                                                 phi_temp_squared,
                                                 phi_temp_squared_total,
                                                 leakage_temp,
                                                 leakage_squared_temp,
                                                 iterations_temp);
                }

                SS.psi_to_phi_stochastic(phi_temp,
                                         psi_temp);

                SS.psi_to_phi(phi_temp_total,
                              psi_temp_total);
                
                vector<double> std_deviation_leakage_temp(2, 0.0);
                vector<double> std_deviation_psi_temp(psi);
                vector<double> std_deviation_psi_temp_total(psi_total);
                vector<double> std_deviation_phi_temp(phi_temp);
                vector<double> std_deviation_phi_temp_total(phi_temp_total);
                                
                SS.calculate_std_deviation(psi_temp,
                                           psi_temp_squared,
                                           std_deviation_psi_temp);

                SS.calculate_std_deviation(psi_temp_total,
                                           psi_temp_squared_total,
                                           std_deviation_psi_temp_total);

                SS.calculate_std_deviation(phi_temp,
                                           phi_temp_squared,
                                           std_deviation_phi_temp);

                SS.calculate_std_deviation(phi_temp_total,
                                           phi_temp_squared_total,
                                           std_deviation_phi_temp_total);
                
                SS.calculate_std_deviation(leakage_temp,
                                           leakage_squared_temp,
                                           std_deviation_leakage_temp);
                
                if (d == 0)
                {
                    phi_benchmark = phi_temp;
                    phi_benchmark_total = phi_temp_total;
                    leakage_benchmark = leakage_temp;
                }
                
                SS.calculate_error_phi(error_phi,
                                       error_phi_total,
                                       phi_benchmark,
                                       phi_benchmark_total,
                                       phi_temp,
                                       phi_temp_total);
                
                SS.calculate_error_leakage(error_leakage,
                                           leakage_benchmark,
                                           leakage_temp);
                
                output_data(input_folder + "/" + methods[d] + "/psi_squared", psi_temp_squared);
                output_data(input_folder + "/" + methods[d] + "/psi_squared_total", psi_temp_squared_total);
                output_data(input_folder + "/" + methods[d] + "/leakage_squared", leakage_squared_temp);
                output_data(input_folder + "/" + methods[d] + "/std_deviation_psi", std_deviation_psi_temp);
                output_data(input_folder + "/" + methods[d] + "/std_deviation_psi_total", std_deviation_psi_temp_total);
                output_data(input_folder + "/" + methods[d] + "/std_deviation_leakage", std_deviation_leakage_temp);
                output_data(input_folder + "/" + methods[d] + "/iterations", iterations_temp);
                output_data(input_folder + "/" + methods[d] + "/phi_squared", phi_temp_squared);
                output_data(input_folder + "/" + methods[d] + "/phi_squared_total", phi_temp_squared_total);
                output_data(input_folder + "/" + methods[d] + "/std_deviation_phi", std_deviation_phi_temp);
                output_data(input_folder + "/" + methods[d] + "/std_deviation_phi_total", std_deviation_phi_temp_total);
            }
            else if (methods[d] == "atomic" || methods[d] == "lp" || methods[d] == "mesh" || methods[d] == "dist" || methods[d] =="lp_mesh" || methods[d].find("skip") != string::npos)
            {
                unsigned iterations_temp = 0;
            
                if (methods[d] == "atomic")
                {
                    SS.run_atomic_mix(psi_temp,
                                      psi_temp_total,
                                      leakage_temp,
                                      iterations_temp);
                }
                else if (methods[d] == "lp")
                {
                    SS.run_levermore_pomraning(psi_temp,
                                               psi_temp_total,
                                               leakage_temp,
                                               iterations_temp);
                }
                else if (methods[d] == "lp_mesh")
                {
                    SS.run_lp_mesh(psi_temp,
                                   psi_temp_total,
                                   leakage_temp,
                                   iterations_temp);
                }
                else if (methods[d] == "mesh")
                {
                    vector<double> transition_probability_dist_temp(number_of_cells * number_of_materials, 0);
                
                    SS.run_constant_mesh(psi_temp,
                                         psi_temp_total,
                                         leakage_temp,
                                         transition_probability_dist_temp,
                                         iterations_temp);

                    output_data(input_folder + "/" + methods[d] + "/transition_probability", transition_probability_dist_temp);
                }
                else if (methods[d] == "dist")
                {
                    vector<double> transition_probability_dist_temp(number_of_cells * number_of_materials, 0);
                    bool distribution = true;
                
                    SS.run_constant_mesh(psi_temp,
                                         psi_temp_total,
                                         leakage_temp,
                                         transition_probability_dist_temp,
                                         iterations_temp,
                                         distribution);

                    output_data(input_folder + "/" + methods[d] + "/transition_probability", transition_probability_dist_temp);
                }
                else if (methods[d].find("skip") != string::npos)
                {
                    unsigned cell_mixing;

                    if (methods[d] == "skip5")
                    {
                        cell_mixing = 5;
                    }
                    else if (methods[d] == "skip10")
                    {
                        cell_mixing = 10;
                    }
                    else if (methods[d] == "skip20")
                    {
                        cell_mixing = 20;
                    }
                    else if (methods[d] == "skip50")
                    {
                        cell_mixing = 50;
                    }
                    else if (methods[d] == "skip100")
                    {
                        cell_mixing = 100;
                    }
                    else if (methods[d] == "skip200")
                    {
                        cell_mixing = 200;
                    }

                    vector<double> transition_probability_dist_temp(number_of_cells * number_of_materials, 0);
                    
                    SS.run_constant_mesh(psi_temp,
                                         psi_temp_total,
                                         leakage_temp,
                                         transition_probability_dist_temp,
                                         iterations_temp,
                                         false,
                                         cell_mixing);

                    output_data(input_folder + "/" + methods[d] + "/transition_probability", transition_probability_dist_temp);
                }

                SS.psi_to_phi_stochastic(phi_temp,
                                         psi_temp);

                SS.psi_to_phi(phi_temp_total,
                              psi_temp_total);

                if (d == 0)
                {
                    phi_benchmark = phi_temp;
                    phi_benchmark_total = phi_temp_total;
                    leakage_benchmark = leakage_temp;
                }
                
                SS.calculate_error_phi(error_phi,
                                       error_phi_total,
                                       phi_benchmark,
                                       phi_benchmark_total,
                                       phi_temp,
                                       phi_temp_total);

                SS.calculate_error_leakage(error_leakage,
                                           leakage_benchmark,
                                           leakage_temp);

                output_data(input_folder + "/" + methods[d] + "/iterations", iterations_temp);
            
            }
            else if (methods[d] == "benchmark")
            {
                string benchmark_folder;
                
                parse_data(input_folder + "/benchmark_folder", benchmark_folder);
                
                if (d == 0)
                {
                    parse_data(benchmark_folder + "phi", phi_benchmark);
                    parse_data(benchmark_folder + "phi_total", phi_benchmark_total);
                    parse_data(benchmark_folder + "leakage", leakage_benchmark);

                    parse_data(benchmark_folder + "psi", psi_temp);
                    parse_data(benchmark_folder + "psi_total", psi_temp_total);

                    phi_temp = phi_benchmark;
                    phi_temp_total = phi_benchmark_total;
                    leakage_temp = leakage_benchmark;
                }
            }

            output_data(input_folder + "/" + methods[d] + "/phi", phi_temp);
            output_data(input_folder + "/" + methods[d] + "/phi_total", phi_temp_total);
            output_data(input_folder + "/" + methods[d] + "/psi", psi_temp);
            output_data(input_folder + "/" + methods[d] + "/psi_total", psi_temp_total);
            output_data(input_folder + "/" + methods[d] + "/leakage", leakage_temp);
            output_data(input_folder + "/" + methods[d] + "/error_phi", error_phi);
            output_data(input_folder + "/" + methods[d] + "/error_phi_total", error_phi_total);
            output_data(input_folder + "/" + methods[d] + "/error_leakage", error_leakage);
            
            output_data(input_folder + "/ordinates", ordinates);
            output_data(input_folder + "/weights", weights);
        }
    }
    else
    {
        parse_folder(input_folder,
                     number_of_cells,
                     number_of_groups,
                     number_of_ordinates,
                     max_iterations,
                     tolerance,
                     psi,
                     internal_source,
                     boundary_sources,
                     boundary_conditions,
                     cell_length,
                     sigma_t,
                     sigma_s,
                     nu_sigma_f,
                     chi);

        Gauss_Legendre ordinate_set(number_of_ordinates);
        ordinate_set.get_ordinates(ordinates,weights);
        
        Sn_Transport ST(internal_source,
                        boundary_sources,
                        boundary_conditions,
                        cell_length,
                        sigma_t,
                        sigma_s,
                        nu_sigma_f,
                        chi,
                        ordinates,
                        weights,
                        number_of_cells,
                        number_of_groups,
                        number_of_ordinates,
                        max_iterations,
                        tolerance);
    
        ST.lumped_linear_discontinuous(psi,
                                       leakage,
                                       iterations);

        output_data(input_folder + "/psi", psi);
        output_data(input_folder + "/leakage", leakage);
        output_data(input_folder + "/iterations", iterations);
        output_data(input_folder + "/ordinates", ordinates);
        output_data(input_folder + "/weights", weights);
    }
    
    return 0;
}
