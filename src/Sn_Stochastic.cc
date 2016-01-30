#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <random>
#include <algorithm>

#include "Sn_Transport.hh"
#include "Sn_Stochastic.hh"

using namespace std;

Sn_Stochastic::
Sn_Stochastic(vector<double> &internal_source,
              vector<double> &boundary_sources,
              vector<string> &boundary_conditions,
              vector<double> &cell_length,
              vector<double> &sigma_t,
              vector<double> &sigma_s,
              vector<double> &nu_sigma_f,
              vector<double> &chi,
              vector<double> &ordinates,
              vector<double> &weights,
              vector<double> &chord_length,
              unsigned &number_of_cells,
              unsigned &number_of_groups,
              unsigned &number_of_ordinates,
              unsigned &number_of_materials,
              unsigned &number_of_benchmarks,
              unsigned &max_iterations,
              double &tolerance)
:
    Sn_Transport(internal_source,
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
                 tolerance),
    chord_length(chord_length),
    number_of_materials(number_of_materials),
    number_of_benchmarks(number_of_benchmarks)
{
    // run checks here?

    // calculate transition probability from chord length
    transition_probability.resize(number_of_materials, 0.0);

    double denominator = 0;
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        denominator += chord_length[m];
    }

    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        transition_probability[m] = chord_length[m] / denominator;
    }
}

void Sn_Stochastic::
calculate_error_phi(vector<double> &error_phi,
                    double &error_phi_total,
                    vector<double> &phi_benchmark,
                    vector<double> &phi_benchmark_total,
                    vector<double> &phi_temp,
                    vector<double> &phi_temp_total)
{
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        double sum = 0;

        for (unsigned i = 0; i < number_of_cells; ++i)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                for (unsigned n = 0; n < number_of_nodes; ++n)
                {
                    unsigned k1 = n + number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                    sum += pow((phi_temp[k1] - phi_benchmark[k1]) / phi_benchmark[k1], 2);
                }
            }
        }

        error_phi[m] = sqrt(sum / (number_of_cells * number_of_nodes));
    }
    
    
    double sum = 0;
    
    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        for (unsigned g = 0; g < number_of_groups; ++g)
        {
            for (unsigned n = 0; n < number_of_nodes; ++n)
            {
                unsigned k1 = n + number_of_nodes * (g + number_of_groups * i);
                sum += pow((phi_temp_total[k1] - phi_benchmark_total[k1]) / phi_benchmark_total[k1], 2);
            }
        }
    }

    error_phi_total = sqrt(sum / (number_of_cells * number_of_nodes));
}

void Sn_Stochastic::
calculate_error_leakage(vector<double> &error_leakage,
                        vector<double> &leakage_benchmark,
                        vector<double> &leakage_temp)
{
    for (unsigned i = 0; i < 2; ++i)
    {
        error_leakage[i] = (leakage_temp[i] - leakage_benchmark[i]) / leakage_benchmark[i];
    }
}

void Sn_Stochastic::
calculate_stochastic_source(vector<double> &q,
                            vector<double> &phi,
                            vector<double> &transition_probability_dist)
{
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        for (unsigned i = 0; i < number_of_cells; ++i)
        {
            for (unsigned n = 0; n < number_of_nodes; ++n)
            {
                for (unsigned g = 0; g < number_of_groups; ++g)
                {
                    double sum = 0;
                    for (unsigned g2 = 0; g2 < number_of_groups; ++g2)
                    {
                        unsigned k1 = g2 + number_of_groups * (i + number_of_cells * m);
                        unsigned k2 = g2 + number_of_groups * (g + number_of_groups * (i + number_of_cells * m));
                        unsigned k3 = n + number_of_nodes * (g2 + number_of_groups * i);
                        unsigned k4 = g + number_of_groups * (i + number_of_cells * m);
                            
                        sum += (chi[k4] * nu_sigma_f[k1] + sigma_s[k2]) * phi[k3];
                    }
                    unsigned k1 = n + number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                    unsigned k2 = g + number_of_groups * i;
                    unsigned k3 = i + number_of_cells * m;
			
                    q[k1] = (internal_source[k2] + sum) / 2 * transition_probability_dist[k3];
                }
            }
        }
    }
}

void Sn_Stochastic::
constant_mesh_sweep(vector<double> &psi,
                    vector<double> &q,
                    vector<double> &transition_probability_dist,
                    unsigned cell_mixing)
{
    vector<double> psi_boundary_sources (number_of_groups * number_of_ordinates, 0);
    
    // boundary condition, x=0
    if (boundary_conditions[0] == "reflected")
    {
        for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                double sum = 0;
                for (unsigned m = 0; m < number_of_materials; ++m)
                {
                    unsigned k1 = 0 + number_of_nodes * (g + number_of_groups * ((number_of_ordinates - 1 - o) + number_of_ordinates * (0 + number_of_cells * m)));

                    sum += psi[k1];
                }
                unsigned k1 = g + number_of_groups * o;

                psi_boundary_sources[k1] = sum;
            }
        }
    }
    else if (boundary_conditions[0] == "vacuum")
    {
        for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned k1 = g + number_of_groups * o;
                                        
                psi_boundary_sources[k1] = 0;
            }
        }
    }
        
    // sweep first cell
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned i = 0;
                unsigned k1 = g + number_of_groups * (i + number_of_cells * m);
                unsigned k2 = g + number_of_groups * o;
                unsigned k3 = number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
                unsigned k5 = i + number_of_cells * m;
                    
                double a1 = ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                double a2 = ordinates[o] / 2;
                double a3 = -ordinates[o] / 2;
                double a4 = ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                double s1 = cell_length[i] * (q[k3] / 2) + ordinates[o] *  (boundary_sources[k2] + psi_boundary_sources[k2]) * transition_probability_dist[k5];
                double s2 = cell_length[i] * (q[k3+1] / 2);
		
                psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
            }
        }
    }

    // sweep right over cells
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        for (unsigned i = 1; i < number_of_cells; ++i)
        {
            for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
            {
                for (unsigned g = 0; g < number_of_groups; ++g)
                {
                    unsigned k1 = g + number_of_groups * (i + number_of_cells * m);
                    unsigned k3 = number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                    unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
                    unsigned k6 = i + number_of_cells * m;
                        
                    double stochastic_edge_source = 0; // part of previous cell's flux to be sweeped in current material

                    for (unsigned m1 = 0; m1 < number_of_materials; ++m1)
                    {
                        unsigned k5 = 1 + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * ((i-1) + number_of_cells * m1)));

                        stochastic_edge_source += psi[k5] * transition_probability_dist[k6] / cell_mixing;
                    }
                    
                    unsigned k5 = 1 + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * ((i-1) + number_of_cells * m)));
                    
                    stochastic_edge_source = psi[k5] * (1 - 1 / cell_mixing) + stochastic_edge_source;
                    
                    // if ((i + 1) % cell_mixing == 0)
                    // {
                    //     for (unsigned m1 = 0; m1 < number_of_materials; ++m1)
                    //     {
                    //         unsigned k5 = 1 + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * ((i-1) + number_of_cells * m1)));
                            
                    //         stochastic_edge_source += psi[k5] * transition_probability_dist[k6];
                    //     }
                    // }
                    // else
                    // {
                    //     unsigned k5 = 1 + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * ((i-1) + number_of_cells * m)));
                        
                    //     stochastic_edge_source = psi[k5];
                    // }
		    
                    double a1 = ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                    double a2 = ordinates[o] / 2;
                    double a3 = -ordinates[o] / 2;
                    double a4 = ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                    double s1 = cell_length[i] * (q[k3] / 2) + ordinates[o] * stochastic_edge_source;
                    double s2 = cell_length[i] * (q[k3+1] / 2);
				
                    psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                    psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
                }
            }
        }
    }

    // boundary condition, x=X
    if (boundary_conditions[1] == "reflected")
    {
        for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {

                double sum = 0;
                for (unsigned m = 0; m < number_of_materials; ++m)
                {
                    unsigned k1 = 1 + number_of_nodes * (g + number_of_groups * ((number_of_ordinates - 1 - o) + number_of_ordinates * ((number_of_cells - 1) + number_of_cells * m)));

                    sum += psi[k1];
                }
                unsigned k1 = g + number_of_groups * o;

                psi_boundary_sources[k1] = sum;
            }
        }
    }
    else if (boundary_conditions[1] == "vacuum")
    {
        for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned k1 = g + number_of_groups * o;
                                        
                psi_boundary_sources[k1] = 0;
            }
        }
    }
   
    // sweep final cell
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned i = number_of_cells - 1;
                unsigned k1 = g + number_of_groups * (i + number_of_cells * m);
                unsigned k2 = g + number_of_groups * o;
                unsigned k3 = number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
                unsigned k5 = i + number_of_cells * m;
                    
                double a1 = -ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                double a2 = ordinates[o] / 2;
                double a3 = -ordinates[o] / 2;
                double a4 = -ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                double s1 = cell_length[i] * (q[k3] / 2);
                double s2 = cell_length[i] * (q[k3+1] / 2) - ordinates[o] * (boundary_sources[k2] + psi_boundary_sources[k2]) * transition_probability_dist[k5];
				
                psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);		
            }
        }
    }
        
    // sweep left over cells
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        for (int  i = number_of_cells - 2; i >= 0; --i)
        {
            for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
            {
                for (unsigned g = 0; g < number_of_groups; ++g)
                {
                    unsigned k1 = g + number_of_groups * (i + number_of_cells * m);
                    unsigned k3 = number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                    unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
                    unsigned k6 = i + number_of_cells * m;
                        
                    double stochastic_edge_source = 0; // part of previous cell's flux to be sweeped in current material

                    for (unsigned m1 = 0; m1 < number_of_materials; ++m1)
                    {
                        unsigned k5 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * ((i+1) + number_of_cells * m1)));

                        stochastic_edge_source += psi[k5] * transition_probability_dist[k6];
                    }

                    unsigned k5 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * ((i+1) + number_of_cells * m)));
                    
                    stochastic_edge_source = psi[k5] * (1 - 1 / cell_mixing) + stochastic_edge_source / cell_mixing; 
                    
                    // if ((i+1) % cell_mixing == 0)
                    // {
                    //     for (unsigned m1 = 0; m1 < number_of_materials; ++m1)
                    //     {
                    //         unsigned k5 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * ((i+1) + number_of_cells * m1)));
			
                    //         stochastic_edge_source += psi[k5] * transition_probability_dist[k6];
                    //     }
                    // }
                    // else
                    // {
                    //     unsigned k5 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * ((i+1) + number_of_cells * m)));
			
                    //     stochastic_edge_source = psi[k5];
                    // }
		    
                    double a1 = -ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                    double a2 = ordinates[o] / 2;
                    double a3 = -ordinates[o] / 2;
                    double a4 = -ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                    double s1 = cell_length[i] * (q[k3] / 2);
                    double s2 = cell_length[i] * (q[k3+1] / 2) - ordinates[o] * stochastic_edge_source;
				
                    psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                    psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);		
                }
            }
        }
    }
}                            

void Sn_Stochastic::
psi_to_phi_stochastic(vector<double> &phi,
                      vector<double> &psi)
{
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        for (unsigned i = 0; i < number_of_cells; ++i)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                for (unsigned n = 0; n < number_of_nodes; ++n)
                {
                    double sum = 0;
		
                    for (unsigned o = 0; o < number_of_ordinates; ++o)
                    {
                        unsigned k1 = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
			
                        sum += weights[o] * psi[k1];
                    }
                    unsigned k1 = n + number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
		
                    phi[k1] = sum;
                }
            }
        }
    }
}

void Sn_Stochastic::
psi_stochastic_to_psi(vector<double> &psi_total,
                      vector<double> &psi)
{
    for (unsigned n = 0; n < number_of_nodes; ++n)
    {
        for (unsigned i = 0; i < number_of_cells; ++i)
        {
            for (unsigned o = 0; o < number_of_ordinates; ++o)
            {
                for (unsigned g = 0; g < number_of_groups; ++g)
                {
                    double sum = 0;
	    
                    for (unsigned m = 0; m < number_of_materials; ++m)
                    {
                        unsigned k1 = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
                        
                        sum += psi[k1];
                    }
	    
                    unsigned k1 = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
	    
                    psi_total[k1] = sum;
                }
            }
        }
    }
}

void Sn_Stochastic::
run_constant_mesh(vector<double> &psi,
                  vector<double> &psi_total,
                  vector<double> &leakage,
                  vector<double> &transition_probability_dist,
                  unsigned &iterations,
                  bool distribution,
                  unsigned cell_mixing)
{
    using namespace std;
    
    // temporary variables
    vector<double> phi (number_of_cells * number_of_groups * number_of_nodes, 0); // average scalar flux
    vector<double> phi_old (number_of_cells * number_of_groups * number_of_nodes, 0); // average scalar flux from previous iteration
    vector<double> q (number_of_materials * number_of_cells * number_of_groups * number_of_nodes, 0); // total source, including fission and scattering
    vector<double> error_phi(number_of_cells * number_of_groups * number_of_nodes, 1);
    vector<double> error_phi_old(number_of_cells * number_of_groups * number_of_nodes, 1);
    bool converged = false;
    
    get_transition_probability_dist(transition_probability_dist,
                                    distribution);
    
    // begin iterations
    for (unsigned it = 0; it < max_iterations; ++it)
    {
        calculate_stochastic_source(q,
                                    phi,
                                    transition_probability_dist);
        
        constant_mesh_sweep(psi,
                            q,
                            transition_probability_dist,
                            cell_mixing);

        phi_old = phi;
        psi_stochastic_to_psi(psi_total,
                              psi);
        psi_to_phi(phi,
                   psi_total);
        
        error_phi_old = error_phi;
        check_convergence(converged,
                          phi,
                          phi_old,
                          error_phi,
                          error_phi_old);
        
        if (converged)
        {
            iterations = it + 1;
            break;
        }
        else if (it==max_iterations-1)
        {
            iterations = it + 1;
            cout << "constant mesh diverged" << endl;
        }
    }

    psi_stochastic_to_psi(psi_total,
                          psi);
    
    calculate_leakage(psi_total,
                      leakage);
}

void Sn_Stochastic::
get_transition_probability_dist(vector<double> &transition_probability_dist,
                                bool &distribution)
{
    if (distribution)
    {
        vector<unsigned> cell_material(number_of_cells * number_of_groups, 0);
        
        for (unsigned b = 0; b < number_of_benchmarks; ++b)
        {
            unsigned current_cell = 0;
            unsigned current_material = number_of_materials; //
            unsigned end_cell = 0;
    
            while (current_cell < number_of_cells)
            {
                get_randomized_slab(current_material,
                                    current_cell,
                                    end_cell);

                uniform_int_distribution<int> distribution(0, 1);
                bool reverse_data = distribution(generator) == 1;
        
                for (unsigned i = current_cell; i < end_cell; ++i)
                {
                    unsigned i1 = i;
                    if (reverse_data)
                    {
                        i1 = number_of_cells - i - 1;
                    }
            
                    for (unsigned g = 0; g < number_of_groups; ++g)
                    {
                        unsigned k1 = g + number_of_groups * i1;

                        cell_material[k1] = current_material;
                    }
                }

                current_cell = end_cell;
            }

            for (unsigned i = 0; i < number_of_cells; ++i)
            {
                unsigned k1 = i + number_of_cells * cell_material[i];

                transition_probability_dist[k1] += 1;
            }
        }

        for (unsigned i = 0; i < number_of_cells; ++i)
        {
            unsigned sum = 0;
            
            for (unsigned m = 0; m < number_of_materials; ++m)
            {
                unsigned k1 = i + number_of_cells * m;
                
                sum += transition_probability_dist[k1];
            }

            for (unsigned m = 0; m < number_of_materials; ++m)
            {
                unsigned k1 = i + number_of_cells * m;

                transition_probability_dist[k1] /= sum;
            }
        }
    }
    else
    {
        for (unsigned m = 0; m < number_of_materials; ++m)
        {
            for (unsigned i = 0; i < number_of_cells; ++i)
            {
                unsigned k1 = i + number_of_cells * m;
                transition_probability_dist[k1] = transition_probability[m];
            }
        }
    }
}

void Sn_Stochastic::
run_percentage_benchmark(vector<double> &psi_benchmark,
                         vector<double> &psi_benchmark_total,
                         vector<double> &psi_benchmark_squared,
                         vector<double> &psi_benchmark_squared_total,
                         vector<double> &phi_benchmark_squared,
                         vector<double> &phi_benchmark_squared_total,
                         vector<double> &leakage,
                         vector<double> &leakage_squared,
                         vector<unsigned> &iterations)
{
    vector<unsigned> number_of_material_cells(number_of_materials, 0.0);
    vector<double> sigma_t_percentage(number_of_cells * number_of_groups, 0.0);
    vector<double> sigma_s_percentage(number_of_cells * number_of_groups * number_of_groups, 0.0);
    vector<double> nu_sigma_f_percentage(number_of_cells * number_of_groups, 0.0);
    vector<double> chi_percentage(number_of_cells * number_of_groups, 0.0);
    vector<unsigned> number_of_psi_averages(number_of_nodes * number_of_cells * number_of_materials * number_of_groups, 0);
    
    unsigned sum = 0;
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        number_of_material_cells[m] = transition_probability[m] * number_of_cells;
        sum += number_of_material_cells[m];
    }
    
    number_of_material_cells[0] += number_of_cells - sum;
    
    vector<unsigned> cell_material(number_of_cells * number_of_groups, 0);
    unsigned current_cell = 0;

    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        unsigned start_cell = current_cell;
        for (unsigned i = start_cell; i < start_cell + number_of_material_cells[m]; ++i)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned k1 = g + number_of_groups * i;
                
                cell_material[k1] = m;
            }
            current_cell += 1;
        }
    }

    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        for (unsigned o = 0; o < number_of_ordinates; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned k1 = g + number_of_groups * i;
                unsigned k2 = g + number_of_groups * (i + number_of_cells * cell_material[k1]);
                    
                sigma_t_percentage[k1] = sigma_t[k2];
                nu_sigma_f_percentage[k1] = nu_sigma_f[k2];
                chi_percentage[k1] = chi[k2];

                for (unsigned g1 = 0; g1 < number_of_groups; ++g1)
                {
                    unsigned k3 = g + number_of_groups * (g1 + number_of_groups * i);
                    unsigned k4 = g + number_of_groups * (g1 + number_of_groups * (i + number_of_cells * cell_material[k1]));;
                        
                    sigma_s_percentage[k3] = sigma_s[k4];
                }
            }
        }
    }

    sigma_t_percentage.swap(sigma_t);
    sigma_s_percentage.swap(sigma_s);
    nu_sigma_f_percentage.swap(nu_sigma_f);
    chi_percentage.swap(chi);

    vector<double> psi_temp(number_of_nodes * number_of_cells * number_of_groups * number_of_ordinates, 0);
    vector<double> leakage_temp(2, 0);
    unsigned iterations_temp = 0;
    unsigned number_unconverged = 0;
    vector<double> phi_temp(number_of_nodes * number_of_cells * number_of_groups, 0);
    
    for (unsigned b = 0; b < number_of_benchmarks; ++b)
    {
        const unsigned seed = b;
            
        shuffle(sigma_t.begin(), sigma_t.end(), default_random_engine(seed));
        shuffle(sigma_s.begin(), sigma_s.end(), default_random_engine(seed));
        shuffle(nu_sigma_f.begin(), nu_sigma_f.end(), default_random_engine(seed));
        shuffle(chi.begin(), chi.end(), default_random_engine(seed));
        shuffle(cell_material.begin(), cell_material.end(), default_random_engine(seed));
        
        lumped_linear_discontinuous(psi_temp,
                                    leakage_temp,
                                    iterations_temp);
        
        iterations[b] = iterations_temp;

        psi_to_phi(phi_temp,
                   psi_temp);
        
        if (iterations_temp < max_iterations)
        {
            unsigned b1 = b - number_unconverged;
        
            for (unsigned n = 0; n < number_of_nodes; ++n)
            {
                for (unsigned i = 0; i < number_of_cells; ++i)
                {
                    for (unsigned g = 0; g < number_of_groups; ++g)
                    {
                        unsigned k2 = g + number_of_groups * i;
                        unsigned k4 = n + number_of_nodes * (g + number_of_groups * (i + number_of_cells * cell_material[k2]));
                        unsigned k5 = n + number_of_nodes * (g + number_of_groups * i);
                        
                        for (unsigned o = 0; o < number_of_ordinates; ++o)
                        {
                            unsigned k1 = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                            unsigned k3 = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * cell_material[k2]))); 
                    
                            psi_benchmark_total[k1] = (psi_temp[k1] + psi_benchmark_total[k1] * b1) / (b1 + 1);
                            psi_benchmark_squared_total[k1] = (pow(psi_temp[k1], 2) + psi_benchmark_squared_total[k1] * b1) / (b1 + 1);

                            psi_benchmark[k3] = (psi_temp[k1] + psi_benchmark[k3] * number_of_psi_averages[k4]) / (number_of_psi_averages[k4] + 1);
                            psi_benchmark_squared[k3] = (pow(psi_temp[k1], 2) + psi_benchmark_squared[k3] * number_of_psi_averages[k4]) / (number_of_psi_averages[k4] + 1);
                        }

                        phi_benchmark_squared_total[k5] = (pow(phi_temp[k5], 2) + phi_benchmark_squared_total[k5] * b1) / (b1 + 1);
                        phi_benchmark_squared[k4] = (pow(phi_temp[k5], 2) + phi_benchmark_squared[k4] * number_of_psi_averages[k4]) / (number_of_psi_averages[k4] + 1);

                        number_of_psi_averages[k4] += 1;
                    }
                }
            }

            for (unsigned s = 0; s < 2; ++s)
            {
                leakage[s] = (leakage_temp[s] + leakage[s] * b1) / (b1 + 1);
                leakage_squared[s] = (pow(leakage_temp[s], 2) + leakage_squared[s] * b1) / (b1 + 1);
            }
        }
        else
        {
            number_unconverged += 1;
        }
    }

    sigma_t_percentage.swap(sigma_t);
    sigma_s_percentage.swap(sigma_s);
    nu_sigma_f_percentage.swap(nu_sigma_f);
    chi_percentage.swap(chi);

    if (number_unconverged > 0)
    {
        cout << "percentage, number unconverged: " << number_unconverged << endl;
    }
}

void Sn_Stochastic::
run_exponential_benchmark(vector<double> &psi_benchmark,
                          vector<double> &psi_benchmark_total,
                          vector<double> &psi_benchmark_squared,
                          vector<double> &psi_benchmark_squared_total,
                          vector<double> &phi_benchmark_squared,
                          vector<double> &phi_benchmark_squared_total,
                          vector<double> &leakage,
                          vector<double> &leakage_squared,
                          vector<unsigned> &iterations)
{
    vector<double> sigma_t_chord(number_of_cells * number_of_groups, 0);
    vector<double> sigma_s_chord(number_of_cells * number_of_groups * number_of_groups, 0);
    vector<double> nu_sigma_f_chord(number_of_cells * number_of_groups, 0);
    vector<double> chi_chord(number_of_cells * number_of_groups, 0);
    vector<unsigned> cell_material(number_of_cells * number_of_groups, 0);
    vector<unsigned> number_of_psi_averages(number_of_nodes * number_of_cells * number_of_materials * number_of_groups, 0);
    vector<double> psi_temp(number_of_nodes * number_of_cells * number_of_groups * number_of_ordinates, 0);
    vector<double> leakage_temp(2, 0);
    unsigned iterations_temp = 0;
    unsigned number_unconverged = 0;
    vector<double> phi_temp(number_of_nodes * number_of_cells * number_of_groups, 0);
    
    for (unsigned b = 0; b < number_of_benchmarks; ++b)
    {
        const unsigned seed = b;
        
        get_randomized_data(sigma_t_chord,
                            sigma_s_chord,
                            nu_sigma_f_chord,
                            chi_chord,
                            cell_material);

        sigma_t_chord.swap(sigma_t);
        sigma_s_chord.swap(sigma_s);
        nu_sigma_f_chord.swap(nu_sigma_f);
        chi_chord.swap(chi);
        
        shuffle(sigma_t.begin(), sigma_t.end(), default_random_engine(seed));
        shuffle(sigma_s.begin(), sigma_s.end(), default_random_engine(seed));
        shuffle(nu_sigma_f.begin(), nu_sigma_f.end(), default_random_engine(seed));
        shuffle(chi.begin(), chi.end(), default_random_engine(seed));
        shuffle(cell_material.begin(), cell_material.end(), default_random_engine(seed));
        
        lumped_linear_discontinuous(psi_temp,
                                    leakage_temp,
                                    iterations_temp);
        
        iterations[b] = iterations_temp;
        
        sigma_t_chord.swap(sigma_t);
        sigma_s_chord.swap(sigma_s);
        nu_sigma_f_chord.swap(nu_sigma_f);
        chi_chord.swap(chi);

        psi_to_phi(phi_temp,
                   psi_temp);
        
        if (iterations_temp < max_iterations)
        {
            unsigned b1 = b - number_unconverged;
            
            for (unsigned n = 0; n < number_of_nodes; ++n)
            {
                for (unsigned i = 0; i < number_of_cells; ++i)
                {
                    for (unsigned g = 0; g < number_of_groups; ++g)
                    {
                        unsigned k2 = g + number_of_groups * i;
                        unsigned k4 = n + number_of_nodes * (g + number_of_groups * (i + number_of_cells * cell_material[k2]));
                        unsigned k5 = n + number_of_nodes * (g + number_of_groups * i);
                        
                        for (unsigned o = 0; o < number_of_ordinates; ++o)
                        {
                            unsigned k1 = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                            unsigned k3 = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * cell_material[k2]))); 
                    
                            psi_benchmark_total[k1] = (psi_temp[k1] + psi_benchmark_total[k1] * b1) / (b1 + 1);
                            psi_benchmark_squared_total[k1] = (pow(psi_temp[k1], 2) + psi_benchmark_squared_total[k1] * b1) / (b1 + 1);

                            psi_benchmark[k3] = (psi_temp[k1] + psi_benchmark[k3] * number_of_psi_averages[k4]) / (number_of_psi_averages[k4] + 1);
                            psi_benchmark_squared[k3] = (pow(psi_temp[k1], 2) + psi_benchmark_squared[k3] * number_of_psi_averages[k4]) / (number_of_psi_averages[k4] + 1);
                        }

                        phi_benchmark_squared_total[k5] = (pow(phi_temp[k5], 2) + phi_benchmark_squared_total[k5] * b1) / (b1 + 1);
                        phi_benchmark_squared[k4] = (pow(phi_temp[k5], 2) + phi_benchmark_squared[k4] * number_of_psi_averages[k4]) / (number_of_psi_averages[k4] + 1);
                        
                        number_of_psi_averages[k4] += 1;
                    }
                }
            }

            for (unsigned s = 0; s < 2; ++s)
            {
                leakage[s] = (leakage_temp[s] + leakage[s] * b1) / (b1 + 1);
                leakage_squared[s] = (pow(leakage_temp[s], 2) + leakage_squared[s] * b1) / (b1 + 1);
            }
        }
        else
        {
            number_unconverged += 1;
        }
    }

    if (number_unconverged > 0)
    {
        cout << "exponential, number unconverged: " << number_unconverged << endl;
    }
}

void Sn_Stochastic::
run_chord_benchmark(vector<double> &psi_benchmark,
                    vector<double> &psi_benchmark_total,
                    vector<double> &psi_benchmark_squared,
                    vector<double> &psi_benchmark_squared_total,
                    vector<double> &phi_benchmark_squared,
                    vector<double> &phi_benchmark_squared_total,
                    vector<double> &leakage,
                    vector<double> &leakage_squared,
                    vector<unsigned> &iterations)
{
    vector<double> psi_temp(number_of_nodes * number_of_cells * number_of_groups * number_of_ordinates, 0);
    vector<double> sigma_t_chord(number_of_cells * number_of_groups, 0);
    vector<double> sigma_s_chord(number_of_cells * number_of_groups * number_of_groups, 0);
    vector<double> nu_sigma_f_chord(number_of_cells * number_of_groups, 0);
    vector<double> chi_chord(number_of_cells * number_of_groups, 0);
    vector<unsigned> cell_material(number_of_cells * number_of_groups, 0);
    vector<unsigned> number_of_psi_averages(number_of_nodes * number_of_cells * number_of_materials * number_of_groups, 0);
    vector<double> leakage_temp(2, 0);
    unsigned iterations_temp = 0;
    unsigned number_unconverged = 0;
    vector<double> phi_temp(number_of_nodes * number_of_cells * number_of_groups, 0);
    
    for (unsigned b = 0; b < number_of_benchmarks; ++b)
    {
        get_randomized_data(sigma_t_chord,
                            sigma_s_chord,
                            nu_sigma_f_chord,
                            chi_chord,
                            cell_material);
        
        sigma_t_chord.swap(sigma_t);
        sigma_s_chord.swap(sigma_s);
        nu_sigma_f_chord.swap(nu_sigma_f);
        chi_chord.swap(chi);
        
        lumped_linear_discontinuous(psi_temp,
                                    leakage_temp,
                                    iterations_temp);
        
        iterations[b] = iterations_temp;
        
        sigma_t_chord.swap(sigma_t);
        sigma_s_chord.swap(sigma_s);
        nu_sigma_f_chord.swap(nu_sigma_f);
        chi_chord.swap(chi);

        psi_to_phi(phi_temp,
                   psi_temp);
        
        if (iterations_temp < max_iterations)
        {
            unsigned b1 = b - number_unconverged;
            
            for (unsigned n = 0; n < number_of_nodes; ++n)
            {
                for (unsigned i = 0; i < number_of_cells; ++i)
                {
                    for (unsigned g = 0; g < number_of_groups; ++g)
                    {
                        unsigned k2 = g + number_of_groups * i;
                        unsigned k4 = n + number_of_nodes * (g + number_of_groups * (i + number_of_cells * cell_material[k2]));
                        unsigned k5 = n + number_of_nodes * (g + number_of_groups * i);
                        
                        for (unsigned o = 0; o < number_of_ordinates; ++o)
                        {
                            unsigned k1 = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
                            unsigned k3 = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * cell_material[k2]))); 
                    
                            psi_benchmark_total[k1] = (psi_temp[k1] + psi_benchmark_total[k1] * b1) / (b1 + 1);
                            psi_benchmark_squared_total[k1] = (pow(psi_temp[k1], 2) + psi_benchmark_squared_total[k1] * b1) / (b1 + 1);

                            psi_benchmark[k3] = (psi_temp[k1] + psi_benchmark[k3] * number_of_psi_averages[k4]) / (number_of_psi_averages[k4] + 1);
                            psi_benchmark_squared[k3] = (pow(psi_temp[k1], 2) + psi_benchmark_squared[k3] * number_of_psi_averages[k4]) / (number_of_psi_averages[k4] + 1);
                        }
                        
                        phi_benchmark_squared_total[k5] = (pow(phi_temp[k5], 2) + phi_benchmark_squared_total[k5] * b1) / (b1 + 1);
                        phi_benchmark_squared[k4] = (pow(phi_temp[k5], 2) + phi_benchmark_squared[k4] * number_of_psi_averages[k4]) / (number_of_psi_averages[k4] + 1);
                        
                        number_of_psi_averages[k4] += 1;
                    }
                }
            }

            for (unsigned s = 0; s < 2; ++s)
            {
                leakage[s] = (leakage_temp[s] + leakage[s] * b1) / (b1 + 1);
                leakage_squared[s] = (pow(leakage_temp[s], 2) + leakage_squared[s] * b1) / (b1 + 1);
            }
        }
        else
        {
            number_unconverged += 1;
        }
    }

    if (number_unconverged > 0)
    {
        cout << "chord, number unconverged: " << number_unconverged << endl;
    }
}

void Sn_Stochastic::
get_randomized_data(vector<double> &sigma_t_chord,
                    vector<double> &sigma_s_chord,
                    vector<double> &nu_sigma_f_chord,
                    vector<double> &chi_chord,
                    vector<unsigned> &cell_material)
{
    unsigned current_cell = 0;
    unsigned current_material = number_of_materials; //
    unsigned end_cell = 0;
    
    while (current_cell < number_of_cells)
    {
        get_randomized_slab(current_material,
                            current_cell,
                            end_cell);

        uniform_int_distribution<int> distribution(0, 1);
        bool reverse_data = distribution(generator) == 1;
        
        for (unsigned i = current_cell; i < end_cell; ++i)
        {
            unsigned i1 = i;
            if (reverse_data)
            {
                i1 = number_of_cells - i - 1;
            }
            
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned k1 = g + number_of_groups * i1;
                unsigned k2 = g + number_of_groups * (i + number_of_cells * current_material);

                cell_material[k1] = current_material;
                
                sigma_t_chord[k1] = sigma_t[k2];
                nu_sigma_f_chord[k1] = nu_sigma_f[k2];
                chi_chord[k1] = chi[k2];
                
                for (unsigned g1 = 0; g1 < number_of_groups; ++g1)
                {
                    unsigned k3 = g + number_of_groups * (g1 + number_of_groups * i1);
                    unsigned k4 = g + number_of_groups * (g1 + number_of_groups * (i + number_of_cells * current_material));;
                        
                    sigma_s_chord[k3] = sigma_s[k4];
                }

            }
        }

        current_cell = end_cell;
    }
}

void Sn_Stochastic::
get_randomized_slab(unsigned &current_material,
                    unsigned &current_cell,
                    unsigned &end_cell)
{
    get_randomized_material(current_material);
    
    double rate_of_occurence = 1 / chord_length[current_material];

    exponential_distribution<double> distribution(rate_of_occurence);

    double slab_length = distribution(generator);


    double sum = 0;
    for (unsigned i = current_cell; i < number_of_cells; ++i)
    {
        sum += cell_length[i];

        if (slab_length < sum || i == number_of_cells -1)
        {
            end_cell = i + 1;
            break;
        }
    }
}

void Sn_Stochastic::
get_randomized_material(unsigned &current_material)
{
    uniform_real_distribution<double> distribution(0.0, 1.0);
    double random_0_1 = distribution(generator);

    double sum = 0.0;
    
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        sum += transition_probability[m];

        if (random_0_1 < sum)
        {
            if (m != current_material)
            {
                current_material = m;
                break;
            }
            else
            {
                get_randomized_material(current_material);
                break;
            }
        }
    }
}

void Sn_Stochastic::
run_atomic_mix(vector<double> &psi_benchmark,
               vector<double> &psi_benchmark_total,
               vector<double> &leakage,
               unsigned &iterations)
{
    vector<double> sigma_t_percentage(number_of_cells * number_of_groups, 0.0);
    vector<double> sigma_s_percentage(number_of_cells * number_of_groups * number_of_groups, 0.0);
    vector<double> nu_sigma_f_percentage(number_of_cells * number_of_groups, 0.0);
    vector<double> chi_percentage(number_of_cells * number_of_groups, 0.0);

    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        for (unsigned g = 0; g < number_of_groups; ++g)
        {
            for (unsigned m = 0; m < number_of_materials; ++m)
            {
                unsigned k1 = g + number_of_groups * i;
                unsigned k2 = g + number_of_groups * (i + number_of_cells * m);
                    
                sigma_t_percentage[k1] += sigma_t[k2] * transition_probability[m];
                nu_sigma_f_percentage[k1] += nu_sigma_f[k2] * transition_probability[m];
                chi_percentage[k1] += chi[k2] * transition_probability[m];

                for (unsigned g1 = 0; g1 < number_of_groups; ++g1)
                {
                    unsigned k3 = g + number_of_groups * (g1 + number_of_groups * i);
                    unsigned k4 = g + number_of_groups * (g1 + number_of_groups * (i + number_of_cells * m));
                        
                    sigma_s_percentage[k3] += sigma_s[k4] * transition_probability[m];
                }
            }
        }
    }

    sigma_t_percentage.swap(sigma_t);
    sigma_s_percentage.swap(sigma_s);
    nu_sigma_f_percentage.swap(nu_sigma_f);
    chi_percentage.swap(chi);
            
    lumped_linear_discontinuous(psi_benchmark_total,
                                leakage,
                                iterations);
    
    sigma_t_percentage.swap(sigma_t);
    sigma_s_percentage.swap(sigma_s);
    nu_sigma_f_percentage.swap(nu_sigma_f);
    chi_percentage.swap(chi);

    if (iterations == max_iterations)
    {
        cout << "atomic mix diverged" << endl;
    }
}

void Sn_Stochastic::
run_levermore_pomraning(vector<double> &psi,
                        vector<double> &psi_total,
                        vector<double> &leakage,
                        unsigned &iterations)					
{
    using namespace std;
        
    // temporary variables
    vector<double> phi (number_of_cells * number_of_groups * number_of_nodes, 0); // average scalar flux
    vector<double> phi_old (number_of_cells * number_of_groups * number_of_nodes, 0); // average scalar flux from previous iteration
    vector<double> q (number_of_materials * number_of_cells * number_of_groups * number_of_nodes, 0); // total source, including fission and scattering
    vector<double> stochastic_source (number_of_materials * number_of_cells * number_of_ordinates * number_of_groups * number_of_nodes, 0);
    vector<double> psi_boundary_sources (number_of_groups * number_of_ordinates, 0);
    vector<double> error_phi(number_of_cells * number_of_groups * number_of_nodes, 1);
    vector<double> error_phi_old(number_of_cells * number_of_groups * number_of_nodes, 1);
    bool converged = false;
    
    // begin iterations
    for (unsigned it = 0; it < max_iterations; ++it)
    {
        // calculate the source
        for (unsigned m = 0; m < number_of_materials; ++m)
        {
            for (unsigned i = 0; i < number_of_cells; ++i)
            {
                for (unsigned n = 0; n < number_of_nodes; ++n)
                {
                    for (unsigned g = 0; g < number_of_groups; ++g)
                    {
                        double sum = 0;
                        for (unsigned g2 = 0; g2 < number_of_groups; ++g2)
                        {
                            unsigned k1 = g2 + number_of_groups * (i + number_of_cells * m);
                            unsigned k2 = g2 + number_of_groups * (g + number_of_groups * (i + number_of_cells * m));
                            unsigned k3 = n + number_of_nodes * (g2 + number_of_groups * i);
			    unsigned k4 = g + number_of_groups * (i + number_of_cells * m);
                            
                            sum += (chi[k4] * nu_sigma_f[k1] + sigma_s[k2]) * phi[k3];
                        }
                        unsigned k1 = n + number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                        unsigned k2 = g + number_of_groups * i;
			
                        q[k1] = (internal_source[k2] + sum) / 2 * transition_probability[m];
                    }
                }
            }
        }

        // calculate stochastic part of source for next iteration
        for (unsigned m = 0; m < number_of_materials; ++m)
        {
            for (unsigned i = 0; i < number_of_cells; ++i)
            {
                for (unsigned o = 0; o < number_of_ordinates; ++o)
                {
                    for (unsigned g = 0; g < number_of_groups; ++g)
                    {
                        for (unsigned n = 0; n < number_of_nodes; ++n)
                        {
                            double sum = 0; 
                            for (unsigned m1 = 0; m1 < number_of_materials; ++m1)
                            {
                                if (m1 != m)
                                {
                                    double sum1 = 0;
                                    
                                    for (unsigned m2 = 0; m2 < number_of_materials; ++m2)
                                    {
                                        if (m2 != m1)
                                        {
                                            sum1 += chord_length[m2];
                                        }
                                    }
                                    
                                    unsigned k2 = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m1)));
                                   
                                    sum += psi[k2] * abs(ordinates[o]) / chord_length[m1] * chord_length[m] / sum1;
                                }
                            }

                            unsigned k1 = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
                            stochastic_source[k1] = sum;// - psi[k1] / chord_length[m];
                        }
                    }
                }
            }
        }

        // boundary condition, x=0
        if (boundary_conditions[0] == "reflected")
        {
            for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
            {
                for (unsigned g = 0; g < number_of_groups; ++g)
                {
                    double sum = 0;
                    for (unsigned m = 0; m < number_of_materials; ++m)
                    {
                        unsigned k1 = 0 + number_of_nodes * (g + number_of_groups * ((number_of_ordinates - 1 - o) + number_of_ordinates * (0 + number_of_cells * m)));

                        sum += psi[k1];
                    }
                    unsigned k1 = g + number_of_groups * o;

                    psi_boundary_sources[k1] = sum;
                }
            }
        }
        else if (boundary_conditions[0] == "vacuum")
        {
            for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
            {
                for (unsigned g = 0; g < number_of_groups; ++g)
                {
                    unsigned k1 = g + number_of_groups * o;
                                        
                    psi_boundary_sources[k1] = 0;
                }
            }
        }
        
        // sweep first cell
        for (unsigned m = 0; m < number_of_materials; ++m)
        {
            for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
            {
                for (unsigned g = 0; g < number_of_groups; ++g)
                {
                    unsigned i = 0;
                    unsigned k1 = g + number_of_groups * (i + number_of_cells * m);
                    unsigned k2 = g + number_of_groups * o;
                    unsigned k3 = number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                    unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
                            
                    double a1 = ordinates[o] / 2 + (sigma_t[k1] + abs(ordinates[o]) / chord_length[m]) * cell_length[i] / 2;
                    double a2 = ordinates[o] / 2;
                    double a3 = -ordinates[o] / 2;
                    double a4 = ordinates[o] / 2 + (sigma_t[k1] + abs(ordinates[o]) / chord_length[m]) * cell_length[i] / 2;
                    double s1 = cell_length[i] * (q[k3] + stochastic_source[k4]) / 2 + ordinates[o] *  (boundary_sources[k2] + psi_boundary_sources[k2]) * transition_probability[m];
                    double s2 = cell_length[i] * (q[k3+1] + stochastic_source[k4+1]) / 2;
		
                    psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                    psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
                }
            }
        }

        // sweep right over cells
        for (unsigned m = 0; m < number_of_materials; ++m)
        {
            for (unsigned i = 1; i < number_of_cells; ++i)
            {
                for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
                {
                    for (unsigned g = 0; g < number_of_groups; ++g)
                    {
                        unsigned k1 = g + number_of_groups * (i + number_of_cells * m);
                        unsigned k3 = number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                        unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
                        unsigned k5 = 1 + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * ((i-1) + number_of_cells * m))); 
                        
                        double stochastic_edge_source = psi[k5]; 

                        double a1 = ordinates[o] / 2 + (sigma_t[k1] + abs(ordinates[o]) / chord_length[m]) * cell_length[i] / 2;
                        double a2 = ordinates[o] / 2;
                        double a3 = -ordinates[o] / 2;
                        double a4 = ordinates[o] / 2 + (sigma_t[k1] + abs(ordinates[o]) / chord_length[m]) * cell_length[i] / 2;
                        double s1 = cell_length[i] * (q[k3] + stochastic_source[k4]) / 2 + ordinates[o] * stochastic_edge_source;
                        double s2 = cell_length[i] * (q[k3+1] + stochastic_source[k4+1]) / 2;
				
                        psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                        psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
                    }
                }
            }
        }

        // boundary condition, x=X
        if (boundary_conditions[1] == "reflected")
        {
            for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
            {
                for (unsigned g = 0; g < number_of_groups; ++g)
                {

                    double sum = 0;
                    for (unsigned m = 0; m < number_of_materials; ++m)
                    {
                        unsigned k1 = 1 + number_of_nodes * (g + number_of_groups * ((number_of_ordinates - 1 - o) + number_of_ordinates * ((number_of_cells - 1) + number_of_cells * m)));

                        sum += psi[k1];
                    }
                    unsigned k1 = g + number_of_groups * o;

                    psi_boundary_sources[k1] = sum;
                }
            }
        }
        else if (boundary_conditions[1] == "vacuum")
        {
            for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
            {
                for (unsigned g = 0; g < number_of_groups; ++g)
                {
                    unsigned k1 = g + number_of_groups * o;
                                        
                    psi_boundary_sources[k1] = 0;
                }
            }
        }
   
        // sweep final cell
        for (unsigned m = 0; m < number_of_materials; ++m)
        {
            for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
            {
                for (unsigned g = 0; g < number_of_groups; ++g)
                {
                    unsigned i = number_of_cells - 1;
                    unsigned k1 = g + number_of_groups * (i + number_of_cells * m);
                    unsigned k2 = g + number_of_groups * o;
                    unsigned k3 = number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                    unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));

                    double a1 = -ordinates[o] / 2 + (sigma_t[k1] + abs(ordinates[o]) / chord_length[m]) * cell_length[i] / 2;
                    double a2 = ordinates[o] / 2;
                    double a3 = -ordinates[o] / 2;
                    double a4 = -ordinates[o] / 2 + (sigma_t[k1] + abs(ordinates[o]) / chord_length[m]) * cell_length[i] / 2;
                    double s1 = cell_length[i] * (q[k3] + stochastic_source[k4]) / 2;
                    double s2 = cell_length[i] * (q[k3+1] + stochastic_source[k4+1]) / 2 - ordinates[o] * (boundary_sources[k2] + psi_boundary_sources[k2]) * transition_probability[m];
				
                    psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                    psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);		
                }
            }
        }
        
        // sweep left over cells
        for (unsigned m = 0; m < number_of_materials; ++m)
        {
            for (int i = number_of_cells - 2; i >= 0; --i)
            {
                for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
                {
                    for (unsigned g = 0; g < number_of_groups; ++g)
                    {
                        unsigned k1 = g + number_of_groups * (i + number_of_cells * m);
                        unsigned k3 = number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                        unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
                        unsigned k5 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * ((i+1) + number_of_cells * m)));
                        
                        double stochastic_edge_source = psi[k5]; 

                        double a1 = -ordinates[o] / 2 + (sigma_t[k1] + abs(ordinates[o]) / chord_length[m]) * cell_length[i] / 2;
                        double a2 = ordinates[o] / 2;
                        double a3 = -ordinates[o] / 2;
                        double a4 = -ordinates[o] / 2 + (sigma_t[k1] + abs(ordinates[o]) / chord_length[m]) * cell_length[i] / 2;
                        double s1 = cell_length[i] * (q[k3] + stochastic_source[k4]) / 2;
                        double s2 = cell_length[i] * (q[k3+1] + stochastic_source[k4+1]) / 2 - ordinates[o] * stochastic_edge_source;
				
                        psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                        psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);		
                    }
                }
            }
        }
        
        // calculate new scalar flux
        phi_old = phi;
        psi_stochastic_to_psi(psi_total,
                              psi);
        psi_to_phi(phi,
                   psi_total);

        // check for convergence
        error_phi_old = error_phi;
        check_convergence(converged,
                          phi,
                          phi_old,
                          error_phi,
                          error_phi_old);
        
        if (converged)
        {
            iterations = it + 1;
            break;
        }
        else if (it==max_iterations-1)
        {
            iterations = it + 1;
            cout << "lp diverged" << endl;
        }
    }

    psi_stochastic_to_psi(psi_total,
                          psi);
    
    calculate_leakage(psi_total,
                      leakage);
}						

void Sn_Stochastic::
calculate_std_deviation(vector<double> &psi,
                        vector<double> &psi_squared,
                        vector<double> &std_deviation)
{
    using namespace std;
    
    for (unsigned i = 0; i < psi.size(); ++i)
    {
        std_deviation[i] = sqrt(psi_squared[i] - pow(psi[i], 2));
    }
}

void Sn_Stochastic::
run_lp_mesh(vector<double> &psi,
            vector<double> &psi_total,
            vector<double> &leakage,
            unsigned &iterations)
{
    using namespace std;
    
    // temporary variables
    vector<double> phi (number_of_cells * number_of_groups * number_of_nodes, 0); // average scalar flux
    vector<double> phi_old (number_of_cells * number_of_groups * number_of_nodes, 0); // average scalar flux from previous iteration
    vector<double> q (number_of_materials * number_of_cells * number_of_groups * number_of_nodes, 0); // total source, including fission and scattering
    vector<double> error_phi(number_of_cells * number_of_groups * number_of_nodes, 1);
    vector<double> error_phi_old(number_of_cells * number_of_groups * number_of_nodes, 1);
    vector<double> transition_probability_dist(number_of_cells * number_of_materials, 0);
    bool distribution = false;

    bool converged = false;
    
    get_transition_probability_dist(transition_probability_dist,
                                    distribution);
    
    // begin iterations
    for (unsigned it = 0; it < max_iterations; ++it)
    {
        calculate_stochastic_source(q,
                                    phi,
                                    transition_probability_dist);
    
        lp_mesh_sweep(psi,
                      q);

        phi_old = phi;
        psi_stochastic_to_psi(psi_total,
                              psi);
        psi_to_phi(phi,
                   psi_total);
        
        error_phi_old = error_phi;
        check_convergence(converged,
                          phi,
                          phi_old,
                          error_phi,
                          error_phi_old);
        
        if (converged)
        {
            iterations = it + 1;
            break;
        }
        else if (it==max_iterations-1)
        {
            iterations = it + 1;
            cout << "lp mesh diverged" << endl;
        }
    }

    psi_stochastic_to_psi(psi_total,
                          psi);
    
    calculate_leakage(psi_total,
                      leakage);
}


void Sn_Stochastic::
lp_mesh_sweep(vector<double> &psi,
              vector<double> &q)
{
    vector<double> psi_boundary_sources (number_of_groups * number_of_ordinates, 0);
    
    // boundary condition, x=0
    if (boundary_conditions[0] == "reflected")
    {
        for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                double sum = 0;
                for (unsigned m = 0; m < number_of_materials; ++m)
                {
                    unsigned k1 = 0 + number_of_nodes * (g + number_of_groups * ((number_of_ordinates - 1 - o) + number_of_ordinates * (0 + number_of_cells * m)));

                    sum += psi[k1];
                }
                unsigned k1 = g + number_of_groups * o;

                psi_boundary_sources[k1] = sum;
            }
        }
    }
    else if (boundary_conditions[0] == "vacuum")
    {
        for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned k1 = g + number_of_groups * o;
                                        
                psi_boundary_sources[k1] = 0;
            }
        }
    }
        
    // sweep first cell
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned i = 0;
                unsigned k1 = g + number_of_groups * (i + number_of_cells * m);
                unsigned k2 = g + number_of_groups * o;
                unsigned k3 = number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
                    
                double a1 = ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                double a2 = ordinates[o] / 2;
                double a3 = -ordinates[o] / 2;
                double a4 = ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                double s1 = cell_length[i] * (q[k3] / 2) + ordinates[o] *  (boundary_sources[k2] + psi_boundary_sources[k2]) * transition_probability[m];
                double s2 = cell_length[i] * (q[k3+1] / 2);
		
                psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
            }
        }
    }

    // sweep right over cells
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        for (unsigned i = 1; i < number_of_cells; ++i)
        {
            for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
            {
                for (unsigned g = 0; g < number_of_groups; ++g)
                {
                    unsigned k1 = g + number_of_groups * (i + number_of_cells * m);
                    unsigned k3 = number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                    unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
                        
                    double stochastic_edge_source = 0; // part of previous cell's flux to be sweeped in current material

                    for (unsigned m1 = 0; m1 < number_of_materials; ++m1)
                    {
                        unsigned k5 = 1 + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * ((i-1) + number_of_cells * m1)));
                        if (m1 == m)
                        {
                            stochastic_edge_source += (1 - cell_length[i] / chord_length[m1]) * psi[k5];
                        }
                        else
                        {
                            stochastic_edge_source += cell_length[i] / chord_length[m1] * psi[k5];
                        }
                    }
		    
                    double a1 = ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                    double a2 = ordinates[o] / 2;
                    double a3 = -ordinates[o] / 2;
                    double a4 = ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                    double s1 = cell_length[i] * (q[k3] / 2) + ordinates[o] * stochastic_edge_source;
                    double s2 = cell_length[i] * (q[k3+1] / 2);
				
                    psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                    psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
                }
            }
        }
    }

    // boundary condition, x=X
    if (boundary_conditions[1] == "reflected")
    {
        for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {

                double sum = 0;
                for (unsigned m = 0; m < number_of_materials; ++m)
                {
                    unsigned k1 = 1 + number_of_nodes * (g + number_of_groups * ((number_of_ordinates - 1 - o) + number_of_ordinates * ((number_of_cells - 1) + number_of_cells * m)));

                    sum += psi[k1];
                }
                unsigned k1 = g + number_of_groups * o;

                psi_boundary_sources[k1] = sum;
            }
        }
    }
    else if (boundary_conditions[1] == "vacuum")
    {
        for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned k1 = g + number_of_groups * o;
                                        
                psi_boundary_sources[k1] = 0;
            }
        }
    }
   
    // sweep final cell
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned i = number_of_cells - 1;
                unsigned k1 = g + number_of_groups * (i + number_of_cells * m);
                unsigned k2 = g + number_of_groups * o;
                unsigned k3 = number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
                    
                double a1 = -ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                double a2 = ordinates[o] / 2;
                double a3 = -ordinates[o] / 2;
                double a4 = -ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                double s1 = cell_length[i] * (q[k3] / 2);
                double s2 = cell_length[i] * (q[k3+1] / 2) - ordinates[o] * (boundary_sources[k2] + psi_boundary_sources[k2]) * transition_probability[m];
				
                psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);		
            }
        }
    }
        
    // sweep left over cells
    for (unsigned m = 0; m < number_of_materials; ++m)
    {
        for (int  i = number_of_cells - 2; i >= 0; --i)
        {
            for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
            {
                for (unsigned g = 0; g < number_of_groups; ++g)
                {
                    unsigned k1 = g + number_of_groups * (i + number_of_cells * m);
                    unsigned k3 = number_of_nodes * (g + number_of_groups * (i + number_of_cells * m));
                    unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i + number_of_cells * m)));
                        
                    double stochastic_edge_source = 0; // part of previous cell's flux to be sweeped in current material

                    for (unsigned m1 = 0; m1 < number_of_materials; ++m1)
                    {
                        unsigned k5 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * ((i+1) + number_of_cells * m1)));

                        if (m1 == m)
                        {
                            stochastic_edge_source += (1 - cell_length[i] / chord_length[m1]) * psi[k5];
                        }
                        else
                        {
                            stochastic_edge_source += cell_length[i] / chord_length[m1] * psi[k5];
                        }
                    }
		    
                    double a1 = -ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                    double a2 = ordinates[o] / 2;
                    double a3 = -ordinates[o] / 2;
                    double a4 = -ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                    double s1 = cell_length[i] * (q[k3] / 2);
                    double s2 = cell_length[i] * (q[k3+1] / 2) - ordinates[o] * stochastic_edge_source;
				
                    psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                    psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);		
                }
            }
        }
    }
}                            

