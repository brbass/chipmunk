#ifndef Sn_Stochastic_hh
#define Sn_Stochastic_hh

#include <vector>
#include <random>

#include "Sn_Transport.hh"

using namespace std;

class Sn_Stochastic: public Sn_Transport
{
public:

    // creator
    
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
                  double &tolerance);

    // implementation

    void run_constant_mesh(vector<double> &psi,
                           vector<double> &psi_total,
                           vector<double> &leakage,
                           vector<double> &transition_probability_dist,
                           unsigned &iterations,
                           bool distribution = false,
                           unsigned cell_mixing = 1);

    void run_lp_mesh(vector<double> &psi,
                     vector<double> &psi_total,
                     vector<double> &leakage,
                     unsigned &iterations);
    
    void get_transition_probability_dist(vector<double> &transition_probability_dist,
                                         bool &distribution);

    
    void run_percentage_benchmark(vector<double> &psi_benchmark,
                                  vector<double> &psi_benchmark_total,
                                  vector<double> &psi_benchmark_squared,
                                  vector<double> &psi_benchmark_squared_total,
                                  vector<double> &phi_benchmark_squared,
                                  vector<double> &phi_benchmark_squared_total,
                                  vector<double> &leakage,
                                  vector<double> &leakage_squared,
                                  vector<unsigned> &iterations);

    void run_exponential_benchmark(vector<double> &psi_benchmark,
                                   vector<double> &psi_benchmark_total,
                                   vector<double> &psi_benchmark_squared,
                                   vector<double> &psi_benchmark_squared_total,
                                   vector<double> &phi_benchmark_squared,
                                   vector<double> &phi_benchmark_squared_total,
                                   vector<double> &leakage,
                                   vector<double> &leakage_squared,
                                   vector<unsigned> &iterations);
    
    void run_chord_benchmark(vector<double> &psi_benchmark,
                             vector<double> &psi_benchmark_total,
                             vector<double> &psi_benchmark_squared,
                             vector<double> &psi_benchmark_squared_total,
                             vector<double> &phi_benchmark_squared,
                             vector<double> &phi_benchmark_squared_total,
                             vector<double> &leakage,
                             vector<double> &leakage_squared,
                             vector<unsigned> &iterations);

    void run_atomic_mix(vector<double> &psi_benchmark,
                        vector<double> &psi_benchmark_total,
                        vector<double> &leakage,
                        unsigned &iterations);

    void run_levermore_pomraning(vector<double> &psi_benchmark,
                                 vector<double> &psi_benchmark_total,
                                 vector<double> &leakage,
                                 unsigned &iterations);

    void calculate_std_deviation(vector<double> &psi,
                                 vector<double> &psi_squared_total,
                                 vector<double> &std_deviation);
    
    void psi_to_phi_stochastic(vector<double> &phi,
                               vector<double> &psi);

    void calculate_error_phi(vector<double> &error_phi,
                             double &error_phi_total,
                             vector<double> &phi_benchmark,
                             vector<double> &phi_benchmark_total,
                             vector<double> &phi_temp,
                             vector<double> &phi_temp_total);

    void calculate_error_leakage(vector<double> &error_leakage,
                                 vector<double> &leakage_benchmark,
                                 vector<double> &leakage_temp);
protected:

    // functions

    void get_randomized_data(vector<double> &sigma_t_chord,
                             vector<double> &sigma_s_chord,
                             vector<double> &nu_sigma_f_chord,
                             vector<double> &chi_chord,
                             vector<unsigned> &cell_material);
    
    void get_randomized_slab(unsigned &current_material,
                             unsigned &current_cell,
                             unsigned &end_cell);

    void get_randomized_material(unsigned &current_material);

    // vector operators

    void calculate_stochastic_source(vector<double> &q,
                                     vector<double> &phi,
                                     vector<double> &transition_probability_dist);
    
    void constant_mesh_sweep(vector<double> &psi,
                             vector<double> &q,
                             vector<double> &transition_probability_dist,
                             unsigned cell_mixing = 1);

    void lp_mesh_sweep(vector<double> &psi,
                       vector<double> &q);
    
    void psi_stochastic_to_psi(vector<double> &psi_total,
                               vector<double> &psi);
    
    // data

    vector<double> &chord_length;
    vector<double> transition_probability;
    unsigned &number_of_materials;
    unsigned &number_of_benchmarks;
    default_random_engine generator;
};

#endif
