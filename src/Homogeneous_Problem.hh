#ifndef Homogeneous_Problem_hh
#define Homogeneous_Problem_hh

#include <string>
#include <vector>

#include "Sn_Stochastic.hh"

using namespace std;

class Homogeneous_Problem
{
    
public:

    Homogeneous_Problem(unsigned number_of_cells,
                        unsigned number_of_groups,
                        unsigned number_of_ordinates,
                        unsigned number_of_materials,
                        unsigned number_of_benchmarks,
                        unsigned max_iterations,
                        double tolerance);

    ~Homogeneous_Problem();
    
    void set_sigma_t(vector<double> sigma_t);

    void set_scattering_percentage(vector<double> scattering_percentage);

    void set_chord_length(vector<double> chord_length);

    void set_length(double slab_length);
    
    void set_source(double internal_source,
                    double boundary_source);

    void run_mesh(string method);
    void run_standard(string method);
    void run_benchmark(string method);
    void run_skip(unsigned cell_mixing);

    void calculate_statistics(string output_path, string description_path);
    
private:

    unsigned
        number_of_cells_,
        number_of_groups_,
        number_of_ordinates_,
        number_of_materials_,
        number_of_benchmarks_,
        max_iterations_,
        length_psi_,
        length_psi_total_,
        length_phi_,
        length_phi_total_,
        number_of_nodes_ = 2;
    
    double tolerance_;
    
    vector<double>
        internal_source_,
        boundary_sources_,
        cell_length_,
        sigma_t_,
        sigma_s_,
        nu_sigma_f_,
        chi_,
        ordinates_,
        weights_,
        chord_length_;
    
    vector<string>
        methods_,
        boundary_conditions_;

    Sn_Stochastic *sn_stochastic_;
    
    vector<vector<double> > phi_store_;
    vector<vector<double> > phi_store_total_;
    vector<vector<double> > leakage_store_;
    vector<string> method_store_;
    
    void print_average(vector<double> &a);
    void print_average(vector<unsigned> &a);
};

#endif
