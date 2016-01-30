#include "Homogeneous_Problem.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Gauss_Legendre.hh"
#include "Sn_Stochastic.hh"

Homogeneous_Problem::
Homogeneous_Problem(unsigned number_of_cells,
                    unsigned number_of_groups,
                    unsigned number_of_ordinates,
                    unsigned number_of_materials,
                    unsigned number_of_benchmarks,
                    unsigned max_iterations,
                    double tolerance,
                    string description_path,
                    string output_path):
    number_of_cells_(number_of_cells),
    number_of_groups_(number_of_groups),
    number_of_ordinates_(number_of_ordinates),
    number_of_materials_(number_of_materials),
    number_of_benchmarks_(number_of_benchmarks),
    max_iterations_(max_iterations),
    tolerance_(tolerance),
    description_file_(description_path, fstream::trunc),
    output_file_(output_path, fstream::app)
{
    vector<double> sigma_t(number_of_groups * number_of_materials, 1.0);
    set_sigma_t(sigma_t, false);
    
    vector<double> scattering_percentage(number_of_groups * number_of_materials * number_of_groups, 0.1);
    set_scattering_percentage(scattering_percentage, false);
    
    vector<double> chord_length (number_of_materials, 1.0);
    set_chord_length(chord_length, false);
    set_source(0.0, 0.0, false);
    set_length(10.0, false);
    
    nu_sigma_f_.assign(number_of_materials_ * number_of_cells_ * number_of_groups_, 0.0);
    chi_.assign(number_of_materials_ * number_of_cells_ * number_of_groups_, 0.0);
    boundary_conditions_.assign(2, "vacuum");
    
    Gauss_Legendre ordinate_set(number_of_ordinates_);
    ordinate_set.get_ordinates(ordinates_,
                               weights_);
    
    sn_stochastic_ = new Sn_Stochastic (internal_source_,
                                        boundary_sources_,
                                        boundary_conditions_,
                                        cell_length_,
                                        sigma_t_,
                                        sigma_s_,
                                        nu_sigma_f_,
                                        chi_,
                                        ordinates_,
                                        weights_,
                                        chord_length_,
                                        number_of_cells_,
                                        number_of_groups_,
                                        number_of_ordinates_,
                                        number_of_materials_,
                                        number_of_benchmarks_,
                                        max_iterations_,
                                        tolerance_);

    length_psi_ = number_of_materials_ * number_of_cells_ * number_of_groups_ * number_of_ordinates_ * number_of_nodes_;
    length_psi_total_ = number_of_cells_ * number_of_groups_ * number_of_ordinates_ * number_of_nodes_;
    length_phi_ = number_of_materials_ * number_of_cells_ * number_of_groups_ * number_of_nodes_;
    length_phi_total_ = number_of_cells_ * number_of_groups_ * number_of_nodes_;
    
    description_file_ << "number_of_cells" << "\t";
    output_file_ << number_of_cells_ << "\t";
    
    description_file_ << "number_of_groups" << "\t";
    output_file_ << number_of_groups_ << "\t";

    description_file_ << "number_of_ordinates" << "\t";
    output_file_ << number_of_ordinates_ << "\t";

    description_file_ << "number_of_materials" << "\t";
    output_file_ << number_of_materials_ << "\t";

    description_file_ << "number_of_benchmarks" << "\t";
    output_file_ << number_of_benchmarks_ << "\t";
}

Homogeneous_Problem::
~Homogeneous_Problem()
{
    delete sn_stochastic_;
}

void Homogeneous_Problem::
set_sigma_t(vector<double> sigma_t, bool write)
{
    sigma_t_.resize(number_of_materials_ * number_of_cells_ * number_of_groups_, 0);

    for (unsigned m = 0; m < number_of_materials_; ++m)
    {
        for (unsigned g = 0; g < number_of_groups_; ++g)
        {
            unsigned k = g + number_of_groups_ * m;

            for (unsigned i = 0; i < number_of_cells_; ++i)
            {
                unsigned k1 = g + number_of_groups_ * (i + number_of_cells_ * m);
                
                sigma_t_[k1] = sigma_t[k];
            }

            if (write)
            {
                description_file_ << "sigma_t" << m << "," << g << "\t";
                output_file_ << sigma_t[k] << "\t";
            }
        }
    }
}

void Homogeneous_Problem::
set_scattering_percentage(vector<double> scattering_percentage, bool write)
{
    sigma_s_.resize(number_of_materials_ * number_of_cells_ * number_of_groups_ * number_of_groups_, 0);
    
    for (unsigned m = 0; m < number_of_materials_; ++m)
    {
        for (unsigned g = 0; g < number_of_groups_; ++g)
        {
            unsigned k = g + number_of_groups_ * m;

            for (unsigned g1 = 0; g1 < number_of_groups_; ++g1)
            {
                for (unsigned i = 0; i < number_of_cells_; ++i)
                {
                    unsigned k1 = g + number_of_groups_ * (i + number_of_cells_ * m);
                    unsigned k2 = g + number_of_groups_ * (g1 + number_of_groups_ * (i + number_of_cells_ * m));
                    
                    sigma_s_[k2] = (scattering_percentage[k] * sigma_t_[k1] / number_of_groups_);
                }
            }

            if (write)
            {
                description_file_ << "c" << m << "," << g << "\t";
                output_file_ << scattering_percentage[k] << "\t";
            }
        }
    }
}

void Homogeneous_Problem::
set_chord_length(vector<double> chord_length, bool write)
{
    chord_length_.resize(number_of_materials_, 0.0);

    for (unsigned m = 0; m < number_of_materials_; ++m)
    {
        chord_length_[m] = chord_length[m];

        if (write)
        {
            description_file_ << "chord_length" << m << "\t";
            output_file_ << chord_length[m] << "\t";
        }
    }
}

void Homogeneous_Problem::
set_length(double slab_length, bool write)
{
    cell_length_.assign(number_of_cells_, slab_length / number_of_cells_);

    if (write)
    {
        description_file_ << "slab_length" << "\t";
        output_file_ << slab_length << "\t";
    }
}

void Homogeneous_Problem::
set_source(double internal_source,
           double boundary_source,
           bool write)
{
    internal_source_.assign(number_of_cells_ * number_of_groups_, internal_source);
    boundary_sources_.assign(number_of_ordinates_ * number_of_groups_, 0.0);
    
    for (unsigned o = 0; o < number_of_ordinates_ / 2; ++o)
    {
        for (unsigned g = 0; g < number_of_groups_; ++g)
        {
            boundary_sources_[g + number_of_groups_ * o] = boundary_source;
        }
    }

    if (write)
    {
        description_file_ << "internal_source" << "\t";
        output_file_ << internal_source << "\t";
        
        description_file_ << "boundary_source" << "\t";
        output_file_ << boundary_source << "\t";
    }
}

void Homogeneous_Problem::
run_benchmark(string method)
{
    vector<double> psi_temp(length_psi_, 1.0);
    vector<double> psi_temp_total(length_psi_total_, 1.0);
    vector<double> psi_temp_squared(length_psi_, 0.0);
    vector<double> psi_temp_squared_total(length_psi_total_, 0.0);
    vector<double> phi_temp_squared(length_phi_, 0.0);
    vector<double> phi_temp_squared_total(length_phi_total_, 0.0);
    vector<double> leakage_temp(2, 0.0);
    vector<double> leakage_squared_temp(2, 0.0);
    vector<unsigned> iterations_temp(number_of_benchmarks_, 0);

    if (method == "chord")
    {
        sn_stochastic_->run_chord_benchmark(psi_temp,
                                            psi_temp_total,
                                            psi_temp_squared,
                                            psi_temp_squared_total,
                                            phi_temp_squared,
                                            phi_temp_squared_total,
                                            leakage_temp,
                                            leakage_squared_temp,
                                            iterations_temp);

    }
    else if (method == "percentage")
    {
        sn_stochastic_->run_percentage_benchmark(psi_temp,
                                                 psi_temp_total,
                                                 psi_temp_squared,
                                                 psi_temp_squared_total,
                                                 phi_temp_squared,
                                                 phi_temp_squared_total,
                                                 leakage_temp,
                                                 leakage_squared_temp,
                                                 iterations_temp);
    }
    else if (method == "exponential")
    {
        sn_stochastic_->run_exponential_benchmark(psi_temp,
                                                  psi_temp_total,
                                                  psi_temp_squared,
                                                  psi_temp_squared_total,
                                                  phi_temp_squared,
                                                  phi_temp_squared_total,
                                                  leakage_temp,
                                                  leakage_squared_temp,
                                                  iterations_temp);
    }
    
    vector<double> phi_temp(length_phi_, 0.0);
    vector<double> phi_temp_total(length_phi_total_, 0.0);
    
    sn_stochastic_->psi_to_phi_stochastic(phi_temp,
                                          psi_temp);
    
    sn_stochastic_->psi_to_phi(phi_temp_total,
                               psi_temp_total);

    vector<double> std_deviation_leakage_temp(2, 0.0);
    vector<double> std_deviation_psi_temp(length_psi_, 0.0);
    vector<double> std_deviation_psi_temp_total(length_psi_total_, 0.0);
    vector<double> std_deviation_phi_temp(length_phi_, 0.0);
    vector<double> std_deviation_phi_temp_total(length_phi_total_, 0.0);
    
    sn_stochastic_->calculate_std_deviation(psi_temp,
                                            psi_temp_squared,
                                            std_deviation_psi_temp);

    sn_stochastic_->calculate_std_deviation(psi_temp_total,
                                            psi_temp_squared_total,
                                            std_deviation_psi_temp_total);

    sn_stochastic_->calculate_std_deviation(phi_temp,
                                            phi_temp_squared,
                                            std_deviation_phi_temp);

    sn_stochastic_->calculate_std_deviation(phi_temp_total,
                                            phi_temp_squared_total,
                                            std_deviation_phi_temp_total);
                
    sn_stochastic_->calculate_std_deviation(leakage_temp,
                                            leakage_squared_temp,
                                            std_deviation_leakage_temp);

    phi_store_.push_back(phi_temp);
    phi_store_total_.push_back(phi_temp_total);
    leakage_store_.push_back(leakage_temp);

    double iterations_av_temp = 0;
    compute_iterations_average(iterations_av_temp, iterations_temp);
    iterations_store_.push_back(iterations_av_temp);
    method_store_.push_back(method);
    
    // print_average(phi_temp_total);
}

void Homogeneous_Problem::
run_standard(string method)
{
    unsigned iterations_temp = 0;
    vector<double> psi_temp(length_psi_, 1.0);
    vector<double> psi_temp_total(length_psi_total_, 1.0);
    vector<double> leakage_temp(2, 0.0);
    
    if (method == "atomic")
    {
        sn_stochastic_->run_atomic_mix(psi_temp,
                                       psi_temp_total,
                                       leakage_temp,
                                       iterations_temp);
    }
    else if (method == "lp")
    {
        sn_stochastic_->run_levermore_pomraning(psi_temp,
                                                psi_temp_total,
                                                leakage_temp,
                                                iterations_temp);
    }
    else if (method == "lp_mesh")
    {
        sn_stochastic_->run_lp_mesh(psi_temp,
                                    psi_temp_total,
                                    leakage_temp,
                                    iterations_temp);
    }

    vector<double> phi_temp(length_phi_, 0.0);
    vector<double> phi_temp_total(length_phi_total_, 0.0);

    
    sn_stochastic_->psi_to_phi_stochastic(phi_temp,
                                          psi_temp);
    
    sn_stochastic_->psi_to_phi(phi_temp_total,
                               psi_temp_total);

    phi_store_.push_back(phi_temp);
    phi_store_total_.push_back(phi_temp_total);
    leakage_store_.push_back(leakage_temp);
    method_store_.push_back(method);
    
    // print_average(phi_temp_total);
}

void Homogeneous_Problem::
run_mesh(string method)
{
    unsigned iterations_temp = 0;
    vector<double> psi_temp(length_psi_, 1.0);
    vector<double> psi_temp_total(length_psi_total_, 1.0);
    vector<double> leakage_temp(2, 0.0);
    vector<double> transition_probability_dist_temp(number_of_cells_ * number_of_materials_, 0);

    if (method == "mesh")
    {
        sn_stochastic_->run_constant_mesh(psi_temp,
                                          psi_temp_total,
                                          leakage_temp,
                                          transition_probability_dist_temp,
                                          iterations_temp);
    }
    else if (method == "dist")
    {
        bool distribution = true;
                
        sn_stochastic_->run_constant_mesh(psi_temp,
                                          psi_temp_total,
                                          leakage_temp,
                                          transition_probability_dist_temp,
                                          iterations_temp,
                                          distribution);
    }

    vector<double> phi_temp(length_phi_, 0.0);
    vector<double> phi_temp_total(length_phi_total_, 0.0);

    sn_stochastic_->psi_to_phi_stochastic(phi_temp,
                                          psi_temp);

    sn_stochastic_->psi_to_phi(phi_temp_total,
                               psi_temp_total);

    phi_store_.push_back(phi_temp);
    phi_store_total_.push_back(phi_temp_total);
    leakage_store_.push_back(leakage_temp);
    method_store_.push_back(method);
    
    // print_average(phi_temp_total);
}

void Homogeneous_Problem::
run_skip(unsigned cell_mixing)
{
    unsigned iterations_temp = 0;
    vector<double> psi_temp(length_psi_, 1.0);
    vector<double> psi_temp_total(length_psi_total_, 1.0);
    vector<double> leakage_temp(2, 0.0);
    vector<double> transition_probability_dist_temp(number_of_cells_ * number_of_materials_, 0);
                    
    sn_stochastic_->run_constant_mesh(psi_temp,
                                      psi_temp_total,
                                      leakage_temp,
                                      transition_probability_dist_temp,
                                      iterations_temp,
                                      false,
                                      cell_mixing);
    
    vector<double> phi_temp(length_phi_, 0.0);
    vector<double> phi_temp_total(length_phi_total_, 0.0);

    sn_stochastic_->psi_to_phi_stochastic(phi_temp,
                                          psi_temp);
    
    sn_stochastic_->psi_to_phi(phi_temp_total,
                               psi_temp_total);
    
    phi_store_.push_back(phi_temp);
    phi_store_total_.push_back(phi_temp_total);
    leakage_store_.push_back(leakage_temp);
    method_store_.push_back("skip" + to_string(cell_mixing));
    
    // print_average(phi_temp_total);
}

void Homogeneous_Problem::
print_average(vector<double> &a)
{
    unsigned num = 0;
    double sum = 0;

    for (unsigned i = 0; i < a.size(); ++i)
    {
        num += 1;
        sum += a[i];
    }

    cout << sum / num << endl;
}

void Homogeneous_Problem::
print_average(vector<unsigned> &a)
{
    unsigned num = 0;
    double sum = 0;

    for (unsigned i = 0; i < a.size(); ++i)
    {
        num += 1;
        sum += a[i];
    }

    cout << 1.0 * sum / num << endl;
}

void Homogeneous_Problem::
calculate_statistics()
{
    for (unsigned m1 = 0; m1 < method_store_.size(); ++m1)
    {
        string description = method_store_[m1];

        for (unsigned a = 0; a < 2; ++a)
        {
            description_file_ << description + "|leakage_" << a << "\t";
            output_file_ << leakage_store_[m1][a] << "\t";
        }

        double av_tot = 0;
        vector<double> av(number_of_materials_, 0);

        compute_phi_average(av_tot, phi_store_total_[m1]);
        compute_phi_average(av, phi_store_[m1]);
        
        description_file_ << description + "|phi_average_total" << "\t";
        output_file_ << av_tot << "\t";
        
        for (unsigned a = 0; a < number_of_materials_; ++a)
        {
            description_file_ << description + "|phi_average_" << a << "\t";
            output_file_ << av[a] << "\t";
        }
        
        for (unsigned m2 = m1 + 1; m2 < method_store_.size(); ++m2)
        {
            vector<double> error_phi(number_of_materials_, 0);
            double error_phi_total = 0;
            vector<double> error_leakage(2, 0);
            
            sn_stochastic_->calculate_error_phi(error_phi,
                                                error_phi_total,
                                                phi_store_[m1],
                                                phi_store_total_[m1],
                                                phi_store_[m2],
                                                phi_store_total_[m2]);
                
            sn_stochastic_->calculate_error_leakage(error_leakage,
                                                    leakage_store_[m1],
                                                    leakage_store_[m2]);
            
            description = method_store_[m1] + "/" + method_store_[m2];
            
            description_file_ << description + "|error_phi_total" << "\t";
            output_file_ << error_phi_total << "\t";

            for (unsigned a = 0; a < number_of_materials_; ++a)
            {
                description_file_ << description + "|error_phi_" << a << "\t";
                output_file_ << error_phi[a] << "\t";
            }

            for (unsigned a = 0; a < 2; ++a)
            {
                description_file_ << description + "|error_leakage_" << a << "\t";
                output_file_ << error_leakage[a] << "\t";
            }
        }
    }
    
    description_file_ << endl;
    output_file_ << endl;
}

void Homogeneous_Problem::
compute_phi_average(vector<double> &average, vector<double> &phi)
{
    average.assign(number_of_materials_, 0);
    unsigned num = phi.size() / number_of_materials_;
    
    for (unsigned m = 0; m < number_of_materials_; ++m)
    {
        for (unsigned i = 0; i < num; ++i)
        {
            average[m] += phi[i];
        }
        
        average[m] /= num;
    }
}

void Homogeneous_Problem::
compute_phi_average(double &average, vector<double> &phi)
{
    average = 0;
    unsigned num = phi.size();
    
    for (unsigned i = 0; i < num; ++i)
    {
        average += phi[i];
    }

    average /= num;
}

void Homogeneous_Problem::
compute_iterations_average(double &average, vector<unsigned> &it)
{
    average = 0;
    unsigned num = it.size();

    for (unsigned i = 0; i < num; ++i)
    {
        average += it[i];
    }

    average /= num;
}
