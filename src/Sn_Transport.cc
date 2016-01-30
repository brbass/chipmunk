#include <vector>
#include <cmath>
#include <iostream>
#include <string>

#include "Sn_Transport.hh"

Sn_Transport::
Sn_Transport(vector<double> &internal_source,
             vector<double> &boundary_sources,
             vector<string> &boundary_conditions,
             vector<double> &cell_length,
             vector<double> &sigma_t,
             vector<double> &sigma_s,
             vector<double> &nu_sigma_f,
             vector<double> &chi,
             vector<double> &ordinates,
             vector<double> &weights,
             unsigned &number_of_cells,
             unsigned &number_of_groups,
             unsigned &number_of_ordinates,
             unsigned &max_iterations,
             double &tolerance)
:
    internal_source(internal_source),
    boundary_sources(boundary_sources),
    boundary_conditions(boundary_conditions),
    cell_length(cell_length),
    sigma_t(sigma_t),
    sigma_s(sigma_s),
    nu_sigma_f(nu_sigma_f),
    chi(chi),
    ordinates(ordinates),
    weights(weights),
    number_of_cells(number_of_cells),
    number_of_groups(number_of_groups),
    number_of_ordinates(number_of_ordinates),
    max_iterations(max_iterations),
    tolerance(tolerance)
{
    number_of_nodes = 2;
}

void Sn_Transport::
calculate_source(vector<double> &q,
                 vector<double> &phi)
{
    // calculate the source
    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        for (unsigned n = 0; n < number_of_nodes; ++n)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                double sum = 0;
                for (unsigned g2 = 0; g2 < number_of_groups; ++g2)
                {
                    unsigned k1 = g2 + number_of_groups * i;
                    unsigned k2 = g2 + number_of_groups * (g + number_of_groups * i);
                    unsigned k3 = n + number_of_nodes * (g2 + number_of_groups * i);
                    unsigned k4 = g + number_of_groups * i;
                    sum += (chi[k4] * nu_sigma_f[k1] + sigma_s[k2]) * phi[k3];
                }
                unsigned k1 = n + number_of_nodes * (g + number_of_groups * i);
                unsigned k2 = g + number_of_groups * i;
			
                q[k1] = (internal_source[k2] + sum) / 2;
            }
        }
    }
}

void Sn_Transport::
lumped_linear_discontinuous_sweep(vector<double> &psi,
                                  vector<double> &q)
{
    vector<double> psi_boundary_sources(number_of_groups * number_of_ordinates, 0);
    
    // boundary condition, x=0
    if (boundary_conditions[0] == "reflected")
    {
        for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned k1 = g + number_of_groups * o;
                unsigned k2 = 0 + number_of_nodes * (g + number_of_groups * ((number_of_ordinates - 1 - o) + number_of_ordinates * 0));
                    
                psi_boundary_sources[k1] = psi[k2];
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
    for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
    {
        for (unsigned g = 0; g < number_of_groups; ++g)
        {
            unsigned i = 0;
            unsigned k1 = g + number_of_groups * i;
            unsigned k2 = g + number_of_groups * o;
            unsigned k3 = number_of_nodes * (g + number_of_groups * i);
            unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
				
            double a1 = ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
            double a2 = ordinates[o] / 2;
            double a3 = -ordinates[o] / 2;
            double a4 = ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
            double s1 = cell_length[i] * (q[k3] / 2) + ordinates[o] * (boundary_sources[k2] + psi_boundary_sources[k2]);
            double s2 = cell_length[i] * (q[k3+1] / 2);
				
            psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
            psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
        }
    }

    // sweep right over cells
    for (unsigned i = 1; i < number_of_cells; ++i)
    {
        for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned k1 = g + i * number_of_groups;
                unsigned k2 = 1 + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i-1)));
                unsigned k3 = number_of_nodes * (g + number_of_groups * i);
                unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
				
                double a1 = ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                double a2 = ordinates[o] / 2;
                double a3 = -ordinates[o] / 2;
                double a4 = ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                double s1 = cell_length[i] * (q[k3] / 2) + ordinates[o] * psi[k2];
                double s2 = cell_length[i] * (q[k3+1] / 2);
				
                psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
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
                unsigned k1 = g + number_of_groups * o;
                unsigned k2 = 1 + number_of_nodes * (g + number_of_groups * ((number_of_ordinates - 1 - o) + number_of_ordinates * (number_of_cells - 1)));
                
                psi_boundary_sources[k1] = psi[k2];
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
    for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
    {
        for (unsigned g = 0; g < number_of_groups; ++g)
        {
            unsigned i = number_of_cells - 1;
            unsigned k1 = g + number_of_groups * i;
            unsigned k2 = g + number_of_groups * o;
            unsigned k3 = number_of_nodes * (g + number_of_groups * i);
            unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
				
            double a1 = -ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
            double a2 = ordinates[o] / 2;
            double a3 = -ordinates[o] / 2;
            double a4 = -ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
            double s1 = cell_length[i] * (q[k3] / 2);
            double s2 = cell_length[i] * (q[k3+1] / 2) - ordinates[o] * (boundary_sources[k2] + psi_boundary_sources[k2]);
				
            psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
            psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);		
        }
    }

    // sweep left over cells
    for (int i = number_of_cells - 2; i >= 0; --i)
    {
        for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
        {
            for (unsigned g = 0; g < number_of_groups; ++g)
            {
                unsigned k1 = g + number_of_groups * i;
                unsigned k2 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * (i+1))); // n=1
                unsigned k3 = number_of_nodes * (g + number_of_groups * i);
                unsigned k4 = number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
				
                double a1 = -ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                double a2 = ordinates[o] / 2;
                double a3 = -ordinates[o] / 2;
                double a4 = -ordinates[o] / 2 + sigma_t[k1] * cell_length[i] / 2;
                double s1 = cell_length[i] * (q[k3] / 2);
                double s2 = cell_length[i] * (q[k3+1] / 2) - ordinates[o] * psi[k2];
				
                psi[k4] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
                psi[k4+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);		
            }
        }
    }
}

void Sn_Transport::
psi_to_phi(vector<double> &phi,
           vector<double> &psi)
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
                    unsigned k1 = n + number_of_nodes * (g + number_of_groups * (o + number_of_ordinates * i));
			
                    sum += weights[o] * psi[k1];
                }
                unsigned k1 = n + number_of_nodes * (g + number_of_groups * i);
		
                phi[k1] = sum;
            }
        }
    }
}

void Sn_Transport::
diffusion_synthetic_acceleration(vector<double> &phi,
                                 vector<double> &phi_old)
{
    if (number_of_groups > 1)
    {
        cout << "Only one group implemented in DSA" << endl;
        return;
    }

    vector<double> phi_half = phi;
    vector<double> f0(number_of_cells * number_of_nodes, 0); // correction to scalar flux
    vector<double> f0_old(number_of_cells * number_of_nodes, 0); // correction from previous iteration
    vector<double> q(number_of_cells * number_of_nodes, 0); // source for dsa calculation
    
    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        for (unsigned n = 0; n < number_of_nodes; ++n)
        {
            double sum = 0;
            unsigned k = n + number_of_nodes * i;
            
            sum += (nu_sigma_f[i] + sigma_s[i]) * (phi_half[k] - phi_old[k]);
            
            q[k] = sum;
        }
    }

    bool converged = false;
    
    for (unsigned it = 0; it < 2; ++it)
    {
        f0_old = f0;

        diffusion_synthetic_acceleration_sweep(f0,
                                               q);
        
        check_convergence(converged,
                          f0,
                          f0_old);
        
        if (converged)
        {
            break;
        }
        else if (it==max_iterations-1)
        {
            cout << "DSA failed to converge" << endl;
        }
    }
    
    
    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        for (unsigned n = 0; n < number_of_nodes; ++n)
        {
            unsigned k = n + number_of_nodes * i;
            
            phi[k] = phi_half[k] + f0[k];
        }
    }
}

void Sn_Transport::
diffusion_synthetic_acceleration_sweep(vector<double> &f0,
                                       vector<double> &q)
{
    // calculate incoming current
    vector<double> jinc(2, 0);
    for (unsigned o = 0; o < number_of_ordinates / 2; ++o)
    {
        jinc[0] += weights[o]*ordinates[o]*boundary_sources[o];
    }
    for (unsigned o = number_of_ordinates / 2; o < number_of_ordinates; ++o)
    {
        jinc[1] += weights[o]*ordinates[o]*boundary_sources[o];
    }

    // boundary condition, x=0
    double alpha = 0;
    if (boundary_conditions[0] == "reflected")
    {
        alpha = 1;
    }
    else if (boundary_conditions[0] == "vacuum")
    {
        alpha = 0;
    }
        
    // sweep first cell
    {
        unsigned i = 0;
        unsigned k = number_of_nodes * i;

        double diff = 1 / (3*sigma_t[i]);
        
        double a1 = (1+alpha)*cell_length[i] / 2 + (1-alpha)*diff + pow(cell_length[i], 2) * (sigma_t[i] - sigma_s[i]);
        double a2 = -(1-alpha)*diff;
        double a3 = 0;
        double a4 = pow(cell_length[i], 2) * (sigma_t[i] - sigma_s[i]);
        double s1 = pow(cell_length[i], 2) * q[k] + 2 * cell_length[i] * jinc[0];
        double s2 = pow(cell_length[i], 2) * q[k+1];
        
        f0[k] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
        f0[k+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
    }
    
    // sweep right over cells
    for (unsigned i = 1; i < number_of_cells; ++i)
    {
        unsigned k1 = number_of_nodes * i;
        unsigned k2 = number_of_nodes * (i-1);

        double diff = 1 / (3*sigma_t[i]);
        double diffm1 = 1 / (3*sigma_t[i-1]);
        
        double a1 = cell_length[i] / 2 + diff + pow(cell_length[i], 2) * (sigma_t[i] - sigma_s[i]);
        double a2 = -diff;
        double a3 = 0;
        double a4 = pow(cell_length[i], 2) * (sigma_t[i] - sigma_s[i]);
        double s1 = pow(cell_length[i], 2) * q[k1] + diffm1*f0[k2] + (cell_length[i-1] / 2 - diffm1)*f0[k2+1];  
        double s2 = pow(cell_length[i], 2) * q[k1+1];
        
        f0[k1] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
        f0[k1+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
    }

    // boundary condition, x=X
    if (boundary_conditions[1] == "reflected")
    {
        alpha = 1;
    }
    else if (boundary_conditions[1] == "vacuum")
    {
        alpha = 0;
    }
    
    // sweep final cell
    {
        unsigned i = number_of_cells - 1;
        unsigned k = number_of_nodes * i;
				
        double diff = 1 / (3*sigma_t[i]);

        double a1 = pow(cell_length[i], 2)*(sigma_t[i]-sigma_s[i]);
        double a2 = 0;
        double a3 = -(1-alpha)*diff;
        double a4 = (1+alpha)/2*cell_length[i] + (1-alpha)*diff + pow(cell_length[i], 2)*(sigma_t[i]-sigma_s[i]);
        double s1 = pow(cell_length[i], 2)*q[k];
        double s2 = pow(cell_length[i], 2)*q[k+1] - 2*cell_length[i] * jinc[1];
        
        f0[k] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
        f0[k+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);		
    }

    // sweep left over cells
    for (int i = number_of_cells - 2; i >= 0; --i)
    {
        unsigned k1 = number_of_nodes * i;
        unsigned k2 = number_of_nodes * (i+1);
        
        double diff = 1 / (3*sigma_t[i]);
        double diffp1 = 1 / (3*sigma_t[i+1]);
        
        double a1 = pow(cell_length[i], 2) * (sigma_t[i] - sigma_s[i]);
        double a2 = 0;
        double a3 = - diff;
        double a4 = cell_length[i]/2 + diff+ pow(cell_length[i], 2)*(sigma_t[i] - sigma_s[i]);
        double s1 = pow(cell_length[i], 2) * q[k1];
        double s2 = pow(cell_length[i], 2) * q[k1+1] - (cell_length[i] / 2 + diffp1) * f0[k2] + diffp1*f0[k2+1];
        
        f0[k1] = (a4 * s1 - a2 * s2) / (a1 * a4 - a2 * a3);
        f0[k1+1] = (a3 * s1 - a1 * s2) / (a2 * a3 - a1 * a4);
    }
}

void Sn_Transport::
check_convergence(bool &converged,
                  vector<double> &phi,
                  vector<double> &phi_old)
{
    converged = true;
    
    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        for (unsigned g = 0; g < number_of_groups; ++g)
        {
            for (unsigned n = 0; n < number_of_nodes; ++n)
            {
                unsigned k1 = n + number_of_nodes * (g + number_of_groups * i);
                
                double error_phi = abs(phi[k1] - phi_old[k1]);
                
                if (error_phi > tolerance)
                {
                    converged = false;
                    return;
                }
            }
        }
    }
}

void Sn_Transport::
check_convergence(bool &converged,
                  vector<double> &phi,
                  vector<double> &phi_old,
                  vector<double> &error_phi,
                  vector<double> &error_phi_old)
{
    error_phi_old = error_phi;
    
    converged = true;
    
    for (unsigned i = 0; i < number_of_cells; ++i)
    {
        for (unsigned g = 0; g < number_of_groups; ++g)
        {
            for (unsigned n = 0; n < number_of_nodes; ++n)
            {
                unsigned k1 = n + number_of_nodes * (g + number_of_groups * i);
                
                error_phi[k1] = abs(phi[k1] - phi_old[k1]) / phi[k1];
                
                if (error_phi[k1] > tolerance * (1 - error_phi[k1] / error_phi_old[k1]))
                {
                    converged = false;
                }
            }
        }
    }
}

void Sn_Transport::
lumped_linear_discontinuous(vector<double> &psi,
                            vector<double> &leakage,
                            unsigned &iterations)
{
    using namespace std;
        
    // temporary variables
    vector<double> phi(number_of_cells * number_of_groups * number_of_nodes, 0); // average scalar flux
    vector<double> phi_old(number_of_cells * number_of_groups * number_of_nodes, 0); // average scalar flux from previous iteration
    vector<double> q(number_of_cells * number_of_groups * number_of_nodes, 0); // total source, including fission and scattering
    vector<double> error_phi(number_of_cells * number_of_groups * number_of_nodes, 1);
    vector<double> error_phi_old(number_of_cells * number_of_groups * number_of_nodes, 1);
    bool converged = false;
    
    // begin iterations
    for (unsigned it = 0; it < max_iterations; ++it)
    {
        calculate_source(q,
                         phi);
        
        lumped_linear_discontinuous_sweep(psi,
                                          q);
        
        phi_old = phi;
        psi_to_phi(phi,
                   psi);
        diffusion_synthetic_acceleration(phi,
                                         phi_old);

        check_convergence(converged,
                          phi,
                          phi_old,
                          error_phi,
                          error_phi_old);
        
        if (converged)
        {
            cout << "check phi: " << phi[0] << endl;
            iterations = it + 1;
            cout << "LLD iterations: " << iterations << endl;
            break;
        }
        else if (it==max_iterations-1)
        {
            iterations = it + 1;
            cout << "LLD failed to converge" << endl;
        }

    }
    
    calculate_leakage(psi,
                      leakage);
}

void Sn_Transport::
calculate_leakage(vector<double> &psi,
                  vector<double> &leakage)
{
    unsigned il = 0;
    unsigned ir = number_of_cells - 1;

    leakage.assign(leakage.size(), 0.0);
        
    for (unsigned op = 0; op < number_of_ordinates / 2; ++op)
    {
        unsigned on = number_of_ordinates - op - 1;
        for (unsigned g = 0; g < number_of_groups; ++g)
        {
            unsigned kl = 0 + number_of_nodes * (g + number_of_groups * (on + number_of_ordinates * il));
            unsigned kr = 1 + number_of_nodes * (g + number_of_groups * (op + number_of_ordinates * ir));
                
            leakage[0] += ordinates[on] * weights[on] * psi[kl];
            leakage[1] += ordinates[op] * weights[op] * psi[kr];
        }
    }
}
