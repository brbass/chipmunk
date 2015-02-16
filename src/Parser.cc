#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include "Parser.hh"
#include "Check.hh"

using namespace std;

void parse_folder(string input_folder,
                  unsigned &number_of_cells,
                  unsigned &number_of_groups,
                  unsigned &number_of_ordinates,
                  unsigned &max_iterations,
                  double &tolerance,
                  vector<double> &psi,
                  vector<double> &internal_source,
                  vector<double> &boundary_sources,
                  vector<string> &boundary_conditions,
                  vector<double> &cell_length,
                  vector<double> &sigma_t,
                  vector<double> &sigma_s,
                  vector<double> &nu_sigma_f,
                  vector<double> &chi,
                  unsigned number_of_materials)
{
    unsigned number_of_nodes = 2;
    
    parse_data(input_folder + "/number_of_cells", number_of_cells);
    parse_data(input_folder + "/number_of_groups", number_of_groups);
    parse_data(input_folder + "/number_of_ordinates", number_of_ordinates);
    parse_data(input_folder + "/max_iterations", max_iterations);
    parse_data(input_folder + "/tolerance", tolerance);
    parse_data(input_folder + "/psi", psi);
    parse_data(input_folder + "/internal_source", internal_source);
    parse_data(input_folder + "/boundary_sources", boundary_sources);
    parse_data(input_folder + "/boundary_conditions", boundary_conditions);
    parse_data(input_folder + "/cell_length", cell_length);
    parse_data(input_folder + "/sigma_t", sigma_t);
    parse_data(input_folder + "/sigma_s", sigma_s);
    parse_data(input_folder + "/nu_sigma_f", nu_sigma_f);
    parse_data(input_folder + "/chi", chi);

    unsigned materials_nodes_cells_groups_ordinates = number_of_materials * number_of_nodes * number_of_cells * number_of_groups * number_of_ordinates;
    unsigned cells_groups = number_of_cells * number_of_groups;
    unsigned groups_ordinates = number_of_groups * number_of_ordinates;
    unsigned materials_cells_groups_groups = number_of_materials * number_of_cells * number_of_groups * number_of_groups;
    unsigned number_of_edges = 2;
    unsigned materials_cells_groups = number_of_materials * number_of_cells * number_of_groups;
    
    check(materials_nodes_cells_groups_ordinates, psi, "psi");
    check(cells_groups, internal_source, "internal_source");
    check(groups_ordinates, boundary_sources, "boundary_source");
    check(number_of_edges, boundary_conditions, "boundary_conditions");
    check(number_of_cells, cell_length, "cell_length");
    check(materials_cells_groups, sigma_t, "sigma_t");
    check(materials_cells_groups_groups, sigma_s, "sigma_s");
    check(materials_cells_groups, nu_sigma_f, "nu_sigma_f");
    check(materials_cells_groups, chi, "chi");
}

void parse_folder(string input_folder,
                  unsigned &number_of_cells,
                  unsigned &number_of_groups,
                  unsigned &number_of_ordinates,
                  unsigned &number_of_materials,
                  unsigned &number_of_benchmarks,
                  unsigned &max_iterations,
                  double &tolerance,
                  vector<double> &psi,
                  vector<double> &psi_total,
                  vector<double> &internal_source,
                  vector<double> &boundary_sources,
                  vector<string> &boundary_conditions,
                  vector<double> &cell_length,
                  vector<double> &sigma_t,
                  vector<double> &sigma_s,
                  vector<double> &nu_sigma_f,
                  vector<double> &chi,
                  vector<double> &chord_length,
                  vector<string> &methods)
{
    unsigned number_of_nodes = 2;
    
    parse_data(input_folder + "/chord_length", chord_length);
    parse_data(input_folder + "/number_of_materials", number_of_materials);
    parse_data(input_folder + "/psi_total", psi_total);
    parse_data(input_folder + "/number_of_benchmarks", number_of_benchmarks);
    parse_data(input_folder + "/methods", methods);
    
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
                 chi,
                 number_of_materials);

    unsigned cells_nodes_groups_ordinates = number_of_cells * number_of_nodes * number_of_groups * number_of_ordinates;
    
    check(number_of_materials, chord_length, "chord_length");
    check(cells_nodes_groups_ordinates, psi_total, "psi_total");
}

void parse_data(string path,
                vector<double> &data)
{
    data.resize(0);

    ifstream input_file(path.c_str());
    
    if (input_file.is_open())
    {
        double temp;
        while (input_file >> temp)
        {
            data.push_back(temp);
        }
    }
}

void parse_data(string path,
                vector<string> &data)
{
    data.resize(0);
    
    ifstream input_file(path.c_str());
    
    if (input_file.is_open())
    {
        string temp;
        while (input_file >> temp)
        {
            data.push_back(temp);
        }
    }
}

void parse_data(string path,
                unsigned &data)
{
    ifstream input_file(path.c_str());
    
    if (input_file.is_open())
    {
        input_file >> data;     
    }
}
    
void parse_data(string path,
                double &data)
{
    ifstream input_file(path.c_str());
    
    if (input_file.is_open())
    {
        input_file >> data;     
    }
}

void parse_data(string path,
                string &data)
{
    ifstream input_file(path.c_str());
    
    if (input_file.is_open())
    {
        input_file >> data;     
    }
}


void output_data(string path,
                 vector<double> &data)
{
    ofstream output_file(path.c_str());

    if (output_file.is_open())
    {
        for (unsigned i = 0; i < data.size(); ++i)
        {
            output_file << data[i] << endl;
        }
    }
}

void output_data(string path,
                 vector<unsigned> &data)
{
    ofstream output_file(path.c_str());

    if (output_file.is_open())
    {
        for (unsigned i = 0; i < data.size(); ++i)
        {
            output_file << data[i] << endl;
        }
    }
}

void output_data(string path,
                 double &data)
{
    ofstream output_file(path.c_str());

    if (output_file.is_open())
    {
        output_file << data;
    }
}

void output_data(string path,
                 unsigned &data)
{
    ofstream output_file(path.c_str());

    if (output_file.is_open())
    {
        output_file << data;
    }
}
