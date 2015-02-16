#ifndef Parser_hh
#define Parser_hh

#include <vector>
#include <string>

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
                  unsigned number_of_materials = 1);     

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
                  vector<string> &methods);

void parse_data(string path,
                vector<double> &data);

void parse_data(string path,
                vector<string> &data);

void parse_data(string path,
                unsigned &data);

void parse_data(string path,
                double &data);

void parse_data(string path,
                string &data);

void output_data(string path,
                 vector<double> &data);

void output_data(string path,
                 vector<unsigned> &data);

void output_data(string path,
                 double &data);

void output_data(string path,
                 unsigned &data);

#endif
