#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include "custom.h"
#include <vector>
#include <tuple>
#include <string>

// Boundary Mathematics
bool is_inside(double x, double y, const std::vector<std::vector<double>>& BM_pts);
std::tuple<double, double, double, double, double> project_point_onto_boundary(double x, double y);

// Initialization Shapes
std::vector<std::vector<double>> generate_boundary_shape(double a, double b, double amp, int freq);
void generate_boundary_cells(double a, double b, double amp, int freq, std::string type, double dis, int num_cells);
std::vector<std::vector<double>> generate_circle_boundary();
void generate_circle_cells();
void clustered_cell(double a, double b, double amp, int freq, int num_points, double center_x, double center_y, double radius, std::string type);

#endif