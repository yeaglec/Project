#ifndef __BM_utils_h__
#define __BM_utils_h__

#define _USE_MATH_DEFINES
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include <vector>
#include <string>
#include <tuple>
#include <array>

// Helper function for visualization of boundary
void boundary_to_csv(
    const std::vector<std::vector<double>>& boundary_pts,
    const std::string& filename
);

std::pair<double, double> level_set_gradient(int i, int j, 
	const std::vector<std::vector<double>>& phi,
	double dx, double dy);

std::pair<double,double> level_set_normalize(const std::pair<double, double>& gradient);
std::pair<double,double> voxel_indices(Cell* pCell);

bool is_inside(double x, double y, const std::vector<std::vector<double>>& BM_pts);

std::tuple<double, double, double, double, double> project_point_onto_boundary(double x, double y);
std::vector<std::vector<double>> generate_boundary_shape(double a, double b, double amp, int freq);
std::vector<std::vector<double>> generate_circle_boundary();
void generate_circle_cells();
void generate_boundary_cells(double a, double b, double amp, int freq, std::string type, double dis, int num_cells);

void parallel_cell_division( Cell* parent, Cell* child );







#endif