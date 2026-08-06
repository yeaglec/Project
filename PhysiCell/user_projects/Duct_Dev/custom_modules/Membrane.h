#ifndef __MEMBRANE_H__
#define __MEMBRANE_H__

#include "custom.h"
#include <vector>
#include <utility>

using namespace PhysiCell;

// Active Mechanics
void cell_interactions_cc(Cell* pCell, Phenotype& phenotype, double dt);
std::pair<double,double> basement_membrane_interactions_cc(Cell* pCell, double de, double px, double py);
void update_basement_membrane_deformation(double dt);

// Structural Forces
void BM_Smoothing(std::vector<std::pair<double,double>>& node_forces, double Fx_BM, double Fy_BM, int best_k, double best_px, double best_py, double best_t);
void membrane_strain_lin(std::vector<std::pair<double,double>>& node_forces);
void membrane_strain_exp(std::vector<std::pair<double,double>>& node_forces);
void membrane_restoring_force_lin(std::vector<std::pair<double,double>>& node_forces);
void membrane_restoring_force_exp(std::vector<std::pair<double,double>>& node_forces);
void membrane_bending_stiffness(std::vector<std::pair<double,double>>& node_forces);
void add_membrane_nodes();

// Functions decoupling from Level Set Method
double distance_to_membrane(Cell* pCell, Phenotype& phenotype, double dt);
void rebuild_signed_distance_field();
void initialize_level_set_duct(std::vector<std::vector<double>> boundary_membrane_pts);
std::pair<double,double> voxel_indices(Cell* pCell);


#endif