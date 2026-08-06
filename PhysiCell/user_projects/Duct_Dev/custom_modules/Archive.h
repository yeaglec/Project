#ifndef __ARCHIVE_H__
#define __ARCHIVE_H__

#include "custom.h"
#include <vector>
#include <utility>

using namespace PhysiCell;

// Level Set Math Helpers
std::pair<double, double> level_set_gradient(int i, int j, const std::vector<std::vector<double>>& phi, double dx, double dy);
std::pair<double,double> level_set_normalize(const std::pair<double, double>& gradient);

// Shelved LSM Mechanics
void cell_interactions_LSM(Cell* pCell, Phenotype& phenotype, double dt);
std::pair<double,double> basement_membrane_interactions_LSM(Cell* pCell);
void update_basement_membrane_deformation_LSM(double dt); // Renamed in header to avoid conflict

#endif