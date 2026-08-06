#ifndef __UTILS_H__
#define __UTILS_H__

#include "custom.h"
#include <vector>
#include <string>

void boundary_to_csv(const std::vector<std::vector<double>>& boundary_pts, const std::string& filename);
void parallel_cell_division(Cell* parent, Cell* child);

#endif