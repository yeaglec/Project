#ifndef __BM_tests_h__
#define __BM_tests_h__

#define _USE_MATH_DEFINES
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include <vector>
#include <string>
#include <tuple>
#include <array>

void Test_Line();
void Test_Line_Pert(std::vector<std::pair<double,double>>& node_forces, double current_time);

void Test_Tri();
void Test_Tri_Pert(std::vector<std::pair<double,double>>& node_forces, double current_time);

void Test_Square();
void Test_Square_Pert(std::vector<std::pair<double,double>>& node_forces, double current_time);

void Test_Ring();
void Test_Ring_Pert(std::vector<std::pair<double,double>>& node_forces, double current_time);

#endif