/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#define _USE_MATH_DEFINES
#include <cmath>
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 
#include "fast_marching_method.hpp"

using namespace BioFVM; 
using namespace PhysiCell;

// setup functions to help us along 

void create_cell_types( void );
void setup_tissue( void ); 

// set up the BioFVM microenvironment 
void setup_microenvironment( void ); 

// custom pathology coloring function 

std::vector<std::string> my_coloring_function( Cell* );

// custom functions can go here 

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt );
void custom_function( Cell* pCell, Phenotype& phenotype , double dt );

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt ); 

// functions for the duct

double distance_to_membrane(Cell* pCell, Phenotype& phenotype, double dt);
void basement_membrane_interaction(Cell* pCell, Phenotype& phenotype, double dt);

// _____________________________________
// Defining the FMM

// In custom.cpp


void compute_sdf_from_initial_boundary( std::vector<double>& initial_phi )
{
    // 1. --- Define Grid Properties ---
    const std::array<int, 3> grid_size = {
        (int)microenvironment.mesh.x_coordinates.size(),
        (int)microenvironment.mesh.y_coordinates.size(),
        (int)microenvironment.mesh.z_coordinates.size()
    };
    const std::array<double, 3> grid_spacing = {
        microenvironment.mesh.dx,
        microenvironment.mesh.dy,
        microenvironment.mesh.dz
    };

    // 2. --- Define the Speed ---
    // For a distance function, the speed is uniformly 1.0
    const std::vector<double> speed( microenvironment.number_of_voxels(), 1.0 );

    // 3. --- Identify Known Points (The "Seeds") ---
    std::vector<std::array<int, 3>> known_indices;
    std::vector<double> known_distances;

    // We need to seed the FMM with points on the boundary (dist=0),
    // points just inside (dist<0), and points just outside (dist>0).
    for (int k = 0; k < grid_size[0]; ++k) {
        for (int j = 0; j < grid_size[1]; ++j) {
            for (int i = 0; i < grid_size[2]; ++i) {
                int n = microenvironment.voxel_index(i, j, k);
                
                // If this voxel was marked as being on the boundary by our rasterizer
                if ( initial_phi[n] == 0.0 ) {
                    known_indices.push_back({i, j, k});
                    known_distances.push_back(0.0);
                }
            }
        }
    }
    
    // This is a simplified seeding for the sign. A more robust method would check
    // neighbors. We assume the origin is inside the duct.
    // We check the sign of the center-most voxel relative to the shape's center.
    int center_voxel_index = microenvironment.nearest_voxel_index({0,0,0});
    // A simple test: if the center voxel is far from the boundary, it's inside.
    if( initial_phi[center_voxel_index] > 0 ) {
        // The sign needs to be flipped. Let's make inside negative.
        for( int n=0; n < microenvironment.number_of_voxels(); n++ )
        {
            // A simple heuristic to define inside/outside based on rasterization
            if( initial_phi[n] > 0 ) initial_phi[n] = -1.0; // Mark as inside
            else initial_phi[n] = 1.0; // Mark as outside
        }
    } else {
        // The sign is correct.
        for( int n=0; n < microenvironment.number_of_voxels(); n++ )
        {
            if( initial_phi[n] > 0 ) initial_phi[n] = 1.0;
            else initial_phi[n] = -1.0;
        }
    }
    // Now, initial_phi has -1 for inside, +1 for outside, and 0 on boundary.
    // The FMM library will use the signs of the known_distances to propagate.
    // We can refine the known distances for better accuracy.
    for( size_t i=0; i < known_indices.size() ; ++i )
    {
        int n = microenvironment.voxel_index(known_indices[i], known_indices[i][2], known_indices[i][1]);
        // Check a neighbor to determine sign
        int neighbor_idx = microenvironment.voxel_index(known_indices[i]+1, known_indices[i][2], known_indices[i][1]);
        if( neighbor_idx >= 0 && neighbor_idx < initial_phi.size() && initial_phi[neighbor_idx] < 0)
        {
            // If neighbor is inside, this boundary point should propagate negative values inward
            // But the library handles this based on the initial values. We can just provide 0.
        }
    }
    // The library can handle positive and negative arrival times (distances).
    // We provide the boundary points with distance 0. The library will compute
    // distances, but we need to apply the sign afterwards.

    // 4. --- Run the Fast Marching Method ---
    // The library calculates unsigned distance. We will calculate it twice.
    // Once for the inside, once for the outside.

    // Run 1: Distance inside the duct (negative values)
    std::vector<std::array<int, 3>> inside_seeds;
    std::vector<double> inside_distances;
    // Find boundary voxels with an inside neighbor
    for (const auto& idx_array : known_indices) {
        // Check neighbors of this boundary voxel
        // (A full implementation would check all 6 neighbors)
        int neighbor_i = idx_array + 1;
        if (neighbor_i < grid_size) {
            int neighbor_n = microenvironment.voxel_index(neighbor_i, idx_array[2], idx_array[1]);
            if (initial_phi[neighbor_n] < 0) { // If neighbor is inside
                inside_seeds.push_back(idx_array);
                inside_distances.push_back(0.0);
                break; // Move to next boundary voxel
            }
        }
    }
    auto inside_dist_field = fmm::ArrivalTime(grid_size, grid_spacing, inside_seeds, inside_distances, speed);

    // Run 2: Distance outside the duct (positive values)
    std::vector<std::array<int, 3>> outside_seeds;
    std::vector<double> outside_distances;
    // Find boundary voxels with an outside neighbor
    for (const auto& idx_array : known_indices) {
        int neighbor_i = idx_array + 1;
        if (neighbor_i < grid_size) {
            int neighbor_n = microenvironment.voxel_index(neighbor_i, idx_array[2], idx_array[1]);
            if (initial_phi[neighbor_n] > 0) { // If neighbor is outside
                outside_seeds.push_back(idx_array);
                outside_distances.push_back(0.0);
                break; 
            }
        }
    }
    auto outside_dist_field = fmm::ArrivalTime(grid_size, grid_spacing, outside_seeds, outside_distances, speed);

    // 5. --- Combine and Assign to the Final SDF ---
    for( int n=0; n < microenvironment.number_of_voxels(); n++ )
    {
        if( initial_phi[n] <= 0 ) // If inside or on the boundary
        {
            microenvironment.level_set_phi[n] = -inside_dist_field[n];
        }
        else // If outside
        {
            microenvironment.level_set_phi[n] = outside_dist_field[n];
        }
    }
}
