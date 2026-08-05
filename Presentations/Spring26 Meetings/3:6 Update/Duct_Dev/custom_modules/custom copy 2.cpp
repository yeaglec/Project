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

#include "./custom.h"
#include "./Models.h"
#include "./Utils.h"
#include <cmath>
#include <cfloat>
#include <array>

std::vector<std::vector<double>> boundary_membrane_pts;
std::vector<double> initial_edge_length; 


// ______________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// Testing functions for membrane interactions based on cell centers (not voxel-based)
// ________________________________________________________________________________________________________________________

// ################ Function that computes the BM-to-Cell force ####################
std::pair<double,double> basement_membrane_interactions_cc(Cell* pCell)
{
    double cell_x = pCell->position[0];
    double cell_y = pCell->position[1];

    auto [dist, px, py, k, t] = project_point_onto_boundary(cell_x, cell_y);
    bool inside = is_inside(cell_x, cell_y, boundary_membrane_pts);

    double d = inside ? -dist : dist;
    double R = pCell->phenotype.geometry.radius;
    double de = d - (d < 0 ? -R : R); // Effective distance

	// if (de > 0) {

	// 	return {0.0, 0.0}; // Remove deformation if cells pass through boundary
    // }

    double BM_deadzone = parameters.doubles("membrane_deadzone");
    if (fabs(de) < BM_deadzone) return {0.0, 0.0};

    double strength = parameters.doubles("membrane_spring_constant"); 
    double mag = strength * fabs(de);  // Hooke's law, F = kx

    double nx = px - cell_x; // Normal vector from cell to boundary
    double ny = py - cell_y;
    double norm = sqrt(nx * nx + ny * ny);
	
    if (norm > 0)
    {
        nx /= norm;
        ny /= norm;
    }

	// mag = parameters.doubles("test_velocity_magnitude");
    double Fx =  mag * nx;  // Force components: Hookean force applied in the normal direction
    double Fy =  mag * ny;

    return {Fx, Fy};
}

// ################### Function that implements Cell-to-BM force ####################
void cell_interactions_cc(Cell* pCell,
                                   Phenotype& phenotype,
                                   double dt)
{
    double cell_x = pCell->position[0];
    double cell_y = pCell->position[1];

    auto [dist, px, py, k, t] = project_point_onto_boundary(cell_x, cell_y);

    bool inside = is_inside(cell_x, cell_y, boundary_membrane_pts);
    double d = inside ? -dist : dist;

    double R = pCell->phenotype.geometry.radius;
    double de = d - (d < 0 ? -R : R);

    // --- quick type check ---
    std::string cell_name = cell_definitions_by_index[pCell->type]->name;
    bool isCAF = (cell_name == "CAF");
    bool isEP  = (cell_name == "Epithelial");


    if (isEP)
    {
        // Implement replusion here (original EP behavior)
        if (de > 0)
        {
            double cell_deadzone = parameters.doubles("cell_deadzone");
            double displacement_needed = de + cell_deadzone; // Make cells only move to deadzone

            double nx = cell_x - px; // away from boundary
            double ny = cell_y - py;
            double norm = sqrt(nx * nx + ny * ny);
            if (norm > 1e-16)
            {
                nx /= norm;
                ny /= norm;
            }

            // Calculate the mag of the corrective velocity
            double correction_rate = parameters.doubles("membrane_correction_rate");
            double mag = correction_rate * displacement_needed;

            pCell->velocity[0] += mag * nx;
            pCell->velocity[1] += mag * ny;
        }
        else
        {
            double L = parameters.doubles("membrane_interaction_length");
            if (fabs(de) >= L) return;

            double cell_deadzone = parameters.doubles("cell_deadzone");
            if (fabs(de) < cell_deadzone) return;

            double strength = parameters.doubles("membrane_adhesion_strength");
            double mag = strength * fabs(de);

            double nx = px - cell_x;
            double ny = py - cell_y;
            double norm = sqrt(nx * nx + ny * ny);
            if (norm > 0)
            {
                nx /= norm;
                ny /= norm;
            }

            pCell->velocity[0] += mag * nx;
            pCell->velocity[1] += mag * ny;
        }

        return;
    }

    // --- CAF: treat de < 0 as penetration (they were outside originally) ---
    if (isCAF)
    {
        if (de < 0)
        {
            double cell_deadzone = parameters.doubles("cell_deadzone");
            double displacement_needed = fabs(de) + cell_deadzone; // move out to deadzone

            double nx = cell_x - px; // away from boundary
            double ny = cell_y - py;
            double norm = sqrt(nx * nx + ny * ny);
            if (norm > 1e-16)
            {
                nx /= norm;
                ny /= norm;
            }

            double correction_rate = parameters.doubles("membrane_correction_rate");
            double mag = correction_rate * displacement_needed;

            pCell->velocity[0] += mag * nx;
            pCell->velocity[1] += mag * ny;
        }
        else
        {
            // same adhesion behavior as EPs when not "inside"
            double L = parameters.doubles("membrane_interaction_length");
            if (fabs(de) >= L) return;

            double cell_deadzone = parameters.doubles("cell_deadzone");
            if (fabs(de) < cell_deadzone) return;

            double strength = parameters.doubles("membrane_adhesion_strength");
            double mag = strength * fabs(de);

            double nx = px - cell_x;
            double ny = py - cell_y;
            double norm = sqrt(nx * nx + ny * ny);
            if (norm > 0)
            {
                nx /= norm;
                ny /= norm;
            }

            pCell->velocity[0] += mag * nx;
            pCell->velocity[1] += mag * ny;
        }

        return;
    }

    // For other cell types: do nothing (or keep behavior by removing this return)
    return;
}



void update_basement_membrane_deformation2(double dt)
{
    int Np = (int)boundary_membrane_pts.size();
    std::vector<std::pair<double,double>> node_forces(Np, {0.0, 0.0});

    double sigma = parameters.doubles("membrane_force_smoothing_sigma"); // Smoothing parameter for Gaussian distribution
	double cutoff = 3.0 * sigma;         // ignore nodes beyond 3 sigma 

    // Iterate over all cells and compute membrane forces distributed by Gaussian
    for (Cell* pCell : *all_cells) {

        auto cell_force = basement_membrane_interactions_cc(pCell);  // Cell-to-BM force
        double Fx_cell = cell_force.first;
        double Fy_cell = cell_force.second;

        if (Fx_cell == 0.0 && Fy_cell == 0.0) continue;

        double Fx_BM = -Fx_cell;   // The velocity is being applied to the voxel center for all the boundary nodes (TODO: FIx)
        double Fy_BM = -Fy_cell;

        // Computing the projection of a cell onto the BM TODO: Make helper for proj code
        double cell_x = pCell->position[0];
        double cell_y = pCell->position[1];

        auto [best_dist, best_px, best_py, best_k_d, best_t] = project_point_onto_boundary(cell_x, cell_y);
		int best_k = static_cast<int>(std::round(best_k_d));

		// ###############################################
        // Build Gaussian weights centered at projection point and distribute force (TODO make helper function for this)
        double sumW = 0.0;
    
        std::vector<int> nearby_indices;   // Try storing per-node weights only for nodes within cutoff
        std::vector<double> nearby_weights;
        nearby_indices.reserve(16);
        nearby_weights.reserve(16);

		if (sigma <= 0 || std::isnan(sigma)){

			// Revert to the original 2-node split based on projection t.
            int n1 = best_k;
            int n2 = (best_k + 1) % Np;
            double t_clamped = std::max(0.0, std::min(best_t, 1.0));
            node_forces[n1].first  += (1.0 - t_clamped) * Fx_BM;
            node_forces[n1].second += (1.0 - t_clamped) * Fy_BM;
            node_forces[n2].first  += (t_clamped) * Fx_BM;
            node_forces[n2].second += (t_clamped) * Fy_BM;	
			
		}

		else{
			// Checking if neighboring nodes are close enough to the projection point
			for (int i = 0; i < Np; i++) {
				double dx = boundary_membrane_pts[i][0] - best_px;
				double dy = boundary_membrane_pts[i][1] - best_py;
				double dist2 = dx*dx + dy*dy;
				if (dist2 > cutoff*cutoff ) continue; // ignore distant nodes
				double w = exp(-dist2 / (2.0 * sigma * sigma));   // -d^2 / (2 * sigma^2)
				// if (w <= 0.0) continue;

				nearby_indices.push_back(i); // Store index & weight of nearby node 
				nearby_weights.push_back(w);
				sumW += w;
			}

            // normalize weights and distribute
            for (size_t idx = 0; idx < nearby_indices.size(); ++idx) {
                int node_i = nearby_indices[idx];
                double w = nearby_weights[idx];
                double frac = w / sumW;
                node_forces[node_i].first  += Fx_BM * frac;
                node_forces[node_i].second += Fy_BM * frac;
            }
        }
    } 

	// ##############################################
	// Nonlinear/Fung membrane elasticity 

	// Material parameters given in Li et al. 
	double G = 82.30;    // low-strain elastic modulus 
	double b = 2.10;    // High-strain stiffening parameter

	if (parameters.doubles("elastic_BM")==1.0){ 
		for (int i = 0; i < Np; ++i) {
			int j = (i + 1) % Np;

			double dx = boundary_membrane_pts[j][0] - boundary_membrane_pts[i][0];
			double dy = boundary_membrane_pts[j][1] - boundary_membrane_pts[i][1];

			std::cout << "__dx: " << dx << " dy: " << dy << std::endl;

			double length = std::sqrt(dx*dx + dy*dy);

			std::cout << "LENGTH: " << length << std::endl;
			double rest   = initial_edge_length[i];  // stored initial length
			std::cout << "__rest: " << rest << std::endl;


			if (length <= 1e-12 || rest <= 1e-12){

				continue;
			}

			double fspring = 0.0;
			double lambda = length / rest;
			double inv_lambda = 1.0 / lambda;

			// compute exponent argument

			double high_strain = b * (lambda*lambda + 2.0*inv_lambda - 3.0);
			if (high_strain > 100.0) high_strain = 100.0;
			if (high_strain < -100.0) high_strain = -100.0;

			double Fung_sigma = G * (lambda*lambda - inv_lambda) * std::exp(high_strain); // Fung's Model 
			fspring = - Fung_sigma * rest;

			// direction components
			double fx = fspring * (dx / length);
			double fy = fspring * (dy / length);

			// add to nodal forces (i is pulled toward j and vice versa)
			node_forces[i].first  -= fx;
			node_forces[i].second -= fy;
			node_forces[j].first  += fx;
			node_forces[j].second += fy;
		}
	}

	// Update node positions
    for (int i = 0; i < Np; ++i) {
        boundary_membrane_pts[i][0] += node_forces[i].first  * dt;
        boundary_membrane_pts[i][1] += node_forces[i].second * dt;
    }

    // Rebuild signed distance field after modifications
    rebuild_signed_distance_field();
}

void clustered_cell(double a, double b, double amp, int freq, int num_points, double center_x, double center_y, double radius, std::string type)
{
	Cell_Definition* pBM_def = cell_definitions_by_name[type];   //cell_definitions_by_name[ type_name ] = pCD; 
	std::cout << "Creating clustered BM cells..." << std::endl;
	std::vector<std::vector<double>> pts;
	pts.reserve(num_points);


    for (int j = 0; j < num_points; j++) {
        // Place cells randomly within a disc for a more real look
        double r = radius * sqrt(UniformRandom()); // sqrt for uniform area distribution
        double theta = 2.0 * M_PI * UniformRandom();
        
        double cx = center_x + r * cos(theta);
        double cy = center_y + r * sin(theta);
        
        Cell* pC2 = create_cell(*pBM_def);
        pC2->assign_position({cx, cy, 0.0}); 
    }

}

// ________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// Main PhysiCell Functions
// ________________________________________________________________________________________________________________________

void custom_rule( Cell* pCell, Phenotype& phenotype, double dt )
{	

	return; // This is a custom rule that can be used to implement any custom behavior for the cell, such a proliferation.
}

void create_cell_types( void )
{
	// set the random seed 
	if (parameters.ints.find_index("random_seed") != -1)
	{
		SeedRandom(parameters.ints("random_seed"));
	}
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 

	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.add_cell_basement_membrane_interactions = cell_interactions_cc; 
	cell_defaults.functions.calculate_distance_to_membrane = distance_to_membrane; 

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
 
	
	// std::cout << "Cell type defaults set. (Test)" << std::endl;
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		

	build_cell_definitions_maps(); 

	// Ensure BM helper functions are assigned for types that should have BM interactions

	std::vector<std::string> bm_targets = {"Epithelial","CAF","BM"};
	for (auto &name : bm_targets) {
		auto it = cell_definitions_by_name.find(name);
		if (it == cell_definitions_by_name.end()) continue;
		Cell_Definition* pDef = it->second;
		pDef->functions.add_cell_basement_membrane_interactions = cell_interactions_cc;
		pDef->functions.calculate_distance_to_membrane = distance_to_membrane;
	}


	
	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/*
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_rule; 
	cell_defaults.functions.contact_function = contact_function; 

	// Adding function for parallel division
	Cell_Definition* pCD = cell_definitions_by_name["Epithelial"];
	if (pCD)
	{	
		cell_defaults.functions.cell_division_function = parallel_cell_division;
	}

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();

	//_______________________________________________________________________________________________________________________
	// Placing cells to test the basement membrane deformation

	// Example code for generating a arbitrary boundary

	int num_points = parameters.ints("membrane_num_points");
	double a = 300.0, b = 250.0;
	double amp = 0.1;              // Amplitude of deformation
	int freq = 4;  
	int num_ep = parameters.ints("number_EP_cells");

	boundary_membrane_pts = generate_boundary_shape(a, b, amp, freq);
	generate_boundary_cells(a, b, amp, freq);

	int num_caf = parameters.ints("number_CAF_cells");
	// Cell_Definition* Caf_def = cell_definitions_by_index[2];
	// Cell* Caf = create_cell( *Caf_def );
	// Caf->assign_position( { 225,200, 0.0 } );

	double CAFx = parameters.doubles("CAFx");
	double CAFy = parameters.doubles("CAFy");

	double EPx = parameters.doubles("EPx");
	double EPy = parameters.doubles("EPy");

	double CAF_rad = parameters.doubles("CAF_rad");
	double EP_rad = parameters.doubles("EP_rad");

	// clustered_cell(a, b, amp, freq, num_caf, CAFx, CAFy, CAF_rad, "CAF");
	


	// ##################

	// Example Code for generating a circle boundary

    // int num_points = parameters.ints("membrane_num_points");
	// int num_ep = parameters.ints("number_EP_cells");
	// double circle_radius = parameters.doubles("membrane_circle_radius");

	// boundary_membrane_pts = generate_circle_boundary(circle_radius, num_points);
	// generate_circle_cells(circle_radius, num_ep);


	//_______________________________________________________________________________________________________________________

	for (auto pCell : *all_cells){

		// Stash the volume growth rates in custom data
		pCell->custom_data["default_cyto_rate"] = pCell->phenotype.volume.cytoplasmic_biomass_change_rate;
		pCell->custom_data["default_nuclear_rate"] = pCell->phenotype.volume.nuclear_biomass_change_rate;
		pCell->custom_data["default_fluid_rate"] = pCell->phenotype.volume.fluid_change_rate;

	}
	
	
    initialize_level_set_duct(boundary_membrane_pts);

	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 