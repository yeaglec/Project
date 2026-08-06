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
#include "./Membrane.h"
#include "./Geometry.h"
#include "./Utils.h"
#include "./Tests.h"
#include <cmath>
#include <cfloat>
#include <array>
#include <iostream>

std::vector<std::vector<double>> boundary_membrane_pts;
std::vector<double> initial_edge_length; 
std::vector<std::vector<double>> initial_node_positions; 

std::vector<std::vector<double>> level_set_phi; 
double ls_dx = 0.0, ls_dy = 0.0;
double ls_xmin = 0.0, ls_ymin = 0.0;

int BM_Fx_idx, BM_Fy_idx, BM_k_idx, BM_px_idx, BM_py_idx, BM_t_idx;
void (*test_perb)(std::vector<std::pair<double,double>>&, double) = nullptr;


// ______________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// Testing functions for membrane interactions based on cell centers (not voxel-based)
// ________________________________________________________________________________________________________________________



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
	if (parameters.ints.find_index("random_seed") != -1) // TODO CLEAN: One liner
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

    cell_defaults.custom_data.add_variable( "BM_Fx", "micron/min", 0.0 );
    cell_defaults.custom_data.add_variable( "BM_Fy", "micron/min", 0.0 );
    cell_defaults.custom_data.add_variable( "BM_k", "dimensionless", -1.0 );
    cell_defaults.custom_data.add_variable( "BM_px", "micron", 0.0 );
    cell_defaults.custom_data.add_variable( "BM_py", "micron", 0.0 );
    cell_defaults.custom_data.add_variable( "BM_t", "dimensionless", 0.0 );

    BM_Fx_idx = cell_defaults.custom_data.find_variable_index("BM_Fx");
    BM_Fy_idx = cell_defaults.custom_data.find_variable_index("BM_Fy");
    BM_k_idx = cell_defaults.custom_data.find_variable_index("BM_k");
    BM_px_idx = cell_defaults.custom_data.find_variable_index("BM_px");
    BM_py_idx = cell_defaults.custom_data.find_variable_index("BM_py");
    BM_t_idx = cell_defaults.custom_data.find_variable_index("BM_t");
 
	
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

	std::vector<std::string> bm_targets = {"Epithelial","CAF","BM", "Tumor"};
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
		cell_defaults.functions.cell_division_function = parallel_cell_division;   // TODO CLEAN: Make one liner
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
	//_______________________________________________________________________________________________________________________
    // Initialization

    // Placing cells to test the basement membrane deformation
	// Example code for generating a arbitrary boundary

	int num_points = parameters.ints("membrane_num_points");
	double a = 300.0, b = 250.0;
	double amp = 0.1;              // Amplitude of deformation
	int freq = 4;  
	int num_ep = parameters.ints("number_EP_cells");

	boundary_membrane_pts = generate_boundary_shape(a, b, amp, freq);
	double ep_dis = parameters.doubles("ep_displacement");
	generate_boundary_cells(a, b, amp, freq, "Epithelial", ep_dis, num_ep);

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

	generate_boundary_cells(a, b, amp, freq, "CAF", -5, num_caf);
  	// Cell_Definition* pTumorDef = cell_definitions_by_name["CAF"];
	// Cell* Caf = create_cell( *pTumorDef );
	// Caf->assign_position( { 10,10, 0.0 } );

    // _____________ TESTING  Triangle Membrane Elasticity and Restoring Force __________________

	// int num_caf = parameters.ints("number_CAF_cells");
	// Cell_Definition* Caf_def = cell_definitions_by_index[2];  // Need at least 1 cell or sim gets mad
	// Cell* Caf = create_cell( *Caf_def );
	// Caf->assign_position( { 225,200, 0.0 } );

    // // Put Test Functions Here
    // Test_Ring();
    // test_perb = nullptr;  //nullptr if not testing

	
	// _____________ TESTING Membrane Elasticity and Restoring Force __________________

	// boundary_membrane_pts.push_back({0.0, 50, 0.0}); 
	// boundary_membrane_pts.push_back({0.0, -50, 0.0});

	// int Np = (int)boundary_membrane_pts.size(); // Should be 2
	// std::cout << "Initial Boundary Points: " << Np << std::endl;
	// initial_edge_length.resize(Np);
	// initial_node_positions.resize(Np);

	// int num_caf = parameters.ints("number_CAF_cells");
	// Cell_Definition* Caf_def = cell_definitions_by_index[2];
	// Cell* Caf = create_cell( *Caf_def );
	// Caf->assign_position( { 225,200, 0.0 } );

    // _____________ TESTING Added Membrane Points and Bending Stiffness __________________

	// int num_caf = parameters.ints("number_CAF_cells");
	// Cell_Definition* Caf_def = cell_definitions_by_index[2];  // Need at least 1 cell or sim gets mad
	// Cell* Caf = create_cell( *Caf_def );
	// Caf->assign_position( { 225,200, 0.0 } );

    // // Put Test Functions Here
    // Test_Remesh();
    // test_perb = Test_Remesh_Pert;  //nullptr if not testing

	
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