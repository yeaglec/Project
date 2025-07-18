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
#include <cmath>
#include <cfloat>
#include <array>

std::vector<std::vector<double>> boundary_membrane_pts;

// Helper function to save boundary points to a CSV file
void boundary_to_csv( std::vector<std::vector<double>> const& boundary_pts, 
					  std::string const& filename )
{	
	int index = (int)PhysiCell_globals.full_output_index;
	std::cout << "Saving boundary points to " << filename << " at time " << index << "\n" << std::endl;


	std:: ofstream file(filename);
	if (!file.is_open()) {
		std::cerr << "Error opening file: " << filename << std::endl;
		return;
	}

	for (auto& pt : boundary_pts) {
		file << pt[0] << "," << pt[1] << "\n";
	}
	file.close();

	std::cout << "Boundary points saved to " << filename << " successfully" << std::endl;
}

// Helper function to compute the gradient (TODO: Verify if this is correct, I have doubts)

std::pair<double, double> level_set_gradient(int i, int j, 
										 const std::vector<std::vector<double>>& phi,
										 double dx, double dy)
{
	int Nx = (int)phi.size();
	int Ny = (int)phi[0].size();
	
	double phi_x, phi_y;

	// x-deriv
	if(i>0 && i < Nx-1) {
	phi_x = (phi[i+1][j] - phi[i-1][j])/(2*dx);}

	else if(i==0) {
	phi_x = (phi[i+1][j] - phi[i][j])/(dx);
	}

	else { // i==Nx-1
	phi_x = (phi[i][j] - phi[i-1][j])/(dx);
	}

	// y-deriv
	if(j>0 && j < Ny-1) {
	phi_y = (phi[i][j+1] - phi[i][j-1])/(2*dy);}

	else if(j==0) {
	phi_y = (phi[i][j+1] - phi[i][j])/(dy);}

	else { // j==Ny-1
	phi_y = (phi[i][j] - phi[i][j-1])/(dy);
	}
	
	return { phi_x, phi_y };
}

std::pair<double,double> level_set_normalize(const std::pair<double, double>& gradient)
{
	std::vector<double> normal = { gradient.first, gradient.second };
	double norm = sqrt(gradient.first*gradient.first + gradient.second*gradient.second);
	if(norm > 0) {
		normal[0] /= norm;
		normal[1] /= norm;
	}
	return { normal[0], normal[1] };
}

  // Advect the level‑set field phi by normal speed V over time dt.
  // Incomplete need more testing
 
void advect_level_set(
    std::vector<std::vector<double>>& phi,
    double V, 
    double dt,
    double dx,
    double dy)
{
    int Nx = phi.size();  // Nx is the dimension of the first index (x-direction)
    int Ny = phi[0].size();

    // Make a copy for the update
    std::vector<std::vector<double>> phi_new = phi;

    for(int i = 0; i < Nx; ++i) {
        for(int j = 0; j < Ny; ++j) {
            
            auto grad = level_set_gradient(i, j,phi, dx, dy);
            double grad_mag = std::sqrt(grad.first*grad.first 
                                      + grad.second*grad.second);

            // Upwind Euler step: 
            phi_new[i][j] = phi[i][j] - dt * V * grad_mag;
        }
    }

    phi.swap(phi_new);

    // Optional: reinitialize to maintain signed‑distance property
    // reinitialize_level_set(phi, dx, dy);
}

// TODO: Need to implement a reinitialization function to maintain the signed distance property of the level set function.


// Functions above are needed for the level set method, and require further testing.
//________________________________________________________________________________________________________________________________________
//________________________________________________________________________________________________________________________________________
// Functions below are not specific to the level set method, and should work without error.

// Helper function for getting the voxel indices of a cell
std::pair<double,double> voxel_indices(Cell* pCell,
                             Phenotype& phenotype,
                             double dt)

{
	int v = pCell->get_current_voxel_index();
    auto& vox = microenvironment.mesh.voxels[v];
    double x = vox.center[0];
    double y = vox.center[1];

    // locate grid indices
    int i = (int)floor((x - ls_xmin) / ls_dx);
    int j = (int)floor((y - ls_ymin) / ls_dy);

    // clamp to valid range
    i = std::max(0, std::min(i, (int)level_set_phi.size()-1));
    j = std::max(0, std::min(j, (int)level_set_phi[0].size()-1)); //should solve out of bounds prob

	return {i, j};
}

// Distance should just be the value of SDF at the voxel???
double distance_to_membrane(Cell* pCell,
                             Phenotype& phenotype,
                             double dt)
{
    std::pair<int,int> indices = voxel_indices(pCell, phenotype, dt);
    int i = indices.first;
    int j = indices.second;

	return level_set_phi[i][j]; // SDF value at the voxel
}

void basement_membrane_interaction(Cell* pCell,
                                   Phenotype& phenotype,
                                   double dt)
{
    std::pair<int,int> indices = voxel_indices(pCell, phenotype, dt);
    int i = indices.first;
    int j = indices.second;
    double d = level_set_phi[i][j];
    double L = parameters.doubles("membrane_interaction_length");
    if(fabs(d) >= L) return;

    auto grad   = level_set_gradient( i, j,level_set_phi, ls_dx, ls_dy);
    auto normal = level_set_normalize(grad);
    double sign     = (d < 0.0 ? +1.0 : -1.0);
    double strength = parameters.doubles("membrane_adhesion_strength");
    double mag      = strength * (L - fabs(d));

    pCell->velocity[0] += sign * mag * normal.first;
    pCell->velocity[1] += sign * mag * normal.second;
}


void custom_rule( Cell* pCell, Phenotype& phenotype, double dt )
{	

	return; // This is a custom rule that can be used to implement any custom behavior for the cell, such a proliferation.
}

// Tested and works!
std::vector<std::vector<double>> generate_boundary_shape(double a, double b, double amp, int freq, int num_points){

	std::vector<std::vector<double>> pts;
	pts.reserve(num_points);
    
	// Generating points for a deformed ellipse like shape
    for (int i = 0; i < num_points; i++) {
        
        double theta = 2.0 * M_PI * i / (num_points - 1);
        
        double r_x = a * (1.0 + amp * cos(freq * theta));
        double r_y = b * (1.0 + amp * sin(freq * theta));
        double x = r_x * cos(theta);
        double y = r_y * sin(theta);

        pts.push_back({ x, y, 0.0 });
    }
    return pts;
}

// Helper to test if point (x,y) is inside the polygon 
bool is_inside(double x, double y, const std::vector<std::vector<double>>& BM_pts) {
	// std::cout << "We are inside the is_inside function!!!" << std::endl;
	// std::cout << "Checking if point (" << x << ", " << y << ") is inside the polygon." << std::endl;
	
	bool inside = false;
	int Np = (int)BM_pts.size();

	for(int k=0; k<Np; ++k) {                                         // Looping through each edge of the polygon

		double x1 = BM_pts[k][0], y1 = BM_pts[k][1];           
		double x2 = BM_pts[(k+1)%Np][0], y2 = BM_pts[(k+1)%Np][1];    // (x1,y1) and (x2,y2) are the endpoints of the edge

		// Check if horizontal ray intersects the edge between (x1,y1) and (x2,y2)
		// 
		if(((y1 > y) != (y2 > y)) ) {                                  //This should handle edges cases now

			double xint = x1 + (y - y1)*(x2 - x1)/(y2 - y1);           // Note here xint is the x-coordinate of the intersection point of the horizontal ray with the edge
			
			if(x < xint) inside = !inside;                             // Reverse logic: ray extends to the right
		}
	}
return inside;
};

void initialize_level_set_duct(){
	std::cout << "Initializing level set function for the duct..." << std::endl;

	auto& mesh = get_default_microenvironment()->mesh;
	ls_xmin = mesh.bounding_box[0]; // bounding box is [xmin ymin zmin xmax ymax zmax]
	ls_ymin = mesh.bounding_box[1];
	ls_dx = mesh.dx;
	ls_dy = mesh.dy;     // Uniform grid spacing in x and y directions

	// Compute number of voxels
	int Nx = (int)((mesh.bounding_box[3] - ls_xmin) / ls_dx);
	int Ny = (int)((mesh.bounding_box[4] - ls_ymin) / ls_dy); // Might need +1 here to include bb

	// Create SDF
	level_set_phi.assign(Nx, std::vector<double>(Ny, 0.0)); // Nx x Ny matrix initialized to zero

	// Get initial boundary points 
	// Create new shape for the basement membrane
	// Should probably put in header later
	
	int num_points = parameters.ints("membrane_num_points");
	double a = 300.0, b = 250.0;   // Length of ellipse axes
	double amp = 0.1;              // Amplitude of deformation
	int freq = 4;                 // Number of bumps

	boundary_membrane_pts = generate_boundary_shape(a, b, amp, freq, num_points);

	// Precompute segments for convenience
	// Segments are pairs of points representing the edges of the shape

	int Np = (int)boundary_membrane_pts.size();

	std::vector<std::pair< std::array<double,3>, std::array<double,3> >> segments;     // Segments are edges of the shape

	for(int k=0; k<Np; ++k) {
		auto p1_vec = boundary_membrane_pts[k];
		auto p2_vec = boundary_membrane_pts[(k+1)%Np];
		std::array<double,3> p1 = { p1_vec[0], p1_vec[1], p1_vec[2] };
		std::array<double,3> p2 = { p2_vec[0], p2_vec[1], p2_vec[2] }; // Check if should be array or vectors. For some reason need to be arrays not vectors (not sure why)
		segments.push_back({p1, p2});
}

	// Compute signed distance for each grid point
	// High level overview: This double for loop gets the x and y coordinates of each voxel
	// Then computes the minimum distance to any segment (edge) of the shape

	for(int i=0; i<Nx; i++){
		for(int j=0; j<Ny; j++){

			// Map voxel index to physical coordinate (half for cell center)
			double x = ls_xmin + (i + 0.5)*ls_dx;
			double y = ls_ymin + (j + 0.5)*ls_dy;

			// Compute min distance to any segment of the shape
			double minDist = std::numeric_limits<double>::infinity();
			for(auto &seg : segments){
				double x1 = seg.first[0], y1 = seg.first[1];
				double x2 = seg.second[0], y2 = seg.second[1];

				// Project point onto line segment (similar logic to is_inside)
				double vx = x2 - x1, vy = y2 - y1;              // change from start and end of edge
				double wx = x - x1, wy = y - y1;                // change from start of edge to voxel point (x,y)
				double v2 = vx*vx + vy*vy;
				double t = (v2>0 ? (wx*vx + wy*vy)/v2 : 0.0);   // t is the projection factor along the edge

				double dist;

				if(t <= 0) {
					dist = sqrt(wx*wx + wy*wy);
				} else if(t >= 1) {
					double ux = x - x2, uy = y - y2;
					dist = sqrt(ux*ux + uy*uy);
				} else {
					double px = x1 + t*vx, py = y1 + t*vy;
					double ux = x - px, uy = y - py;
					dist = sqrt(ux*ux + uy*uy);
				}
				minDist = std::min(minDist, dist);
			}
			bool inside = is_inside(x,y, boundary_membrane_pts);

			// Signed distance: negative inside, positive outside
			level_set_phi[i][j] = (inside ? -minDist : +minDist);
		}
	}
	std::cout << "Level set function initialized!!!" << std::endl;
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
		
	cell_defaults.functions.add_cell_basement_membrane_interactions = basement_membrane_interaction; 
	cell_defaults.functions.calculate_distance_to_membrane = distance_to_membrane; 

	build_cell_definitions_maps(); 

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
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();

	// Create new shape for the basement membrane
	// int num_points = parameters.ints("membrane_num_points");
	// double a = 350.0, b = 250.0;
	// int n_pts = 200;
	// double amp = 0.1;   // put in header later
	// int lobes = 4;

	// auto BM_pts = generate_boundary_shape(a, b, amp, lobes, num_points);
	// Cell_Definition* pBM_def = cell_definitions_by_index[1];

	// for (auto &pt : BM_pts){
	// 	Cell* pBM = create_cell( *pBM_def );
	// 	pBM->assign_position( pt );
	// }

	for (auto pCell : *all_cells){

		// Stash the volume growth rates in custom data
		pCell->custom_data["default_cyto_rate"] = pCell->phenotype.volume.cytoplasmic_biomass_change_rate;
		pCell->custom_data["default_nuclear_rate"] = pCell->phenotype.volume.nuclear_biomass_change_rate;
		pCell->custom_data["default_fluid_rate"] = pCell->phenotype.volume.fluid_change_rate;

	}
	
	// Initialize the level set representation of the duct
    // double duct_radius = parameters.doubles("duct_radius"); // Example: get from XML
    // std::vector<double> duct_center = {0.0, 0.0, 0.0};
    
    initialize_level_set_duct();

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