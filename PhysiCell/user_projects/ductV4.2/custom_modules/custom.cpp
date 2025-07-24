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
std::vector<double> initial_edge_length; 

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

// Helper function to compute the gradient (TODO: Verify if this is correct: Good)

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

// Helper function for getting the voxel indices of a cell
std::pair<double,double> voxel_indices(Cell* pCell)

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

std::pair<double,double> compute_cell_force(Cell* pCell)         // TODO: Implement into update_basement_membrane_interactions
{
	std::pair<int,int> indices = voxel_indices(pCell);
    int i = indices.first;
    int j = indices.second;
    double d = level_set_phi[i][j];

    double L = parameters.doubles("membrane_interaction_length"); // 500 now 
	double strength = parameters.doubles("membrane_adhesion_strength"); //.001 right now
	double d0 = 5.0;

    if(fabs(d) <= d0){
		return {0.0, 0.0}; // No force if close to the 
	}

    auto grad   = level_set_gradient( i, j,level_set_phi, ls_dx, ls_dy);
    auto normal = level_set_normalize(grad);
    double sign     = (d < 0.0 ? +1.0 : -1.0);
    
    double mag      = strength * fabs(d-d0);    // ThisS is Hooke's law, F = kx

    double Fx = sign * mag * normal.first;
    double Fy = sign * mag * normal.second;

	return { Fx, Fy };
}

// Distance should just be the value of SDF at the voxel???
double distance_to_membrane(Cell* pCell,
                             Phenotype& phenotype,
                             double dt)
{
    std::pair<int,int> indices = voxel_indices(pCell);
    int i = indices.first;
    int j = indices.second;

	return level_set_phi[i][j]; // SDF value at the voxel
}

void basement_membrane_interaction(Cell* pCell,
                                   Phenotype& phenotype,
                                   double dt)
{
    std::pair<int,int> indices = voxel_indices(pCell);
    int i = indices.first;
    int j = indices.second;
    double d = level_set_phi[i][j];
    double L = parameters.doubles("membrane_interaction_length"); // 500 now 
    if(fabs(d) >= L) return;

	// Make this spring force

    auto grad   = level_set_gradient( i, j,level_set_phi, ls_dx, ls_dy);
    auto normal = level_set_normalize(grad);
    double sign     = (d < 0.0 ? +1.0 : -1.0);
    double strength = parameters.doubles("membrane_adhesion_strength"); //.001 right now
    double mag      = strength * fabs(d);

    pCell->velocity[0] += sign * mag * normal.first;
    pCell->velocity[1] += sign * mag * normal.second;
}

// ________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// Functions for implementing deformations of basement membrane
// ________________________________________________________________________________________________________________________

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

void rebuild_signed_distance_field()
{
	// std::cout << "Rebuilding signed distance field..." << std::endl;
    // Mesh dimensions and spacing (already set in initialize_level_set_duct)
    int Nx = (int) level_set_phi.size();
    int Ny = (int) level_set_phi[0].size();

    // Build segment list from the current boundary points
    int Np = (int) boundary_membrane_pts.size();
    std::vector<std::pair<std::array<double,3>,std::array<double,3>>> segments;
    segments.reserve(Np);
    for(int k = 0; k < Np; ++k)
    {
        auto &p1_vec = boundary_membrane_pts[k];
        auto &p2_vec = boundary_membrane_pts[(k+1) % Np];
        segments.push_back({
            { p1_vec[0], p1_vec[1], p1_vec[2] },
            { p2_vec[0], p2_vec[1], p2_vec[2] }
        });
    }

    // For each grid cell, compute min distance to any segment
    for(int i = 0; i < Nx; ++i)
    {
        double x = ls_xmin + (i + 0.5) * ls_dx;
        for(int j = 0; j < Ny; ++j)
        {
            double y = ls_ymin + (j + 0.5) * ls_dy;

            double minDist = DBL_MAX;
            for(auto &seg : segments)
            {
                double x1 = seg.first[0], y1 = seg.first[1];
                double x2 = seg.second[0], y2 = seg.second[1];

                // projection factor t along segment
                double vx = x2 - x1, vy = y2 - y1;
                double wx = x  - x1, wy = y  - y1;
                double v2 = vx*vx + vy*vy;
                double t  = (v2 > 0 ? (wx*vx + wy*vy) / v2 : 0.0);

                // clamp t to [0,1] and compute distance
                double dist;
                if (t <= 0.0)
                {
                    dist = sqrt(wx*wx + wy*wy);               // to p1
                }
                else if (t >= 1.0)
                {
                    double ux = x - x2, uy = y - y2;
                    dist = sqrt(ux*ux + uy*uy);               // to p2
                }
                else
                {
                    double px = x1 + t * vx;
                    double py = y1 + t * vy;
                    double ux = x - px, uy = y - py;
                    dist = sqrt(ux*ux + uy*uy);               // to projection
                }
                minDist = std::min(minDist, dist);
            }

            //Determine sign via is_inside() and write phi
            bool inside = is_inside(x, y, boundary_membrane_pts);
            level_set_phi[i][j] = inside ? -minDist : +minDist;
        }
    }
}

void update_basement_membrane_deformation(double dt){

	// std::cout << "Updating basement membrane deformation..." << std::endl;

	int Np = (int)boundary_membrane_pts.size();
	std::vector<std::pair<double,double>> node_forces(Np,{0,0});

	for (Cell* pCell : *all_cells){

		double cell_x = pCell->position[0];
		double cell_y = pCell->position[1];

		std::pair <double,double> force = compute_cell_force(pCell);
		double Fx_cell = force.first;
		double Fy_cell = force.second;

		double Fx_BM = -Fx_cell;
		double Fy_BM = -Fy_cell;

		double best_dist = DBL_MAX;
		int best_k = 0;
		double best_t = 0;

		for (int k=0; k<Np; k++){
			
			// Endpoints for the segment from k to k+1
			int j = (k+1) % Np; 
			double x1 = boundary_membrane_pts[k][0], y1 = boundary_membrane_pts[k][1];
			double x2 = boundary_membrane_pts[j][0], y2 = boundary_membrane_pts[j][1];

			// Computing projection of vector c onto vector b

			double bx = x2 - x1, by = y2 - y1;                  // Vector b for segment
			double cx = cell_x-x1, cy = cell_y-y1;              // Vector c from first endpoint to cell position
			double b2 = bx*bx + by*by;
			double t = (b2>0 ? (cx*bx + cy*by) / b2 : 0.0);     // Just the projection formula: c x b / |b|^2

			double dist;
			if (t <= 0.0) {
				dist = sqrt(cx*cx + cy*cy); // Projection is before first point
			} 

			else if (t >= 1.0) {
				dist = sqrt((cell_x-x2)*(cell_x-x2) + (cell_y-y2)*(cell_y-y2)); // After second point
			} 
			
			else {
				double px = x1 + t * bx, py = y1 + t * by; 
				dist = sqrt((cell_x-px)*(cell_x-px) + (cell_y-py)*(cell_y-py)); // Between endpoints
			}

			if (dist < best_dist) {
				best_dist = dist;
				best_k = k;
				best_t = t;
			}
		}

		// What this is doing:  take each cell’s force on the membrane, 
		//find which segment it hits, and split that tug between the two end‑nodes of that segment.
		int n1 = best_k, n2 = (best_k+1) % Np;
		double t_clamped = std::max(0.0, std::min(best_t, 1.0));  // Clamp t to [0,1]
		node_forces[n1].first += (1.0 - t_clamped) * Fx_BM;
		node_forces[n1].second += (1.0 - t_clamped) * Fy_BM;   // Split forces linearly based on t
		node_forces[n2].first += (      t_clamped) * Fx_BM;
		node_forces[n2].second += (      t_clamped) * Fy_BM;

	}

	// // Adding elastic forces to the boundary membrane points, pulling them back towards the original shape
	// double k_spring = parameters.doubles("membrane_spring_constant");
	// for(int i=0; i<Np; i++) {
	// 	int j = (i+1) % Np;

	// 	// Computing the edge vector 
	// 	double dx = boundary_membrane_pts[j][0] - boundary_membrane_pts[i][0];
	// 	double dy = boundary_membrane_pts[j][1] - boundary_membrane_pts[i][1];
	// 	double length = sqrt(dx*dx + dy*dy);
	// 	double rest   = initial_edge_length[i];  // edge's resting length, undeformed

	// 	if (length == 0) continue;  // avoid divide-by-zero

	// 	double fspring = -k_spring * (length - rest);  // Hooke's law: F = -k * (x - x0)

	// 	double Fx = fspring * (dx/length);
	// 	double Fy = fspring * (dy/length);

	// 	node_forces[i].first +=  Fx;  node_forces[i].second +=  Fy;
	// 	node_forces[j].first += -Fx;  node_forces[j].second += -Fy;
	// }

	for(int i=0; i<Np; i++) {
    boundary_membrane_pts[i][0] += node_forces[i].first  * dt;
    boundary_membrane_pts[i][1] += node_forces[i].second * dt;
	}

	// Rebuild the signed distance field after updating boundary points
	rebuild_signed_distance_field();

}


// ________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// These are functions for intializing the level set function phi for the duct
// ________________________________________________________________________________________________________________________


// Tested and works!
std::vector<std::vector<double>> generate_boundary_shape(double a, double b, double amp, int freq, int num_points){

	std::vector<std::vector<double>> pts;
	pts.reserve(num_points);
    
	// Generating points for a deformed ellipse like shape
    for (int i = 0; i < num_points; i++) {
        
        double theta = 2.0 * M_PI * i / (num_points);
        
        double r_x = a * (1.0 + amp * cos(freq * theta));
        double r_y = b * (1.0 + amp * sin(freq * theta));
        double x = r_x * cos(theta);
        double y = r_y * sin(theta);

        pts.push_back({ x, y, 0.0 });
    }
    return pts;
}


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
	initial_edge_length.resize(Np);
	for(int k = 0; k < Np; ++k)
	{
		int j = (k + 1) % Np;  
		double dx = boundary_membrane_pts[j][0] - boundary_membrane_pts[k][0];
		double dy = boundary_membrane_pts[j][1] - boundary_membrane_pts[k][1];
		initial_edge_length[k] = std::sqrt(dx*dx + dy*dy);
	}

	std::vector<std::pair< std::array<double,3>, std::array<double,3> >> segments;     // Segments are edges of the shape

	for(int k=0; k<Np; ++k) {
		auto p1_vec = boundary_membrane_pts[k];
		auto p2_vec = boundary_membrane_pts[(k+1)%Np];
		std::array<double,3> p1 = { p1_vec[0], p1_vec[1], p1_vec[2] };
		std::array<double,3> p2 = { p2_vec[0], p2_vec[1], p2_vec[2] }; // Check if should be array or vectors. For some reason need to be arrays not vectors (not sure why)
		segments.push_back({p1, p2});
}
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
				double t = (v2>0 ? (wx*vx + wy*vy)/v2 : 0.0);   // this is just the projection formula: c x b / |b|^2

				double dist;

				if(t <= 0) {
					dist = sqrt(wx*wx + wy*wy);
				} else if(t >= 1) {
					double ux = x - x2, uy = y - y2;
					dist = sqrt(ux*ux + uy*uy);
				} else {
					double px = x1 + t*vx, py = y1 + t*vy;      // TODO: At some point make this projection code a function
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
	
	// Cell* pC;
	
	// for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	// {
	// 	Cell_Definition* pCD = cell_definitions_by_index[k]; 
	// 	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
	// 	for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
	// 	{
	// 		std::vector<double> position = {0,0,0}; 
	// 		position[0] = Xmin + UniformRandom()*Xrange; 
	// 		position[1] = Ymin + UniformRandom()*Yrange; 
	// 		position[2] = Zmin + UniformRandom()*Zrange; 
			
	// 		pC = create_cell( *pCD ); 
	// 		pC->assign_position( position );
	// 	}
	// }
	// std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml();
	set_parameters_from_distributions();

	//_______________________________________________________________________________________________________________________
	// Placing cells to test the basement membrane deformation

	// // Create new shape for the basement membrane
	// int num_ep_points = parameters.ints("membrane_num_points")/2;
	// double a = 300.0, b = 250.0;
	// int n_pts = 200;
	// double amp = 0.1;   // put in header later
	// int lobes = 4;

	// auto BM_pts = generate_boundary_shape(a, b, amp, lobes, num_ep_points);
	// Cell_Definition* pBM_def = cell_definitions_by_index[0];

	// int n_BM_cells = 0;
	// for (auto &pt : BM_pts){
	// 	if(n_BM_cells < 3){
	// 	Cell* pBM = create_cell( *pBM_def );
	// 	pBM->assign_position( pt );
	// 	}
	// 	n_BM_cells++;
	// }

	double a = 300.0, b = 250.0;
	std::vector<double> cm = { 0, 0, 0 };   // center of ellipse
	Cell_Definition* pBM_def = cell_definitions_by_index[0];
	for( int k=1; k<parameters.ints("number_EP_cells")+1; k++ )
	{
		double theta = 2.0 * M_PI * k / parameters.ints("number_EP_cells");
		// place daughters just *inside* the SDF: r = (a - gamma)
		double gamma = 30.0;                           
		double x = (a - gamma) * cos(theta);
		double y = (b - gamma) * sin(theta);
		Cell* pC = create_cell( *pBM_def );

		if( parameters.ints("number_EP_cells") > 1 ){
			pC->assign_position( { x, y, 0.0 } );
		}
		else{
			pC->assign_position( { parameters.doubles("x"), parameters.doubles("y"), 0.0 } );
		}
	}

	//_______________________________________________________________________________________________________________________

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