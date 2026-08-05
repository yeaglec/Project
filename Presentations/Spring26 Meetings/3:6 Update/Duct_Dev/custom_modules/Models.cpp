
#include "./Models.h" 
#include "./Utils.h"  
#include "./custom.h"
#include <cmath>
#include <cfloat>
#include <iostream>
//___________________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________________
// Models for the Level Set Method
// ___________________________________________________________________________________________________________________________

// ###################### Function that implements Cell-to-BM force ####################
void cell_interactions_LSM(Cell* pCell,Phenotype& phenotype, double dt){

    std::pair<int,int> indices = voxel_indices(pCell);
    int i = indices.first;
    int j = indices.second;
    double d = level_set_phi[i][j];
	

    double L = parameters.doubles("membrane_interaction_length"); // 500 now 
	double R = pCell->phenotype.geometry.radius;
	double de = d - (d < 0 ? -R : R);

    if(fabs(de) >= L) return;

	double cell_deadzone = parameters.doubles("cell_deadzone");
	
	if (fabs(de) <cell_deadzone) return;

	// Make this spring force

    auto grad   = level_set_gradient( i, j,level_set_phi, ls_dx, ls_dy);
    auto normal = level_set_normalize(grad);
    double sign     = (d < 0.0 ? +1.0 : -1.0);
    double strength = parameters.doubles("membrane_adhesion_strength"); //.001 right now
    double mag      = strength * fabs(de);

    pCell->velocity[0] += sign * mag * normal.first;
    pCell->velocity[1] += sign * mag * normal.second;
}

// ________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// Level set functions for membrane interactions (voxel-based)
// ________________________________________________________________________________________________________________________

// Distance should just be the value of SDF at the voxel
double distance_to_membrane(Cell* pCell,Phenotype& phenotype,double dt){
    std::pair<int,int> indices = voxel_indices(pCell);
    int i = indices.first;
    int j = indices.second;

	return level_set_phi[i][j]; // SDF value at the voxel
}

void test_enforce_boundary_repulsion(Cell* pCell,Phenotype& phenotype,double dt){

	double cellx = pCell->position[0];
	double celly = pCell->position[1];

}

// ################ Function that computes the BM-to-Cell force ####################

std::pair<double,double> basement_membrane_interactions_LSM(Cell* pCell){       // TODO: Implement into update_basement_membrane_interactions{
	std::pair<int,int> indices = voxel_indices(pCell);
    int i = indices.first;
    int j = indices.second;

    double d = level_set_phi[i][j];
	// std::cout << "Distance to membrane: " << d << std::endl;
	double R = pCell->phenotype.geometry.radius;
	// std::cout << "Cell radius: " << R << std:: endl;
	double de = d - (d < 0 ? -R : R);                                    // de is the offseted distance to the BM
	// std::cout << "Effective distance to membrane: " << de << std::endl;

    double L = parameters.doubles("membrane_interaction_length");        // 500 now 
	double strength = parameters.doubles("membrane_spring_constant");  //.001 right now

	double BM_deadzone = parameters.doubles("membrane_deadzone");

	if (fabs(de) <BM_deadzone) return {0.0, 0.0}; 

    auto grad   = level_set_gradient( i, j,level_set_phi, ls_dx, ls_dy); // TThis is the gradient from the voxel center not the cell center
    auto normal = level_set_normalize(grad);
    double sign = (d < 0.0 ? +1.0 : -1.0);

    double mag = strength * fabs(de);    // ThisS is Hooke's law, F = kx

    double Fx = sign * mag * normal.first;
    double Fy = sign * mag * normal.second;

	return { Fx, Fy };
}

// ________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// Functions for implementing deformations of basement membrane
// ________________________________________________________________________________________________________________________


void rebuild_signed_distance_field()
{
	// std::cout << "Rebuilding signed distance field..." << std::endl;
    // Mesh dimensions and spacing (already set in initialize duct code, just copied)
	
    int Nx = (int) level_set_phi.size();
    int Ny = (int) level_set_phi[0].size();

    // Build segment list from the current boundary points
    int Np = (int) boundary_membrane_pts.size();

    // For each grid cell, compute min distance to any segment
    for(int i = 0; i < Nx; ++i){
		for(int j = 0; j < Ny; ++j){
            double x = ls_xmin + (i + 0.5)*ls_dx;
			double y = ls_ymin + (j + 0.5)*ls_dy;
			auto [minDist, px, py, k, t] = project_point_onto_boundary(x, y);

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

		std::pair <double,double> force = basement_membrane_interactions_LSM(pCell);
		double Fx_cell = force.first;
		double Fy_cell = force.second;

		double Fx_BM = -Fx_cell;
		double Fy_BM = -Fy_cell;

		auto [best_dist, best_px, best_py, best_k_d, best_t] = project_point_onto_boundary(cell_x, cell_y);
		int best_k = static_cast<int>(std::round(best_k_d));

		// What this is doing:  take each cell’s force on the membrane, 
		//find which segment it hits, and split that tug between the two end‑nodes of that segment.

		int n1 = best_k, n2 = (best_k+1) % Np;
		double t_clamped = std::max(0.0, std::min(best_t, 1.0));  // Clamp t to [0,1]
		node_forces[n1].first += (1.0 - t_clamped) * Fx_BM;
		node_forces[n1].second += (1.0 - t_clamped) * Fy_BM;   // Split forces linearly based on t
		node_forces[n2].first += (t_clamped) * Fx_BM;
		node_forces[n2].second += (t_clamped) * Fy_BM;

	}

	for(int i=0; i<Np; i++) {
    boundary_membrane_pts[i][0] += node_forces[i].first  * dt;
    boundary_membrane_pts[i][1] += node_forces[i].second * dt;
	}

	// Rebuild the signed distance field after updating boundary points
	rebuild_signed_distance_field();

}

void initialize_level_set_duct(std::vector<std::vector<double>> boundary_membrane_pts){
	std::cout << "Initializing level set function for the duct..." << std::endl;

	auto& mesh = get_default_microenvironment()->mesh;
	ls_xmin = mesh.bounding_box[0]; // bounding box is [xmin ymin zmin xmax ymax zmax]
	ls_ymin = mesh.bounding_box[1];
	ls_dx = mesh.dx;
	ls_dy = mesh.dy;     // Uniform grid spacing in x and y directions
	
	int Nx = (int)((mesh.bounding_box[3] - ls_xmin) / ls_dx);
	int Ny = (int)((mesh.bounding_box[4] - ls_ymin) / ls_dy); // Might need +1 here to include bb

	level_set_phi.assign(Nx, std::vector<double>(Ny, 0.0)); // Nx x Ny matrix initialized to zer0

	int Np = (int)boundary_membrane_pts.size();
	initial_edge_length.resize(Np);
	initial_node_positions.resize(Np);

	// Store initial edge lengths for elasticity calculations
	for(int k = 0; k < Np; ++k)
	{
		initial_node_positions[k] = boundary_membrane_pts[k]; // Store initial node positions
		int j = (k + 1) % Np;  
		double dx = boundary_membrane_pts[j][0] - boundary_membrane_pts[k][0];
		double dy = boundary_membrane_pts[j][1] - boundary_membrane_pts[k][1];
		initial_edge_length[k] = std::sqrt(dx*dx + dy*dy);
	}

	// High level overview: This double for loop gets the x and y coordinates of each voxel
	// Then computes the minimum distance to any segment (edge) of the shape

	for(int i=0; i<Nx; i++){
		for(int j=0; j<Ny; j++){

			// Coordinates at voxel center
			double x = ls_xmin + (i + 0.5)*ls_dx;
			double y = ls_ymin + (j + 0.5)*ls_dy;

			auto [minDist, px, py, k, t] = project_point_onto_boundary(x,y); // Project point onto boundary to get min distance
			bool inside = is_inside(x,y, boundary_membrane_pts);

			// Signed distance: negative inside, positive outside
			level_set_phi[i][j] = (inside ? -minDist : +minDist);
		}
	}
	std::cout << "Level set function initialized!!!" << std::endl;
}
// Build Gaussian weights centered at projection point and distribute force 
void BM_Smoothing(std::vector<std::pair<double,double>>& node_forces, double Fx_BM, double Fy_BM, int best_k, double best_px, double best_py, double best_t){

		int Np = (int)boundary_membrane_pts.size();
        double sumW = 0.0;
		double sigma = parameters.doubles("membrane_force_smoothing_sigma"); // Smoothing parameter for Gaussian distribution
		double cutoff = 3.0 * sigma;         // ignore nodes beyond 3 sigma 
    
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

// Segment Elasticity: Add spring forces between adjacent nodes to maintain membrane integrity
void membrane_elasticity(std::vector<std::pair<double,double>>& node_forces)
{
    if (parameters.doubles("Segment_Elasticity") != 1.0) return;

    int Np = (int)boundary_membrane_pts.size();
    double k = parameters.doubles("seg_lin");
    double alpha = parameters.doubles("seg_exp");
    double max_force = 100.0; 

    for(int i=0; i<Np; ++i){
        int j = (i + 1) % Np;

        double dx = boundary_membrane_pts[j][0] - boundary_membrane_pts[i][0];
        double dy = boundary_membrane_pts[j][1] - boundary_membrane_pts[i][1];

        double current_length = sqrt(dx*dx + dy*dy);
        double rest_length = initial_edge_length[i];
        if (current_length <= 1e-12) continue;

        double x = current_length - rest_length;
        double F_mag = k * (std::exp(alpha * x) - 1.0);
        
        if(F_mag > max_force) F_mag = max_force;
        if(F_mag < -max_force) F_mag = -max_force;

        double nx = dx / current_length;
        double ny = dy / current_length;

        node_forces[i].first  += F_mag * nx;
        node_forces[i].second += F_mag * ny;
        node_forces[j].first  -= F_mag * nx;
        node_forces[j].second -= F_mag * ny;
    }
}

// Restoring Force: Force enforcing global structure of membrane
void membrane_restoring_force(std::vector<std::pair<double,double>>& node_forces)
{
    double home_lin = parameters.doubles("home_lin");
	double home_exp = parameters.doubles("home_exp");
    if (home_lin <= 0.0) return;

    int Np = (int)boundary_membrane_pts.size();
    for(int i=0; i<Np; ++i){
        double home_x = initial_node_positions[i][0];
        double home_y = initial_node_positions[i][1];
        double current_x = boundary_membrane_pts[i][0];
        double current_y = boundary_membrane_pts[i][1];

        double dx = home_x - current_x;
        double dy = home_y - current_y;

		double disp = sqrt(dx*dx + dy*dy);
		if (disp < 1e-12) continue;

		double F_mag = home_lin * disp;

		// double F_mag = home_lin * (std::exp(home_exp * disp) - 1.0);

		// std::cout << "Node Forces: " << node_forces[2].second << std::endl;
		// std::cout << "Resotring Force Mag: " << F_mag << std::endl; 

        node_forces[i].first  += F_mag * (dx / disp);
        node_forces[i].second += F_mag * (dy / disp);
    }
}
