
#include "./Membrane.h" 
#include "./Geometry.h"
#include "./Utils.h"
#include "./Tests.h"
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <iostream>


// ________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// TODO: Decouple LSM function from core code
// These 3 functions are still used in custom.cpp, need to be decoupled and archived 
// ________________________________________________________________________________________________________________________

// ########## Helper function for getting the voxel indices of a cell ###########
std::pair<double,double> voxel_indices(Cell* pCell){
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

// Distance should just be the value of SDF at the voxel
double distance_to_membrane(Cell* pCell,Phenotype& phenotype,double dt){
    std::pair<int,int> indices = voxel_indices(pCell);
    int i = indices.first;
    int j = indices.second;

	return level_set_phi[i][j]; // SDF value at the voxel
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



// ################ Function that computes the BM-to-Cell force ####################
std::pair<double,double> basement_membrane_interactions_cc(Cell* pCell, double de, double px, double py)
{
    double cell_x = pCell->position[0];
    double cell_y = pCell->position[1];

	// if (de > 0) {

	// 	return {0.0, 0.0}; // Remove deformation if cells pass through boundary
    // }

    double BM_deadzone = parameters.doubles("membrane_deadzone");
    if (fabs(de) < BM_deadzone) return {0.0, 0.0};

    double L = parameters.doubles("membrane_interaction_length");
    if (fabs(de) >= L) return {0,0};

    double strength = parameters.doubles("membrane_spring_constant"); 
    double mag = strength * fabs(de);  // Hooke's law, F = kx

    double nx = px - cell_x; // Normal vector from cell to boundary
    double ny = py - cell_y;
    double norm = sqrt(nx * nx + ny * ny);
	
    if (norm > 0) { nx /= norm; ny /= norm; }

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

    auto cell_force = basement_membrane_interactions_cc(pCell, de, px, py);
    pCell->custom_data[BM_Fx_idx] = cell_force.first;
    pCell->custom_data[BM_Fy_idx] = cell_force.second;
    pCell->custom_data[BM_k_idx] = k;
    pCell->custom_data[BM_px_idx] = px;
    pCell->custom_data[BM_py_idx] = py;
    pCell->custom_data[BM_t_idx] = t;

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

        if(parameters.doubles("is_lumenal_pressure") == 1){
            // Calculate normalized outward vector for lumenal pressure
            double nx_out = px - cell_x;
            double ny_out = py - cell_y;
            double norm_out = sqrt(nx_out * nx_out + ny_out * ny_out);
            
            if (norm_out > 1e-16)
            {
                nx_out /= norm_out;
                ny_out /= norm_out;
            }

            // double pressure_mag = parameters.doubles("lumenal_pressure_strength");
            // // Apply constant lumenal pressure outward
            // pCell->velocity[0] += pressure_mag * nx_out;
            // pCell->velocity[1] += pressure_mag * ny_out;

            // Apply a progressive lumenal pressure
            // 2. Scale the pressure based on depth into the lumen
            // If 'de' is highly negative when deep in the lumen, we use fabs(de)
            // We can optionally divide by a 'max_lumen_radius' to normalize it between 0 and 1
            
            double pressure_base_mag = parameters.doubles("lumenal_pressure_strength");
            
            // Only apply outward pressure to cells that have detached/stratified (de < 0)
            if (de < 0) 
            {
                // The further into the lumen, the higher the outward push
                double depth_scale = fabs(de); 
                double applied_pressure = pressure_base_mag * depth_scale;

                pCell->velocity[0] += applied_pressure * nx_out;
                pCell->velocity[1] += applied_pressure * ny_out;
            }
        }

        if (de >0)
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

void update_basement_membrane_deformation(double dt)
{
    int Np = (int)boundary_membrane_pts.size();
    std::vector<std::pair<double,double>> node_forces(Np, {0.0, 0.0});

    // Iterate over all cells and compute membrane forces distributed by Gaussian
    for (Cell* pCell : *all_cells) {

        double Fx_cell = pCell->custom_data[BM_Fx_idx];
        double Fy_cell = pCell->custom_data[BM_Fy_idx];

        if (Fx_cell == 0.0 && Fy_cell == 0.0) continue;

        double Fx_BM = -Fx_cell;   
        double Fy_BM = -Fy_cell;

        // Consume the cached projection data
        double best_k_d = pCell->custom_data[BM_k_idx];
        double best_px = pCell->custom_data[BM_px_idx];
        double best_py = pCell->custom_data[BM_py_idx];
        double best_t = pCell->custom_data[BM_t_idx];
		int best_k = static_cast<int>(std::round(best_k_d));

        if (parameters.doubles("is_gaussian_smoothing") == 1) BM_Smoothing(node_forces, Fx_BM, Fy_BM, best_k, best_px, best_py, best_t);
		else{

        // Direct Local Force Transfer (No Smoothing) 
        int k1 = best_k;
        int k2 = (best_k + 1) % Np;

        // We distribute the force to the two segment endpoints based on this distance.
        node_forces[k1].first  += Fx_BM * (1.0 - best_t);
        node_forces[k1].second += Fy_BM * (1.0 - best_t);          // best_t (0.0 to 1.0) represents where on the segment the cell projects.


        node_forces[k2].first  += Fx_BM * best_t;
        node_forces[k2].second += Fy_BM * best_t;
        }
	}

	// Enforcing membrane elasticity between node pairs
    if(parameters.doubles("is_strain_lin")==1) membrane_strain_lin(node_forces);
    else membrane_strain_exp(node_forces);

    if(parameters.doubles("is_bending_stiffness")==1) membrane_bending_stiffness(node_forces);

	// Establishing membrane "memory" or "home" force, keeping the membrane from deforming inwards too much
    if(parameters.doubles("is_restore_lin")==1) membrane_restoring_force_lin(node_forces);
    else membrane_restoring_force_exp(node_forces);

    ///_________
    // Put Membrane Test Functions Here

    if (test_perb != nullptr){
        test_perb(node_forces, PhysiCell_globals.current_time);
    }

	// Update node positions
    for (int i = 0; i < Np; ++i) {
        boundary_membrane_pts[i][0] += node_forces[i].first  * dt;
        boundary_membrane_pts[i][1] += node_forces[i].second * dt;
    }

    if (parameters.doubles("is_add_membrane_nodes") == 1) add_membrane_nodes();

    // Rebuild signed distance field after modifications
    rebuild_signed_distance_field();
}

// Build Gaussian weights centered at projection point and distribute force 
// Possible need to think about revising for effeciency 
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
			// Start at best_k and go "left" (decreasing indices)
			for (int step = 0; step < Np; ++step) {
				int node_i = (best_k - step % Np + Np) % Np;
				double dx = boundary_membrane_pts[node_i][0] - best_px;
				double dy = boundary_membrane_pts[node_i][1] - best_py;
				double dist2 = dx * dx + dy * dy;
				if (dist2 > cutoff * cutoff) break; 
				
				double w = exp(-dist2 / (2.0 * sigma * sigma));
				nearby_indices.push_back(node_i);
				nearby_weights.push_back(w);
				sumW += w;
			}
			
			int num_left = nearby_indices.size();
			// Start at best_k + 1 and go "right" (increasing indices)
			for (int step = 1; step <= Np - num_left; ++step) {
				int node_i = (best_k + step) % Np;
				double dx = boundary_membrane_pts[node_i][0] - best_px;
				double dy = boundary_membrane_pts[node_i][1] - best_py;
				double dist2 = dx * dx + dy * dy;
				if (dist2 > cutoff * cutoff) break; 
				
				double w = exp(-dist2 / (2.0 * sigma * sigma));
				nearby_indices.push_back(node_i);
				nearby_weights.push_back(w);
				sumW += w;
			}

			if (sumW > 0.0) {
				// normalize weights and distribute
				for (size_t idx = 0; idx < nearby_indices.size(); ++idx) {
					int node_i = nearby_indices[idx];
					double w = nearby_weights[idx];
					double frac = w / sumW;
					node_forces[node_i].first  += Fx_BM * frac;
					node_forces[node_i].second += Fy_BM * frac;
				}
			} else {
				// Fallback to 2-node split if no nodes were within cutoff
				int n1 = best_k;
				int n2 = (best_k + 1) % Np;
				double t_clamped = std::max(0.0, std::min(best_t, 1.0));
				node_forces[n1].first  += (1.0 - t_clamped) * Fx_BM;
				node_forces[n1].second += (1.0 - t_clamped) * Fy_BM;
				node_forces[n2].first  += (t_clamped) * Fx_BM;
				node_forces[n2].second += (t_clamped) * Fy_BM;	
			}
        }
    } 


// Segment Elasticity: Add spring forces between adjacent nodes to maintain membrane integrity
void membrane_strain_lin(std::vector<std::pair<double,double>>& node_forces)    
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
        if (current_length <= 1e-16) continue;         //physicell standard is 1e-16

        double x = current_length - rest_length;
        double F_mag = k * x;        // parameterize linear and exponetial cases 
        
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

// Segment Elasticity: Add spring forces between adjacent nodes to maintain membrane integrity
void membrane_strain_exp(std::vector<std::pair<double,double>>& node_forces)    
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
        if (current_length <= 1e-16) continue;         //physicell standard is 1e-16

        double x = current_length - rest_length;
        double F_mag = k * (std::exp(alpha * x) - 1.0);        // parameterize linear and exponetial cases 
        
        if(F_mag > max_force) F_mag = max_force;
        if(F_mag < -max_force) F_mag = -max_force;    // Force cap for stability

        double nx = dx / current_length;
        double ny = dy / current_length;

        node_forces[i].first  += F_mag * nx;
        node_forces[i].second += F_mag * ny;
        node_forces[j].first  -= F_mag * nx;
        node_forces[j].second -= F_mag * ny;
    }
}

// Restoring Force: Force enforcing global structure of membrane
void membrane_restoring_force_lin(std::vector<std::pair<double,double>>& node_forces)
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
		if (disp < 1e-16) continue;

		double F_mag = home_lin * disp;                                        // parameterize linear and exponetial cases 

		// double F_mag = home_lin * (std::exp(home_exp * disp) - 1.0);

		// std::cout << "Node Forces: " << node_forces[2].second << std::endl;
		// std::cout << "Resotring Force Mag: " << F_mag << std::endl; 

        node_forces[i].first  += F_mag * (dx / disp);
        node_forces[i].second += F_mag * (dy / disp);
    }
}
void membrane_restoring_force_exp(std::vector<std::pair<double,double>>& node_forces)
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
		if (disp < 1e-16) continue;

		double F_mag = home_lin * (std::exp(home_exp * disp) - 1.0);

		// std::cout << "Node Forces: " << node_forces[2].second << std::endl;
		// std::cout << "Resotring Force Mag: " << F_mag << std::endl; 

        node_forces[i].first  += F_mag * (dx / disp);
        node_forces[i].second += F_mag * (dy / disp);
    }
}

void membrane_bending_stiffness(std::vector<std::pair<double,double>>& node_forces)
{
    int Np = (int)boundary_membrane_pts.size();
    double kb = parameters.doubles("membrane_bending_constant"); 
	double max_force = 100.0;

    for (int i = 0; i < Np; ++i) {
        int prev_node = (i - 1 + Np) % Np;
        int next_node = (i + 1) % Np;

        // Calculate the midpoint of the  neighbors
        double mid_x = 0.5 * (boundary_membrane_pts[prev_node][0] + boundary_membrane_pts[next_node][0]);
        double mid_y = 0.5 * (boundary_membrane_pts[prev_node][1] + boundary_membrane_pts[next_node][1]);

        // Calculate the vector from the current node to that midpoint (the deviation)
        double dev_x = mid_x - boundary_membrane_pts[i][0];
        double dev_y = mid_y - boundary_membrane_pts[i][1];

        // The bending force is proportional to this deviation
        double F_bx = kb * dev_x;
        double F_by = kb * dev_y;

		if(F_bx > max_force) F_bx = max_force;
        if(F_bx < -max_force) F_bx = -max_force;    // Force cap for stability
		if(F_by > max_force) F_by = max_force;
        if(F_by < -max_force) F_by = -max_force;    // Force cap for stability

        // Apply the restorative bending force to the current node
        node_forces[i].first  += F_bx;
        node_forces[i].second += F_by;

        // Apply equal and opposite reaction forces to the neighbors 
        node_forces[prev_node].first  -= 0.5 * F_bx;
        node_forces[prev_node].second -= 0.5 * F_by;
        
        node_forces[next_node].first  -= 0.5 * F_bx;
        node_forces[next_node].second -= 0.5 * F_by;
    }
}

void add_membrane_nodes()
{
    // Simple remesh mode (default): split any edge longer than max_edge_length.
    // Future toggle point for a cell-body-aware gate:
    // const bool use_cell_body_remesh = (parameters.doubles("use_cell_body_remesh") == 1.0);
    const bool use_cell_body_remesh = false;
    double max_edge_length = parameters.doubles("max_edge_length");

    if (use_cell_body_remesh && all_cells == nullptr) return;

    for (int i = 0; i < (int)boundary_membrane_pts.size(); ++i) {
        int next = (i + 1) % boundary_membrane_pts.size();
        
        double dx = boundary_membrane_pts[next][0] - boundary_membrane_pts[i][0];
        double dy = boundary_membrane_pts[next][1] - boundary_membrane_pts[i][1];
        double dist = sqrt(dx*dx + dy*dy);

        bool should_split = (dist > max_edge_length);

        if (use_cell_body_remesh) {
            // Future cell-body-aware remesh gate (intentionally kept commented for now):
            //
            // double remesh_stretch_factor = parameters.doubles("remesh_stretch_factor");
            // double remesh_cell_overlap_fraction = parameters.doubles("remesh_cell_overlap_fraction");
            // int remesh_min_cell_matches = std::max(1, (int)std::round(parameters.doubles("remesh_min_cell_matches")));
            //
            // double rest_length = initial_edge_length[i];
            // if (rest_length <= 1e-16) should_split = false;
            //
            // double stretch_ratio = dist / rest_length;
            // if (stretch_ratio < remesh_stretch_factor) should_split = false;
            //
            // int supporting_cells = 0;
            // double inv_dist = 1.0 / dist;
            // for (Cell* pCell : *all_cells) {
            //     auto projection = project_point_onto_boundary(pCell->position[0], pCell->position[1]);
            //     int best_k = static_cast<int>(std::round(std::get<3>(projection)));
            //     if (best_k != i) continue;
            //
            //     double t_clamped = std::max(0.0, std::min(std::get<4>(projection), 1.0));
            //     double cell_radius = pCell->phenotype.geometry.radius;
            //     double t_span = cell_radius * inv_dist;
            //     double left = std::max(0.0, t_clamped - t_span);
            //     double right = std::min(1.0, t_clamped + t_span);
            //     double overlap_fraction = std::max(0.0, right - left);
            //
            //     if (overlap_fraction >= remesh_cell_overlap_fraction) {
            //         ++supporting_cells;
            //         if (supporting_cells >= remesh_min_cell_matches) break;
            //     }
            // }
            //
            // if (supporting_cells < remesh_min_cell_matches) should_split = false;
        }

        if (!should_split) continue;

        // Split this edge in half.
        std::vector<double> mid_pt = {
            0.5 * (boundary_membrane_pts[i][0] + boundary_membrane_pts[next][0]),
            0.5 * (boundary_membrane_pts[i][1] + boundary_membrane_pts[next][1]),
            0.0
        };

        std::vector<double> mid_home = {
            0.5 * (initial_node_positions[i][0] + initial_node_positions[next][0]),
            0.5 * (initial_node_positions[i][1] + initial_node_positions[next][1]),
            0.0
        };

        double halved_L0 = initial_edge_length[i] / 2.0;
        initial_edge_length[i] = halved_L0;

        boundary_membrane_pts.insert(boundary_membrane_pts.begin() + i + 1, mid_pt);
        initial_node_positions.insert(initial_node_positions.begin() + i + 1, mid_home);
        initial_edge_length.insert(initial_edge_length.begin() + i + 1, halved_L0);

        i++;
    }
}
