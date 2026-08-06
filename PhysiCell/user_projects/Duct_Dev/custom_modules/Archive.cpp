#include "Utils.h"
#include "./custom.h"
#include <cmath>
#include <cfloat>

// ________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// Level Set Methods Code (Needs to be Reimplemented)
// ________________________________________________________________________________________________________________________


// ######### Helper function to compute the gradient based on voxel center  ##############
std::pair<double, double> level_set_gradient(int i, int j, 
    const std::vector<std::vector<double>>& phi,
	double dx, double dy){
        
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

// ########### Helper function to normalize the gradient ##############
std::pair<double,double> level_set_normalize(const std::pair<double, double>& gradient){
	std::vector<double> normal = { gradient.first, gradient.second };
	double norm = sqrt(gradient.first*gradient.first + gradient.second*gradient.second);
	if(norm > 0) {
		normal[0] /= norm;
		normal[1] /= norm;
	}
	return { normal[0], normal[1] };
}

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

#include "./Archive.h"
#include "./Geometry.h"
#include "./Utils.h"
#include <cmath>
#include <cfloat>
#include <iostream>

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


