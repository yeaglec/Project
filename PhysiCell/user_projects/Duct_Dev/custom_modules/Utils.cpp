#include "./Utils.h"
#include "./Geometry.h"
#include <cmath>
#include <cfloat>
#include <iostream>
#include <array>

// ________________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________________
// Helpers 
// ________________________________________________________________________________________________________________________

// ######### Helper function to save boundary points to a CSV file #########
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



// Cameron Code for division parallel to segments
void parallel_cell_division( Cell* parent, Cell* child ){

	double cell_x = parent->position[0];
	double cell_y = parent->position[1];

	auto [best_dist, best_px, best_py, best_k_d, best_t] = project_point_onto_boundary(cell_x, cell_y);
	int best_k = static_cast<int>(std::round(best_k_d));
	int Np = (int)boundary_membrane_pts.size();

	// Build tangent from best segment
	int k = best_k;
	int kp = (k + 1) % Np;

	double tx = boundary_membrane_pts[kp][0] - boundary_membrane_pts[k][0];
    double ty = boundary_membrane_pts[kp][1] - boundary_membrane_pts[k][1];
    double tz = 0.0;

    double len = sqrt(tx*tx + ty*ty);
    tx /= len; ty /= len; 

	std::vector<double> &orient = parent->state.orientation;
	double polarity = parent -> phenotype.geometry.polarity;

	double dot = (tx * orient[0] + ty * orient[1] + tz * orient[2]); 
	tx -= polarity * dot * orient[0];
	ty -= polarity * dot * orient[1];
	tz -= polarity * dot * orient[2];

	// normalize again
	double norm = sqrt(tx*tx + ty*ty + tz*tz);
	if (norm > 0.0) {
		tx /= norm; ty /= norm; tz /= norm;
	}

	double radius = parent->phenotype.geometry.radius;
	double dx = radius * tx;
	double dy = radius * ty;
	double dz = radius * tz;

	// The core code is inconveniently positioning the parent and child cell using rand_vec
	// We gotta reconstruct the random vector to update the positions using our tangent vector
	// Core code does: child_pos = parent_pos + rand_vec
	// parent_pos = parent_final - 0.5 * rand_orig
	// Ahh, we can just do some algreba and get back rand_vec

	std::array<double,3> parent_final = { parent->position[0], parent->position[1], parent->position[2] };
    std::array<double,3> child_final  = { child->position[0],  child->position[1],  child->position[2] };

	// Recover rand_vec
    std::array<double,3> rand_orig;
    rand_orig[0] = (child_final[0] - parent_final[0]) / 1.5;  // rand_vec = (child_final - parent_final) / 1.5
    rand_orig[1] = (child_final[1] - parent_final[1]) / 1.5;
    rand_orig[2] = (child_final[2] - parent_final[2]) / 1.5;

	// Recover parent_orig
    std::array<double,3> parent_orig;
    parent_orig[0] = parent_final[0] + 0.5 * rand_orig[0];
    parent_orig[1] = parent_final[1] + 0.5 * rand_orig[1];
    parent_orig[2] = parent_final[2] + 0.5 * rand_orig[2];

    // New positions: child_new = parent_orig + d_des ; parent_new = parent_orig - 0.5*d_des
    std::array<double,3> child_new = { parent_orig[0] + dx, parent_orig[1] + dy, parent_orig[2] + dz };
    std::array<double,3> parent_new = { parent_orig[0] - 0.5*dx, parent_orig[1] - 0.5*dy, parent_orig[2] - 0.5*dz };

    child->assign_position( child_new[0], child_new[1], child_new[2] );
    parent->assign_position( parent_new[0], parent_new[1], parent_new[2] );

	// std::cout << "Parallel cell division completed!!!!!!" << std::endl;
};
