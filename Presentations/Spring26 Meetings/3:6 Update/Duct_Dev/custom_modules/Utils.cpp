#include "Utils.h"
#include "./custom.h"
#include <cmath>
#include <cfloat>

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

// ####### Helper to test if point (x,y) is inside the polygon ##########
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
			
			if(x < xint) inside = !inside;                             // Reversed logic: ray extends to the right!
		}
	}
return inside;
}

// ######### Function for finding to projection of a point onto the boundary #########
std::tuple<double, double, double, double, double> project_point_onto_boundary(double x, double y){   // TODO: Replace redundant projection code

	double best_dist = DBL_MAX;
	double best_px = 0.0, best_py = 0.0;
	int best_k = 0;
    double best_t = 0.0;
	int Np = boundary_membrane_pts.size();

	for (int k = 0; k < Np; ++k){

		int j = (k + 1) % Np;
		double x1 = boundary_membrane_pts[k][0], y1 = boundary_membrane_pts[k][1];
		double x2 = boundary_membrane_pts[j][0], y2 = boundary_membrane_pts[j][1];

		// Computing projection of vector c onto vector b

		double bx = x2 - x1, by = y2 - y1;                     // Vector b for segment
		double cx = x - x1, cy = y - y1;             // Vector c from first endpoint to cell position
		double b2 = bx*bx + by*by;
		double t = (b2 > 0.0 ? (cx*bx + cy*by) / b2 : 0.0);    // Just the projection formula: c \cdot b / ||b||^2

		double dist;
		double px, py;   // Our projection

		if (t <= 0.0) {
			// before first endpoint
			dist = sqrt(cx*cx + cy*cy);
			px = x1; py = y1;
		} else if (t >= 1.0) {
			// after second endpoint
			double ux = x - x2, uy = y - y2;
			dist = sqrt(ux*ux + uy*uy);
			px = x2; py = y2;
		} else {
			// projection between endpoints
			px = x1 + t * bx;
			py = y1 + t * by;
			double ux = x - px, uy = y - py;
			dist = sqrt(ux*ux + uy*uy);
		}

		if (dist < best_dist) {
			best_dist = dist;
			best_k = k;
			best_t = t;
			best_px = px;
			best_py = py;
		}
	}

	return {best_dist, best_px, best_py, best_k, best_t};

}

// Rudimentary function to generate duct shape
std::vector<std::vector<double>> generate_boundary_shape(double a, double b, double amp, int freq){

	int num_points = parameters.ints("membrane_num_points");
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

// Generating points for a deformed ellipse like shape
void generate_boundary_cells(double a, double b, double amp, int freq, std::string type, double dis, int num_cells){
	
	Cell_Definition* pBM_def = cell_definitions_by_name[type];   //cell_definitions_by_name[ type_name ] = pCD; 
	;

    for (int i = 0; i < num_cells; i++) {
        
        double theta = 2.0 * M_PI * i / (num_cells);
        
        double r_x = a * (1.0 + amp * cos(freq * theta));
        double r_y = b * (1.0 + amp * sin(freq * theta));
        double x = r_x * cos(theta);
        double y = r_y * sin(theta);

		// compute radial distance and unit‐radial direction
		double r_norm = std::sqrt(x*x + y*y);
		double nx = x / r_norm;    // outward radial unit vector
		double ny = y / r_norm;

		// step back along the normal by ep_dis
		double xi = x - dis * nx;
		double yi = y - dis * ny;

		Cell* pC = create_cell( *pBM_def );
		if( parameters.ints("number_EP_cells") > 1 ){
			pC->assign_position( { xi, yi, 0.0 } );
		}
		else{
			pC->assign_position( { parameters.doubles("x"), parameters.doubles("y"), 0.0 } );
		}
		
		if (i==0){

			if (type =="Epithelial"){
				std::cout << "Setting first cell to proliferate" << std::endl;
				pC->phenotype.cycle.data.exit_rate(0) = parameters.doubles("proliferation_exit_rate");
			}
			
			else{ 
				pC->phenotype.cycle.data.exit_rate(0) = 0;
			}
		}
		
		
	}
}

std::vector<std::vector<double>> generate_circle_boundary(){
	double radius = parameters.doubles("membrane_circle_radius");
	double num_points = parameters.doubles("membrane_num_points");
	std::vector<std::vector<double>> pts;	
	pts.reserve(num_points);

	for (int i = 0; i < num_points; i++) {
		double theta = 2.0 * M_PI * i / num_points;
		double x = radius * cos(theta);
		double y = radius * sin(theta);
		pts.push_back({ x, y, 0.0 });
	}
	return pts;
}

void generate_circle_cells(){
	
		int num_ep = parameters.ints("number_EP_cells");
		Cell_Definition* pBM_def = cell_definitions_by_index[0];
		double ep_dis = parameters.doubles("ep_displacement");

		for (int i=0; i<num_ep; i++) {
			double theta = 2.0 * M_PI * i / num_ep;
			double radius = parameters.doubles("membrane_circle_radius");
			double x = radius * cos(theta);
			double y = radius * sin(theta);

			// compute radial distance and unit‐radial direction
			double r_norm = std::sqrt(x*x + y*y);
			double nx = x / r_norm;    // outward radial unit vector
			double ny = y / r_norm;

			// step back along the normal by ep_dis
			double xi = x - ep_dis * nx;
			double yi = y - ep_dis * ny;

			Cell* pC = create_cell( *pBM_def );

			// Turn on proliferation for first cell

			if( parameters.ints("number_EP_cells") > 1 ){
				pC->assign_position( { xi, yi, 0.0 } );
			}
			else{
				pC->assign_position( { parameters.doubles("x"), parameters.doubles("y"), 0.0 } );
			}
			
			if (i==0){
				std::cout << "Setting first cell to proliferate" << std::endl;
				pC->phenotype.cycle.data.exit_rate(0) = parameters.doubles("proliferation_exit_rate");
			}
			
		}
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
