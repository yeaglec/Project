#include "./Tests.h" 
#include "./Utils.h" 
#include "./Models.h" 
#include "./custom.h"
#include <cmath>
#include <cfloat>
#include <iostream>

void Test_Line(){


    boundary_membrane_pts.push_back({0.0, 50, 0.0}); 
	boundary_membrane_pts.push_back({0.0, -50, 0.0});

	int Np = (int)boundary_membrane_pts.size(); // Should be 2
	std::cout << "Initial Boundary Points: " << Np << std::endl;
	initial_edge_length.resize(Np);
	initial_node_positions.resize(Np);

}

void Test_Line_Pert(std::vector<std::pair<double,double>>& node_forces, double current_time) {
    
    if (current_time < 500.0) {
        node_forces[0].first += 1.0; 
        node_forces[1].first += 1.0; 
    }

    // Do Logging here
    static double next_log_time = 0.0;
    if (current_time >= next_log_time) {
        double current_y = boundary_membrane_pts[1][1];
        double home_y = initial_node_positions[1][1];
        double displacement = current_y - home_y;
        
        std::cout << "Time: " << current_time 
                  << " Node YPos: " << current_y 
                  << " Displacement: " << displacement 
                  << " Net Velocity Y: " << node_forces[1].second << std::endl;
        
        next_log_time += 100.0;
    }
}


void Test_Tri() {

    // Clear any existing points just in case
    boundary_membrane_pts.clear();
    
    boundary_membrane_pts.push_back({-100.0, 0.0, 0.0}); 
    boundary_membrane_pts.push_back({100.0, 0.0, 0.0});
    boundary_membrane_pts.push_back({0.0, 100.0, 0.0});

    int Np = (int)boundary_membrane_pts.size();     
    initial_edge_length.resize(Np);
    initial_node_positions.resize(Np);

}

void Test_Tri_Pert(std::vector<std::pair<double,double>>& node_forces, double current_time) {
    
    if (current_time < 500.0) {
        node_forces[2].second += .5; 
    }

    // Do Logging here
    static double next_log_time = 0.0;
    if (current_time >= next_log_time) {
        double current_y = boundary_membrane_pts[2][1];
        double home_y = initial_node_positions[2][1];
        double displacement = current_y - home_y;
        
        std::cout << "Time: " << current_time 
                  << " Node YPos: " << current_y 
                  << " Displacement: " << displacement 
                  << " Net Velocity Y: " << node_forces[2].second << std::endl;
        
        next_log_time += 100.0;
    }
}

void Test_Square() {
    std::cout << "Initializing Square Test..." << std::endl;
    boundary_membrane_pts.clear();
    
    // 8-node square (Corners and midpoints)
    boundary_membrane_pts.push_back({-50.0, -50.0, 0.0}); // Bottom Left (0)
    boundary_membrane_pts.push_back({  0.0, -50.0, 0.0}); // Bottom Mid (1)
    boundary_membrane_pts.push_back({ 50.0, -50.0, 0.0}); // Bottom Right (2)
    boundary_membrane_pts.push_back({ 50.0,   0.0, 0.0}); // Right Mid (3)
    boundary_membrane_pts.push_back({ 50.0,  50.0, 0.0}); // Top Right (4)
    boundary_membrane_pts.push_back({  0.0,  50.0, 0.0}); // Top Mid (5)
    boundary_membrane_pts.push_back({-50.0,  50.0, 0.0}); // Top Left (6)
    boundary_membrane_pts.push_back({-50.0,   0.0, 0.0}); // Left Mid (7)
}

void Test_Square_Pert(std::vector<std::pair<double,double>>& node_forces, double current_time) {
    // Pull the Top Mid node (index 5) upwards for the first 400 minutes
    if (current_time < 400.0) {
        node_forces[5].second += 2.0; 
    }

    static double next_log_time = 0.0;
    if (current_time >= next_log_time) {
        double current_y = boundary_membrane_pts[5][1];
        double corner_x = boundary_membrane_pts[4][0]; // Watch the top right corner X drift
        
        std::cout << "[Square Test] Time: " << current_time 
                  << " | Top Node Y: " << current_y 
                  << " | Corner Node X: " << corner_x 
                  << std::endl;
        
        next_log_time += 50.0;
    }
}

void Test_Ring() {
    std::cout << "Initializing Ring Test..." << std::endl;
    boundary_membrane_pts.clear();
    
    int Np = 16;
    double radius = 100.0;
    
    for (int i = 0; i < Np; i++) {
        double theta = (2.0 * M_PI * i) / Np;
        double x = radius * std::cos(theta);
        double y = radius * std::sin(theta);
        boundary_membrane_pts.push_back({x, y, 0.0});
    }
}

void Test_Ring_Pert(std::vector<std::pair<double,double>>& node_forces, double current_time) {
    // Apply an outward velocity to the far-right node (index 0) and its two neighbors (15 and 1)
    // to simulate localized tumor/cell proliferation pushing the boundary.
    if (current_time < 600.0) {
        node_forces[0].first  += 3.0; // Pushing right
        node_forces[1].first  += 2.0; // Pushing mostly right, slightly up
        node_forces[1].second += 1.0; 
        node_forces[15].first += 2.0; // Pushing mostly right, slightly down
        node_forces[15].second -= 1.0; 
    }

    static double next_log_time = 0.0;
    if (current_time >= next_log_time) {
        double bulge_x = boundary_membrane_pts[0][0];
        
        std::cout << "[Ring Test] Time: " << current_time 
                  << " | Bulge Node X-Pos (Ideal = 100): " << bulge_x 
                  << std::endl;
        
        next_log_time += 100.0;
    }
}