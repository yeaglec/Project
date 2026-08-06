#include "./Tests.h" 
#include "./Membrane.h"
#include "./Geometry.h"
#include "./Utils.h" 
#include <cmath>
#include <cfloat>
#include <iomanip>
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

void Test_Bending() {
    std::cout << "Initializing Bending Smoothing Test..." << std::endl;
    boundary_membrane_pts.clear();

    // Closed ring test: apply force to one node only and watch bending spread it.
    const int Np = 25;
    const double radius = 120.0;
    for (int i = 0; i < Np; ++i) {
        double theta = (2.0 * M_PI * i) / Np;
        boundary_membrane_pts.push_back({radius * std::cos(theta), radius * std::sin(theta), 0.0});
    }

    initial_edge_length.resize(Np);
    initial_node_positions.resize(Np);

    for (int i = 0; i < Np; ++i) {
        initial_node_positions[i] = boundary_membrane_pts[i];
        int next = (i + 1) % Np;
        double dx = boundary_membrane_pts[next][0] - boundary_membrane_pts[i][0];
        double dy = boundary_membrane_pts[next][1] - boundary_membrane_pts[i][1];
        initial_edge_length[i] = std::sqrt(dx * dx + dy * dy);
    }

}

void Test_Bending_Pert(std::vector<std::pair<double,double>>& node_forces, double current_time) {
    int Np = (int)boundary_membrane_pts.size();
    if (Np < 5) return;

    const int drive_idx = 0;
    const double pulse_duration = 250.0;
    double pulse_force = parameters.doubles("test_pull_force");
    if (pulse_force <= 0.0) pulse_force = 2.5;

    // Point load on a single node, analogous to unsmoothed force transfer.
    if (current_time < pulse_duration) {
        double x = boundary_membrane_pts[drive_idx][0];
        double y = boundary_membrane_pts[drive_idx][1];
        double norm = std::sqrt(x * x + y * y);

        if (norm > 1e-16) {
            node_forces[drive_idx].first  += pulse_force * (x / norm);
            node_forces[drive_idx].second += pulse_force * (y / norm);
        } else {
            node_forces[drive_idx].first += pulse_force;
        }
    }

    auto node_radius = [&](int idx) -> double {
        double x = boundary_membrane_pts[idx][0];
        double y = boundary_membrane_pts[idx][1];
        return std::sqrt(x * x + y * y);
    };

    auto node_kink = [&](int idx) -> double {
        int prev = (idx - 1 + Np) % Np;
        int next = (idx + 1) % Np;
        double d2x = boundary_membrane_pts[prev][0] + boundary_membrane_pts[next][0] - 2.0 * boundary_membrane_pts[idx][0];
        double d2y = boundary_membrane_pts[prev][1] + boundary_membrane_pts[next][1] - 2.0 * boundary_membrane_pts[idx][1];
        return std::sqrt(d2x * d2x + d2y * d2y);
    };

    static bool header_printed = false;
    static double next_log_time = 0.0;
    if (current_time >= next_log_time) {
        int left1 = (drive_idx - 1 + Np) % Np;
        int right1 = (drive_idx + 1) % Np;
        int left2 = (drive_idx - 2 + Np) % Np;
        int right2 = (drive_idx + 2) % Np;

        double drive_radius = node_radius(drive_idx);
        double adjacent_mean_radius = 0.5 * (node_radius(left1) + node_radius(right1));
        double second_neighbor_mean_radius = 0.5 * (node_radius(left2) + node_radius(right2));
        double drive_kink = node_kink(drive_idx);

        double mean_kink = 0.0;
        for (int i = 0; i < Np; ++i) {
            mean_kink += node_kink(i);
        }
        mean_kink /= (double)Np;

        if (!header_printed) {
            std::cout << "[Bending Smoothing Test] One-node pulse load, then relaxation." << std::endl;
            std::cout << "Columns: time, drive radius, neighbors radius, kink at drive node, mean kink." << std::endl;
            header_printed = true;
        }

        std::cout << std::fixed << std::setprecision(2)
                  << "[Bending Smoothing Test] t=" << current_time
                  << " | drive_r=" << drive_radius
                  << " | adj_r=" << adjacent_mean_radius
                  << " | adj2_r=" << second_neighbor_mean_radius
                  << " | drive_kink=" << drive_kink
                  << " | mean_kink=" << mean_kink
                  << " | pulse=" << (current_time < pulse_duration ? "ON" : "OFF")
                  << std::endl;

        next_log_time += 25.0;
    }
}

void Test_Remesh() {
    std::cout << "Initializing Remesh Refinement Test..." << std::endl;
    boundary_membrane_pts.clear();

    // Centered near the single setup CAF so cell-aware remeshing has local support.
    const double cx = 0;
    const double cy = 0;
    const int Np = 12;
    const double radius =200 ;

    for (int i = 0; i < Np; i++) {
        double theta = (2.0 * M_PI * i) / Np;
        double x = cx + radius * std::cos(theta);
        double y = cy + radius * std::sin(theta);
        boundary_membrane_pts.push_back({x, y, 0.0});
    }

    std::cout << "[Remesh Test] Local pull should stretch one region and trigger local node insertion." << std::endl;
}

void Test_Remesh_Pert(std::vector<std::pair<double,double>>& node_forces, double current_time) {
    int Np = (int)boundary_membrane_pts.size();
    if (Np < 4) return;

    const double cx = 0;
    const double cy = 0;
    const int drive_idx = 0;
    const double pull_duration = 350.0;

    double pull_force = parameters.doubles("test_pull_force");
    if (pull_force <= 0.0) pull_force = 3.0;

    // Localized pull: stretches one zone so remeshing is easier to interpret.
    if (current_time < pull_duration) {
        int left = (drive_idx - 1 + Np) % Np;
        int right = (drive_idx + 1) % Np;

        auto add_radial_pull = [&](int idx, double scale) {
            double dx = boundary_membrane_pts[idx][0] - cx;
            double dy = boundary_membrane_pts[idx][1] - cy;
            double norm = std::sqrt(dx * dx + dy * dy);
            if (norm <= 1e-16) return;

            node_forces[idx].first += scale * pull_force * (dx / norm);
            node_forces[idx].second += scale * pull_force * (dy / norm);
        };

        add_radial_pull(drive_idx, 1.0);
        add_radial_pull(left, 0.45);
        add_radial_pull(right, 0.45);
    }

    auto edge_stats = [&]() {
        int n = (int)boundary_membrane_pts.size();
        double max_edge = 0.0;
        double mean_edge = 0.0;

        for (int i = 0; i < n; ++i) {
            int j = (i + 1) % n;
            double dx = boundary_membrane_pts[j][0] - boundary_membrane_pts[i][0];
            double dy = boundary_membrane_pts[j][1] - boundary_membrane_pts[i][1];
            double len = std::sqrt(dx * dx + dy * dy);
            mean_edge += len;
            if (len > max_edge) max_edge = len;
        }

        mean_edge /= (double)n;
        return std::pair<double, double>(max_edge, mean_edge);
    };

    static int last_node_count = -1;
    if (last_node_count < 0) {
        last_node_count = (int)boundary_membrane_pts.size();
    }

    int current_nodes = (int)boundary_membrane_pts.size();
    if (current_nodes > last_node_count) {
        std::cout << "[Remesh Event] t=" << current_time
                  << " | Nodes " << last_node_count << " -> " << current_nodes
                  << std::endl;
        last_node_count = current_nodes;
    }

    static double next_log_time = 0.0;
    if (current_time >= next_log_time) {
        auto [max_edge, mean_edge] = edge_stats();

        std::cout << std::fixed << std::setprecision(2)
                  << "[Remesh Test] t=" << current_time
                  << " | nodes=" << current_nodes
                  << " | max_edge=" << max_edge
                  << " | mean_edge=" << mean_edge
                  << " | threshold=" << parameters.doubles("max_edge_length")
                  << " | stretch_factor=" << parameters.doubles("remesh_stretch_factor")
                  << " | pull=" << (current_time < pull_duration ? "ON" : "OFF")
                  << std::endl;

        next_log_time += 40.0;
    }
}