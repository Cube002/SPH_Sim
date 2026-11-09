#include <stdlib.h>     // rand
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>    // std::max
#include <random>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

#include "../../learnSPH/io.h"
#include "../../learnSPH/sampling.h"

int main()
{
	std::cout << "Welcome to the learnSPH framework!!" << std::endl;
	std::cout << "Generating a sample scene...";

	// Load a obj surface mesh
	const std::vector<learnSPH::TriMesh> meshes = learnSPH::read_tri_meshes_from_obj("../res/box.obj");
	const learnSPH::TriMesh& box = meshes[0];

	// Sample the mesh with particles
	const double sampling_distance = 0.4;
	const double r = 0.05;
	std::vector<Eigen::Vector3d> particles;	
	learnSPH::sampling::triangle_mesh(particles, box.vertices, box.triangles, sampling_distance);
	std::vector<Eigen::Vector3d> particles_velocity(particles.size(), Eigen::Vector3d::Zero());

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist_vel(-2.0, 2.0);
	for (size_t i = 0; i < particles_velocity.size(); ++i) {
		particles_velocity[i] = Eigen::Vector3d(dist_vel(gen), dist_vel(gen), dist_vel(gen));
	}


	// Initialize data vectors
	std::vector<double> particles_scalar_data(particles.size());
	std::vector<Eigen::Vector3d> particles_vector_data(particles.size());
	

	// Scalar data will be the particle_id
	for (int i = 0; i < (int)particles.size(); i++) {
		particles_scalar_data[i] = i;
		particles[i] += Eigen::Vector3d(0,0,10);
	}

	const double search_radius = 5;
	CompactNSearch::NeighborhoodSearch nsearch(search_radius);
	if (!particles.empty()) {
		nsearch.add_point_set(particles[0].data(), particles.size(), true);
		nsearch.find_neighbors();
	}

	// Simulation loop
	for (int time_step = 0; time_step < 1000; time_step++) {
		double deltatime = 0.01;
		const Eigen::Vector3d g = Eigen::Vector3d(0, 0, -9.81);
		if (!particles.empty()) {
			nsearch.find_neighbors(true);
		}
		const auto &ps = nsearch.point_set(0);
		for (int i = 0; i < (int)particles.size(); i++) {
			
			//Force:
			//Gravity sim
			Eigen::Vector3d F = Eigen::Vector3d(0,0,0);
			const unsigned int num_n = static_cast<unsigned int>(ps.n_neighbors(0, static_cast<unsigned int>(i)));
			for (unsigned int k = 0; k < num_n; ++k) {
				const unsigned int j = ps.neighbor(0, static_cast<unsigned int>(i), k);
				if (j == static_cast<unsigned int>(i)) continue;
				const Eigen::Vector3d d = particles[j] - particles[i];
				const double len2 = d.squaredNorm();
				if (len2 < 1e-3) continue;
				F += d.normalized()*0.01 / len2;
			}


			particles_velocity[i] += deltatime * F;
			particles_velocity[i] = particles_velocity[i]*0.99;
			particles[i] += deltatime * particles_velocity[i];
			particles_vector_data[i] = particles_velocity[i]; // geschwindigkeiten werden gespeichert.
			// if(std::abs(particles[i][2]) <= r){
			// 	particles[i][2] = 0;
			// 	particles_velocity[i][2] = - particles_velocity[i][2];
			// }

		}

		if (particles.size() >= 2 && (time_step == 0 || time_step == 50 || time_step == 100 || time_step == 200)) {
			Eigen::Vector3d d = particles[0] - particles[1];
			std::cout << "step " << time_step << " dist(p0,p1)=" << d.norm() << std::endl;
		}

		// Save output
		const std::string filename = "../res/example_" + std::to_string(time_step) + ".vtk";
		learnSPH::write_particles_to_vtk(filename, particles, particles_scalar_data, particles_vector_data);
	}

	std::cout << "completed!" << std::endl;
	std::cout << "The scene files have been saved in the folder `<build_folder>/res`. \nYou can visualize them with Blender by using the Blender Sequence Loaded addon." << std::endl;

	return 0;
}