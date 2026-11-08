#include <iostream>
#include <vector>
#include <string>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

#include "../learnSPH/io.h"
#include "../learnSPH/sampling.h"
#include "../learnSPH/sph_math.h"

int main()
{

	//Values defined:
	const double particle_radius = 0.0025;
	const double particle_diameter = 2.0 * particle_radius;
	const double fluid_sampling_distance = particle_diameter;
	const double smoothing_length = 1.2 * particle_diameter;
	const double compact_support = 2.0 * smoothing_length;

	const double fluid_density = 1000.0;
	const double particle_mass = fluid_density
		* fluid_sampling_distance
		* fluid_sampling_distance
		* fluid_sampling_distance;

	//Define Fluid Domain:
	const int n_per_axis = 20;
	const double box_size = (n_per_axis - 0.5) * fluid_sampling_distance;
	std::vector<Eigen::Vector3d> particles;
	learnSPH::sampling::fluid_box(
		particles,
		Eigen::Vector3d::Zero(),
		Eigen::Vector3d::Constant(box_size),
		fluid_sampling_distance);
	std::vector<Eigen::Vector3d> boundary_particles;//none at the moment
	std::vector<Eigen::Vector3d> particles_velocity(particles.size(), Eigen::Vector3d::Zero());

	//Prepare search module:
	CompactNSearch::NeighborhoodSearch nsearch(compact_support);
	const unsigned int fluid_id = nsearch.add_point_set(
		particles[0].data(), particles.size(), true);

	std::vector<double> densities;
	std::vector<double> masses;

	const int num_time_steps = 100;
	const double dt = 0.01;
	double time = 0.0;
	//Time Loop
	for (int time_step = 0; time_step < num_time_steps; ++time_step)
	{
		nsearch.find_neighbors(true);
		learnSPH::sph_math::compute_density(
			particles, densities, nsearch, fluid_id,
			particle_mass, smoothing_length);
		// learnSPH::sph_math::compute_compensated_mass(
		// 	densities, masses, particle_mass, fluid_density);
		std::vector<Eigen::Vector3d> particle_accelerations(particles.size(),Eigen::Vector3d::Zero());
		learnSPH::sph_math::compute_pressure_acceleration(particles,densities,particle_accelerations,nsearch,fluid_id,particle_mass,fluid_density,3,smoothing_length);
		for(size_t i = 0; i < particles.size(); i++)
		{
			
			
			particles_velocity[i] += dt * particle_accelerations[i];
			particles[i]+= dt * particles_velocity[i];
		}
		
		
		const std::string filename = "../res/masses_" + std::to_string(time_step) + ".vtk";
		learnSPH::write_particles_to_vtk(filename, particles, masses);

		time += dt;
	}

	std::cout << "Simulation finished at t = " << time << " s" << std::endl;
	std::cout << "The scene files have been saved in the folder `<build_folder>/res`." << std::endl;

	return 0;
}
