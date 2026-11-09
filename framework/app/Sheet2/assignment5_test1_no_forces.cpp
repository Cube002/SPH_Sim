#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

#include "../../learnSPH/io.h"
#include "../../learnSPH/sampling.h"
#include "../../learnSPH/sph_math.h"

int main()
{
	// Assignment 5, first baseline test case:
	// 1 m^3 fluid, 3 cm particle radius, B = 0, no viscosity, no gravity.
	const double particle_radius = 0.03;
	const double fluid_sampling_distance = 2.0 * particle_radius;
	const double smoothing_length = 1.2 * fluid_sampling_distance;
	const double compact_support = 2.0 * smoothing_length;

	const double rest_density = 1000.0;
	const double fluid_particle_mass = rest_density
		* fluid_sampling_distance
		* fluid_sampling_distance
		* fluid_sampling_distance;

	const double stiffness = 0.0;
	const double fluid_viscosity = 0.0;
	const double boundary_viscosity = 0.0;
	const Eigen::Vector3d gravity = Eigen::Vector3d::Zero();

	const double max_time_step = 0.001;
	const double simulation_time = 0.1;
	const double frame_interval = 0.02;
	const int num_time_steps = static_cast<int>(
		std::round(simulation_time / max_time_step));

	std::vector<Eigen::Vector3d> fluid_particles;
	learnSPH::sampling::fluid_box(
		fluid_particles,
		Eigen::Vector3d::Zero(),
		Eigen::Vector3d::Ones(),
		fluid_sampling_distance);

	std::vector<Eigen::Vector3d> fluid_velocities(
		fluid_particles.size(), Eigen::Vector3d::Zero());
	std::vector<Eigen::Vector3d> fluid_accelerations(
		fluid_particles.size(), Eigen::Vector3d::Zero());

	// No boundary point set is registered in this first test case. Later
	// Assignment 5 scenes can add one here and use the boundary-aware force
	// overloads, just like assignment4_forces_demo.cpp.
	CompactNSearch::NeighborhoodSearch nsearch(compact_support);
	const unsigned int fluid_id = nsearch.add_point_set(
		fluid_particles[0].data(), fluid_particles.size(), true);

	std::cout << "Assignment 5 baseline test" << std::endl;
	std::cout << "Fluid particles: " << fluid_particles.size() << std::endl;
	std::cout << "B = " << stiffness
		<< ", nu_f = " << fluid_viscosity
		<< ", nu_b = " << boundary_viscosity << std::endl;
	std::cout << "gravity = (" << gravity[0] << ", " << gravity[1]
		<< ", " << gravity[2] << ") m/s^2" << std::endl;
	std::cout << "dt = " << max_time_step << " s, steps = "
		<< num_time_steps << std::endl;

	std::vector<double> fluid_densities;
	std::vector<double> fluid_pressures;
	double time = 0.0;
	double next_frame_time = 0.0;
	int frame_index = 0;

	for (int time_step = 0; time_step < num_time_steps; ++time_step)
	{
		nsearch.find_neighbors(true);

		learnSPH::sph_math::compute_density(
			fluid_particles, fluid_densities, nsearch, fluid_id,
			fluid_particle_mass, smoothing_length);
		learnSPH::sph_math::compute_pressure(
			fluid_densities, fluid_pressures, rest_density, stiffness);

		fluid_accelerations.assign(
			fluid_particles.size(), Eigen::Vector3d::Zero());
		learnSPH::sph_math::compute_pressure_acceleration(
			fluid_particles, fluid_densities, fluid_pressures,
			fluid_accelerations, nsearch, fluid_id,
			fluid_particle_mass, smoothing_length);
		learnSPH::sph_math::compute_friction_acceleration(
			fluid_particles, fluid_velocities, fluid_densities,
			fluid_accelerations, nsearch, fluid_id,
			fluid_particle_mass, fluid_viscosity, smoothing_length);
		learnSPH::sph_math::compute_body_acceleration(
			fluid_accelerations, gravity);

		// Semi-implicit Euler time integration.
		for (unsigned int i = 0; i < fluid_particles.size(); ++i)
		{
			fluid_velocities[i] += max_time_step * fluid_accelerations[i];
			fluid_particles[i] += max_time_step * fluid_velocities[i];
		}

		time += max_time_step;

		// Export positions, density and velocity for visualization.
		if (time + 1e-12 >= next_frame_time)
		{
			const std::string filename = "../res/assignment5_test1_no_forces_"
				+ std::to_string(frame_index) + ".vtk";
			learnSPH::write_particles_to_vtk(
				filename, fluid_particles, fluid_densities,
				fluid_velocities);
			++frame_index;
			next_frame_time += frame_interval;
		}

		// TODO Assignment 5: replace the fixed step with
		// dt = min(dt_CFL, dt_default), add the second fluid block, enable
		// gravity, and add a boundary floor with boundary viscosity.
		// TODO Moving boundaries: update boundary positions and velocities,
		// refresh the point-set positions, and pass boundary_velocities to the
		// boundary-aware friction overload.
	}

	std::cout << "Simulation finished at t = " << time << " s" << std::endl;
	std::cout << "Wrote " << frame_index
		<< " VTK frames to `<build_folder>/res`." << std::endl;

	return 0;
}
