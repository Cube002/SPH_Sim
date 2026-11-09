#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

#include "../../learnSPH/io.h"
#include "../../learnSPH/sampling.h"
#include "../../learnSPH/sph_math.h"

int main()
{
	const double particle_radius = 0.0025;
	const double particle_diameter = 2.0 * particle_radius;
	const double fluid_sampling_distance = particle_diameter;
	const double boundary_sampling_distance = 0.8 * particle_diameter;
	const double smoothing_length = 1.2 * particle_diameter;
	const double compact_support = 2.0 * smoothing_length;

	const double rest_density = 1000.0;
	const double fluid_particle_mass = rest_density
		* fluid_sampling_distance
		* fluid_sampling_distance
		* fluid_sampling_distance;

	const int particles_per_axis = 20;
	const double box_side = particles_per_axis * fluid_sampling_distance;
	const Eigen::Vector3d fluid_bottom = Eigen::Vector3d::Zero();
	const Eigen::Vector3d fluid_top = Eigen::Vector3d::Constant(box_side);

	std::vector<Eigen::Vector3d> fluid_particles;
	learnSPH::sampling::fluid_box(
		fluid_particles, fluid_bottom, fluid_top,
		fluid_sampling_distance);

	const double floor_z =
		0.5 * fluid_sampling_distance - 2.0 * particle_radius;
	const double floor_margin = boundary_sampling_distance;
	const double floor_min = -floor_margin;
	const double floor_max = box_side + floor_margin;
	const int floor_particles_per_axis = static_cast<int>(std::round(
		(floor_max - floor_min) / boundary_sampling_distance)) + 1;

	std::vector<Eigen::Vector3d> boundary_particles;
	boundary_particles.reserve(
		static_cast<std::size_t>(floor_particles_per_axis)
		* floor_particles_per_axis);
	for (int i = 0; i < floor_particles_per_axis; ++i)
	{
		const double x = floor_min + i * boundary_sampling_distance;
		for (int j = 0; j < floor_particles_per_axis; ++j)
		{
			const double y = floor_min + j * boundary_sampling_distance;
			boundary_particles.push_back(Eigen::Vector3d(x, y, floor_z));
		}
	}

	CompactNSearch::NeighborhoodSearch nsearch(compact_support);
	const unsigned int boundary_id = nsearch.add_point_set(
		boundary_particles[0].data(), boundary_particles.size(), true);
	const unsigned int fluid_id = nsearch.add_point_set(
		fluid_particles[0].data(), fluid_particles.size(), true);
	nsearch.find_neighbors(true);

	std::vector<double> boundary_volumes;
	std::vector<double> boundary_masses;
	learnSPH::sph_math::compute_boundary_volumes(
		boundary_particles, boundary_volumes, nsearch, boundary_id,
		smoothing_length);
	learnSPH::sph_math::compute_boundary_masses(
		boundary_volumes, boundary_masses, rest_density);

	std::vector<double> fluid_densities;
	learnSPH::sph_math::compute_density(
		fluid_particles, boundary_particles, fluid_densities,
		nsearch, fluid_id, boundary_id, fluid_particle_mass,
		boundary_volumes, rest_density, smoothing_length);

	learnSPH::write_particles_to_vtk(
		"../res/assignment3_fluid_density.vtk",
		fluid_particles, fluid_densities, "density");
	learnSPH::write_particles_to_vtk(
		"../res/assignment3_floor_boundary.vtk",
		boundary_particles, boundary_volumes, "volume");

	const auto density_range = std::minmax_element(
		fluid_densities.begin(), fluid_densities.end());
	const double density_mean = std::accumulate(
		fluid_densities.begin(), fluid_densities.end(), 0.0)
		/ static_cast<double>(fluid_densities.size());
	const std::size_t close_to_rest = static_cast<std::size_t>(std::count_if(
		fluid_densities.begin(), fluid_densities.end(),
		[rest_density](double density)
		{
			return std::abs(density - rest_density) / rest_density < 0.01;
		}));

	const auto &fluid_point_set = nsearch.point_set(fluid_id);
	std::size_t fluid_with_boundary_neighbors = 0;
	for (unsigned int i = 0; i < fluid_particles.size(); ++i)
	{
		if (fluid_point_set.n_neighbors(boundary_id, i) > 0)
		{
			++fluid_with_boundary_neighbors;
		}
	}

	std::cout << "Assignment 3: fluid density with boundary contribution"
		<< std::endl;
	std::cout << "  fluid particles: " << fluid_particles.size() << std::endl;
	std::cout << "  boundary particles: " << boundary_particles.size()
		<< std::endl;
	std::cout << "  density range: [" << *density_range.first
		<< ", " << *density_range.second << "] kg/m^3" << std::endl;
	std::cout << "  mean density: " << density_mean << " kg/m^3" << std::endl;
	std::cout << "  particles within 1% of rest density: "
		<< close_to_rest << " / " << fluid_densities.size() << std::endl;
	std::cout << "  fluid particles with boundary neighbors: "
		<< fluid_with_boundary_neighbors << std::endl;
	std::cout << "  wrote assignment3_fluid_density.vtk and "
		<< "assignment3_floor_boundary.vtk to `<build_folder>/res`."
		<< std::endl;

	return 0;
}
