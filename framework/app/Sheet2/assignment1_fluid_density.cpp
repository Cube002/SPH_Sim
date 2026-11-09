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
	const double smoothing_length = 1.2 * particle_diameter;
	const double compact_support = 2.0 * smoothing_length;

	const double rest_density = 1000.0;
	const double particle_mass = rest_density
		* fluid_sampling_distance
		* fluid_sampling_distance
		* fluid_sampling_distance;

	const int particles_per_axis = 20;
	const double box_side = particles_per_axis * fluid_sampling_distance;

	std::vector<Eigen::Vector3d> fluid_particles;
	learnSPH::sampling::fluid_box(
		fluid_particles,
		Eigen::Vector3d::Zero(),
		Eigen::Vector3d::Constant(box_side),
		fluid_sampling_distance);

	CompactNSearch::NeighborhoodSearch nsearch(compact_support);
	const unsigned int fluid_id = nsearch.add_point_set(
		fluid_particles[0].data(), fluid_particles.size(), false);
	nsearch.find_neighbors();

	std::vector<double> fluid_densities;
	learnSPH::sph_math::compute_density(
		fluid_particles, fluid_densities, nsearch, fluid_id,
		particle_mass, smoothing_length);

	learnSPH::write_particles_to_vtk(
		"../res/assignment1_fluid_density.vtk",
		fluid_particles, fluid_densities, "density");

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

	std::cout << "Assignment 1: fluid-box density" << std::endl;
	std::cout << "  particles: " << fluid_particles.size() << std::endl;
	std::cout << "  density range: [" << *density_range.first
		<< ", " << *density_range.second << "] kg/m^3" << std::endl;
	std::cout << "  mean density: " << density_mean << " kg/m^3" << std::endl;
	std::cout << "  particles within 1% of rest density: "
		<< close_to_rest << " / " << fluid_densities.size() << std::endl;
	std::cout << "  wrote assignment1_fluid_density.vtk to `<build_folder>/res`."
		<< std::endl;

	return 0;
}
