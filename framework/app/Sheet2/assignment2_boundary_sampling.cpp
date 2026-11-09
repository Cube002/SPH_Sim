#include <algorithm>
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
	const double particle_radius = 0.0025;
	const double particle_diameter = 2.0 * particle_radius;
	const double boundary_sampling_distance = 0.8 * particle_diameter;
	const double smoothing_length = 1.2 * particle_diameter;
	const double compact_support = 2.0 * smoothing_length;
	const double boundary_density = 1000.0;

	const std::vector<learnSPH::TriMesh> meshes =
		learnSPH::read_tri_meshes_from_obj("../obj/IcoSphere.obj");
	if (meshes.empty())
	{
		std::cerr << "Could not load ../obj/IcoSphere.obj" << std::endl;
		return 1;
	}
	const learnSPH::TriMesh &mesh = meshes[0];

	std::vector<Eigen::Vector3d> boundary_particles;
	learnSPH::sampling::triangle_mesh(
		boundary_particles,
		mesh.vertices,
		mesh.triangles,
		boundary_sampling_distance);

	CompactNSearch::NeighborhoodSearch nsearch(compact_support);
	const unsigned int boundary_id = nsearch.add_point_set(
		boundary_particles[0].data(), boundary_particles.size(), false);
	nsearch.find_neighbors();

	std::vector<double> boundary_volumes;
	std::vector<double> boundary_masses;
	learnSPH::sph_math::compute_boundary_volumes(
		boundary_particles, boundary_volumes, nsearch, boundary_id,
		smoothing_length);
	learnSPH::sph_math::compute_boundary_masses(
		boundary_volumes, boundary_masses, boundary_density);

	learnSPH::write_particles_to_vtk(
		"../res/assignment2_boundary_volume.vtk",
		boundary_particles, boundary_volumes, "volume");
	learnSPH::write_particles_to_vtk(
		"../res/assignment2_boundary_mass.vtk",
		boundary_particles, boundary_masses, "mass");

	const auto volume_range = std::minmax_element(
		boundary_volumes.begin(), boundary_volumes.end());
	const auto mass_range = std::minmax_element(
		boundary_masses.begin(), boundary_masses.end());

	std::cout << "Assignment 2: boundary sampling" << std::endl;
	std::cout << "  boundary particles: " << boundary_particles.size()
		<< std::endl;
	std::cout << "  volume range: [" << *volume_range.first
		<< ", " << *volume_range.second << "] m^3" << std::endl;
	std::cout << "  mass range: [" << *mass_range.first
		<< ", " << *mass_range.second << "] kg" << std::endl;
	std::cout << "  wrote assignment2_boundary_volume.vtk and "
		<< "assignment2_boundary_mass.vtk to `<build_folder>/res`."
		<< std::endl;

	return 0;
}
