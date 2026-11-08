#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

#include "../learnSPH/io.h"
#include "../learnSPH/sampling.h"
#include "../learnSPH/sph_math.h"

int main()
{
	// Particle size relations used throughout the course.
	const double particle_radius = 0.0025;
	const double particle_diameter = 2.0 * particle_radius;
	const double fluid_sampling_distance = particle_diameter;
	const double boundary_sampling_distance = 0.8 * particle_diameter;
	const double smoothing_length = 1.2 * particle_diameter;
	const double compact_support = 2.0 * smoothing_length;

	const double fluid_density = 1000.0;
	const double boundary_density = fluid_density;
	const double fluid_particle_mass = fluid_density
		* fluid_sampling_distance
		* fluid_sampling_distance
		* fluid_sampling_distance;

	// Load the boundary surface from the provided OBJ file.
	const std::vector<learnSPH::TriMesh> meshes =
		learnSPH::read_tri_meshes_from_obj("../obj/IcoSphere.obj");
	if (meshes.empty())
	{
		std::cerr << "Could not load ../obj/IcoSphere.obj" << std::endl;
		return 1;
	}
	const learnSPH::TriMesh &sphere = meshes[0];

	// Sample the sphere surface with boundary particles.
	std::vector<Eigen::Vector3d> boundary_particles;
	learnSPH::sampling::triangle_mesh(
		boundary_particles,
		sphere.vertices,
		sphere.triangles,
		boundary_sampling_distance);

	double sphere_top = sphere.vertices[0][2];
	for (const Eigen::Vector3d &vertex : sphere.vertices)
	{
		sphere_top = std::max(sphere_top, vertex[2]);
	}

	// Sample a 20x20x20 fluid box just above the sphere. The vertical gap is
	// one particle diameter, so the lowest fluid particles barely touch the
	// highest boundary particles.
	const int n_per_axis = 20;
	const double half_width =
		0.5 * (n_per_axis - 1) * fluid_sampling_distance;
	const double box_size =
		(n_per_axis - 0.5) * fluid_sampling_distance;

	const Eigen::Vector3d fluid_bottom(
		-half_width,
		-half_width,
		sphere_top + particle_diameter);
	const Eigen::Vector3d fluid_top =
		fluid_bottom + Eigen::Vector3d::Constant(box_size);

	std::vector<Eigen::Vector3d> fluid_particles;
	learnSPH::sampling::fluid_box(
		fluid_particles,
		fluid_bottom,
		fluid_top,
		fluid_sampling_distance);

	// Register both point sets with CompactNSearch. Cross-set neighbors are
	// enabled by default, which lets fluid particles find boundary particles.
	CompactNSearch::NeighborhoodSearch nsearch(compact_support);
	const unsigned int boundary_id = nsearch.add_point_set(
		boundary_particles[0].data(), boundary_particles.size(), true);
	const unsigned int fluid_id = nsearch.add_point_set(
		fluid_particles[0].data(), fluid_particles.size(), true);
	nsearch.find_neighbors(true);

	// Assignment 2: compensate the boundary volume and mass for uneven sampling.
	std::vector<double> boundary_volumes;
	std::vector<double> boundary_masses;
	learnSPH::sph_math::compute_boundary_volumes(
		boundary_particles, boundary_volumes, nsearch, boundary_id,
		smoothing_length);
	learnSPH::sph_math::compute_boundary_masses(
		boundary_volumes, boundary_masses, boundary_density);

	// Assignment 3: fluid density including the boundary contribution.
	std::vector<double> fluid_densities;
	learnSPH::sph_math::compute_density(
		fluid_particles, boundary_particles, fluid_densities,
		nsearch, fluid_id, boundary_id, fluid_particle_mass,
		boundary_volumes, fluid_density, smoothing_length);

	// Assignment 4 buffers. Pressure and acceleration computation will be added
	// here; the current task only initializes the rigid setup.
	std::vector<Eigen::Vector3d> fluid_velocities(
		fluid_particles.size(), Eigen::Vector3d::Zero());
	std::vector<double> fluid_pressures(fluid_particles.size(), 0.0);
	std::vector<Eigen::Vector3d> fluid_accelerations(
		fluid_particles.size(), Eigen::Vector3d::Zero());

	// Export boundary volume and mass for visualization.
	learnSPH::write_particles_to_vtk(
		"../res/boundary_volume.vtk", boundary_particles, boundary_volumes);
	learnSPH::write_particles_to_vtk(
		"../res/boundary_mass.vtk", boundary_particles, boundary_masses);
	learnSPH::write_particles_to_vtk(
		"../res/fluid_density_boundary.vtk", fluid_particles, fluid_densities);

	const auto density_range = std::minmax_element(
		fluid_densities.begin(), fluid_densities.end());
	const auto volume_range = std::minmax_element(
		boundary_volumes.begin(), boundary_volumes.end());
	const auto mass_range = std::minmax_element(
		boundary_masses.begin(), boundary_masses.end());

	const auto &fluid_point_set = nsearch.point_set(fluid_id);
	std::size_t fluid_with_boundary_neighbors = 0;
	for (unsigned int i = 0; i < fluid_particles.size(); ++i)
	{
		if (fluid_point_set.n_neighbors(boundary_id, i) > 0)
		{
			++fluid_with_boundary_neighbors;
		}
	}

	std::cout << "Boundary particles: " << boundary_particles.size() << std::endl;
	std::cout << "Fluid particles: " << fluid_particles.size() << std::endl;
	std::cout << "Boundary volume range: ["
		<< *volume_range.first << ", " << *volume_range.second << "] m^3"
		<< std::endl;
	std::cout << "Boundary mass range: ["
		<< *mass_range.first << ", " << *mass_range.second << "] kg"
		<< std::endl;
	std::cout << "Fluid density range: ["
		<< *density_range.first << ", " << *density_range.second << "] kg/m^3"
		<< std::endl;
	std::cout << "Fluid particles with boundary neighbors: "
		<< fluid_with_boundary_neighbors << std::endl;
	std::cout << "Prepared " << fluid_velocities.size()
		<< " velocity, pressure, and acceleration buffers for Assignment 4."
		<< std::endl;
	std::cout << "Wrote boundary_volume.vtk, boundary_mass.vtk, and "
		<< "fluid_density_boundary.vtk to `<build_folder>/res`." << std::endl;

	(void)fluid_pressures;
	(void)fluid_accelerations;

	return 0;
}
