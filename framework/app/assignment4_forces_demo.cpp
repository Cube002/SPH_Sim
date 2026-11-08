#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

#include "../learnSPH/io.h"
#include "../learnSPH/sampling.h"
#include "../learnSPH/sph_math.h"

namespace
{
	double scalar_min(const std::vector<double> &values)
	{
		return *std::min_element(values.begin(), values.end());
	}

	double scalar_max(const std::vector<double> &values)
	{
		return *std::max_element(values.begin(), values.end());
	}
}

int main()
{
	// Assignment 4 force demo: one IcoSphere boundary and one fluid block.
	const double particle_radius = 0.0025;
	const double particle_diameter = 2.0 * particle_radius;
	const double fluid_sampling_distance = particle_diameter;
	const double boundary_sampling_distance = 0.8 * particle_diameter;
	const double smoothing_length = 1.2 * particle_diameter;
	const double compact_support = 2.0 * smoothing_length;

	const double fluid_density = 1000.0;
	const double boundary_density = fluid_density;
	const double rest_density = fluid_density;
	const double stiffness = 1000.0;
	const double fluid_viscosity = 0.1;
	const double boundary_viscosity = 1.0;
	const double fluid_particle_mass = fluid_density
		* fluid_sampling_distance
		* fluid_sampling_distance
		* fluid_sampling_distance;

	// Load and sample the provided boundary object.
	const std::vector<learnSPH::TriMesh> meshes =
		learnSPH::read_tri_meshes_from_obj("../obj/IcoSphere.obj");
	if (meshes.empty())
	{
		std::cerr << "Could not load ../obj/IcoSphere.obj" << std::endl;
		return 1;
	}
	const learnSPH::TriMesh &sphere = meshes[0];

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

	// Place the fluid block one particle diameter above the boundary.
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

	CompactNSearch::NeighborhoodSearch nsearch(compact_support);
	const unsigned int boundary_id = nsearch.add_point_set(
		boundary_particles[0].data(), boundary_particles.size(), true);
	const unsigned int fluid_id = nsearch.add_point_set(
		fluid_particles[0].data(), fluid_particles.size(), true);
	nsearch.find_neighbors(true);

	// Assignment 2/3: corrected boundary volumes and fluid-boundary density.
	std::vector<double> boundary_volumes;
	std::vector<double> boundary_masses;
	learnSPH::sph_math::compute_boundary_volumes(
		boundary_particles, boundary_volumes, nsearch, boundary_id,
		smoothing_length);
	learnSPH::sph_math::compute_boundary_masses(
		boundary_volumes, boundary_masses, boundary_density);

	std::vector<double> fluid_densities;
	learnSPH::sph_math::compute_density(
		fluid_particles, boundary_particles, fluid_densities,
		nsearch, fluid_id, boundary_id, fluid_particle_mass,
		boundary_volumes, fluid_density, smoothing_length);

	// Assignment 4: pressure, friction and body acceleration.
	std::vector<double> fluid_pressures;
	learnSPH::sph_math::compute_pressure(
		fluid_densities, fluid_pressures, rest_density, stiffness);

	std::vector<Eigen::Vector3d> fluid_velocities(
		fluid_particles.size(), Eigen::Vector3d(0.1, 0.0, 0.0));
	std::vector<Eigen::Vector3d> fluid_accelerations(
		fluid_particles.size(), Eigen::Vector3d::Zero());

	learnSPH::sph_math::compute_pressure_acceleration(
		fluid_particles, boundary_particles, fluid_densities, fluid_pressures,
		boundary_volumes, fluid_accelerations, nsearch,
		fluid_id, boundary_id, fluid_particle_mass, rest_density,
		smoothing_length);

	std::vector<Eigen::Vector3d> boundary_velocities;
	learnSPH::sph_math::compute_friction_acceleration(
		fluid_particles, boundary_particles, fluid_velocities,
		boundary_velocities, fluid_densities, boundary_volumes,
		fluid_accelerations, nsearch, fluid_id, boundary_id,
		fluid_particle_mass, fluid_viscosity, boundary_viscosity,
		smoothing_length);

	learnSPH::sph_math::compute_body_acceleration(
		fluid_accelerations, Eigen::Vector3d(0.0, 0.0, -9.81));

	std::vector<double> acceleration_magnitudes(fluid_particles.size());
	for (std::size_t i = 0; i < fluid_particles.size(); ++i)
	{
		acceleration_magnitudes[i] = fluid_accelerations[i].norm();
	}

	learnSPH::write_particles_to_vtk(
		"../res/assignment4_density.vtk", fluid_particles, fluid_densities);
	learnSPH::write_particles_to_vtk(
		"../res/assignment4_pressure.vtk", fluid_particles, fluid_pressures);
	learnSPH::write_particles_to_vtk(
		"../res/assignment4_acceleration.vtk", fluid_particles,
		acceleration_magnitudes, fluid_accelerations);
	learnSPH::write_particles_to_vtk(
		"../res/assignment4_velocity.vtk", fluid_particles,
		std::vector<double>(), fluid_velocities);

	std::cout << "Assignment 4 force demo" << std::endl;
	std::cout << "Boundary particles: " << boundary_particles.size() << std::endl;
	std::cout << "Fluid particles: " << fluid_particles.size() << std::endl;
	std::cout << "Density range: ["
		<< scalar_min(fluid_densities) << ", "
		<< scalar_max(fluid_densities) << "] kg/m^3" << std::endl;
	std::cout << "Pressure range: ["
		<< scalar_min(fluid_pressures) << ", "
		<< scalar_max(fluid_pressures) << "] Pa" << std::endl;
	std::cout << "Acceleration magnitude range: ["
		<< scalar_min(acceleration_magnitudes) << ", "
		<< scalar_max(acceleration_magnitudes) << "] m/s^2" << std::endl;
	std::cout << "Wrote Assignment 4 VTK files to `<build_folder>/res`."
		<< std::endl;

	return 0;
}
