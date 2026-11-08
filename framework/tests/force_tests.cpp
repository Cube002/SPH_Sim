#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

#include "../learnSPH/io.h"
#include "../learnSPH/kernel.h"
#include "../learnSPH/sampling.h"
#include "../learnSPH/sph_math.h"

namespace
{
	struct CourseParticleSizes
	{
		double particle_radius;
		double particle_diameter;
		double fluid_sampling_distance;
		double boundary_sampling_distance;
		double smoothing_length;
		double compact_support;

		explicit CourseParticleSizes(double radius)
			: particle_radius(radius)
			, particle_diameter(2.0 * particle_radius)
			, fluid_sampling_distance(particle_diameter)
			, boundary_sampling_distance(0.8 * particle_diameter)
			, smoothing_length(1.2 * particle_diameter)
			, compact_support(2.0 * smoothing_length)
		{
		}
	};

	void check_vector(
		const Eigen::Vector3d &actual,
		const Eigen::Vector3d &expected)
	{
		for (int component = 0; component < 3; ++component)
		{
			CHECK(actual[component] == Catch::Approx(expected[component])
				.epsilon(1e-12).margin(1e-12));
		}
	}
}

TEST_CASE("compute_pressure applies and clamps the equation of state", "[forces][pressure]")
{
	const std::vector<double> densities = { 999.0, 1000.0, 1001.0, 1010.0 };
	std::vector<double> pressures;

	learnSPH::sph_math::compute_pressure(densities, pressures, 1000.0, 1000.0);

	REQUIRE(pressures.size() == densities.size());
	CHECK(pressures[0] == Catch::Approx(0.0));
	CHECK(pressures[1] == Catch::Approx(0.0));
	CHECK(pressures[2] == Catch::Approx(1000.0));
	CHECK(pressures[3] == Catch::Approx(10000.0));

	learnSPH::sph_math::compute_pressure(densities, pressures, 1000.0, 0.0);
	for (double pressure : pressures)
	{
		CHECK(pressure == Catch::Approx(0.0));
	}
}

TEST_CASE("pressure acceleration includes both fluid particles and adds to the buffer", "[forces][pressure]")
{
	const CourseParticleSizes sizes(0.0025);
	const double dx = sizes.fluid_sampling_distance;
	const double h = sizes.smoothing_length;
	const double particle_mass = 1000.0 * dx * dx * dx;
	const double density = 1010.0;
	const double pressure = 10000.0;

	std::vector<Eigen::Vector3d> particles = {
		Eigen::Vector3d::Zero(),
		Eigen::Vector3d(dx, 0.0, 0.0)
	};
	const std::vector<double> densities = { density, density };
	std::vector<double> pressures;
	learnSPH::sph_math::compute_pressure(densities, pressures, 1000.0, 1000.0);

	CompactNSearch::NeighborhoodSearch nsearch(sizes.compact_support);
	const unsigned int fluid_id = nsearch.add_point_set(
		particles[0].data(), particles.size(), false);
	nsearch.find_neighbors();

	std::vector<Eigen::Vector3d> accelerations = {
		Eigen::Vector3d(1.0, 2.0, 3.0),
		Eigen::Vector3d(-1.0, -2.0, -3.0)
	};

	learnSPH::sph_math::compute_pressure_acceleration(
		particles, densities, pressures, accelerations,
		nsearch, fluid_id, particle_mass, h);

	const Eigen::Vector3d grad_w = learnSPH::kernel::cubic_grad_spline(
		particles[0], particles[1], h);
	const Eigen::Vector3d contribution = -particle_mass
		* (pressure / (density * density) + pressure / (density * density))
		* grad_w;

	check_vector(accelerations[0], Eigen::Vector3d(1.0, 2.0, 3.0) + contribution);
	check_vector(accelerations[1], Eigen::Vector3d(-1.0, -2.0, -3.0) - contribution);
}

TEST_CASE("pressure acceleration includes boundary volume contribution", "[forces][pressure][boundary]")
{
	const CourseParticleSizes sizes(0.0025);
	const double dx = sizes.fluid_sampling_distance;
	const double h = sizes.smoothing_length;
	const double particle_mass = 1000.0 * dx * dx * dx;
	const double rest_density = 1000.0;
	const double density = 1010.0;
	const double pressure = 10000.0;
	const double boundary_volume = 0.5 * dx * dx * dx;

	std::vector<Eigen::Vector3d> fluid_particles = { Eigen::Vector3d::Zero() };
	std::vector<Eigen::Vector3d> boundary_particles = {
		Eigen::Vector3d(dx, 0.0, 0.0)
	};
	const std::vector<double> densities = { density };
	std::vector<double> pressures;
	learnSPH::sph_math::compute_pressure(
		densities, pressures, rest_density, 1000.0);
	const std::vector<double> boundary_volumes = { boundary_volume };

	CompactNSearch::NeighborhoodSearch nsearch(sizes.compact_support);
	const unsigned int fluid_id = nsearch.add_point_set(
		fluid_particles[0].data(), fluid_particles.size(), false);
	const unsigned int boundary_id = nsearch.add_point_set(
		boundary_particles[0].data(), boundary_particles.size(), false);
	nsearch.find_neighbors();

	std::vector<Eigen::Vector3d> accelerations(
		fluid_particles.size(), Eigen::Vector3d::Zero());
	learnSPH::sph_math::compute_pressure_acceleration(
		fluid_particles, boundary_particles, densities, pressures,
		boundary_volumes, accelerations, nsearch, fluid_id, boundary_id,
		particle_mass, rest_density, h);

	const Eigen::Vector3d grad_w = learnSPH::kernel::cubic_grad_spline(
		fluid_particles[0], boundary_particles[0], h);
	const Eigen::Vector3d expected = -rest_density * boundary_volume
		* (pressure / (density * density)) * grad_w;

	check_vector(accelerations[0], expected);
	CHECK(accelerations[0].norm() > 0.0);
}

TEST_CASE("fluid friction acceleration follows Eq. 13", "[forces][friction]")
{
	const CourseParticleSizes sizes(0.0025);
	const double dx = sizes.fluid_sampling_distance;
	const double h = sizes.smoothing_length;
	const double particle_mass = 1000.0 * dx * dx * dx;
	const double density = 1000.0;
	const double viscosity = 0.1;

	std::vector<Eigen::Vector3d> particles = {
		Eigen::Vector3d::Zero(),
		Eigen::Vector3d(dx, 0.0, 0.0)
	};
	std::vector<Eigen::Vector3d> velocities = {
		Eigen::Vector3d(1.0, 0.0, 0.0),
		Eigen::Vector3d(0.0, 0.0, 0.0)
	};
	const std::vector<double> densities = { density, density };

	CompactNSearch::NeighborhoodSearch nsearch(sizes.compact_support);
	const unsigned int fluid_id = nsearch.add_point_set(
		particles[0].data(), particles.size(), false);
	nsearch.find_neighbors();

	std::vector<Eigen::Vector3d> accelerations(
		particles.size(), Eigen::Vector3d::Zero());
	learnSPH::sph_math::compute_friction_acceleration(
		particles, velocities, densities, accelerations,
		nsearch, fluid_id, particle_mass, viscosity, h);

	const Eigen::Vector3d displacement = particles[0] - particles[1];
	const Eigen::Vector3d grad_w = learnSPH::kernel::cubic_grad_spline(
		particles[0], particles[1], h);
	const double denominator = displacement.squaredNorm() + 0.01 * h * h;
	const double directional_term = displacement.dot(grad_w) / denominator;
	const Eigen::Vector3d expected = 2.0 * viscosity
		* (particle_mass / density) * directional_term
		* (velocities[0] - velocities[1]);

	check_vector(accelerations[0], expected);
	check_vector(accelerations[1], -expected);
}

TEST_CASE("boundary friction supports stationary and moving boundaries", "[forces][friction][boundary]")
{
	const CourseParticleSizes sizes(0.0025);
	const double dx = sizes.fluid_sampling_distance;
	const double h = sizes.smoothing_length;
	const double particle_mass = 1000.0 * dx * dx * dx;
	const double density = 1000.0;
	const double boundary_volume = 0.5 * dx * dx * dx;
	const double boundary_viscosity = 1.0;

	std::vector<Eigen::Vector3d> fluid_particles = { Eigen::Vector3d::Zero() };
	std::vector<Eigen::Vector3d> boundary_particles = {
		Eigen::Vector3d(dx, 0.0, 0.0)
	};
	std::vector<Eigen::Vector3d> fluid_velocities = {
		Eigen::Vector3d(1.0, 0.0, 0.0)
	};
	const std::vector<double> densities = { density };
	const std::vector<double> boundary_volumes = { boundary_volume };

	CompactNSearch::NeighborhoodSearch nsearch(sizes.compact_support);
	const unsigned int fluid_id = nsearch.add_point_set(
		fluid_particles[0].data(), fluid_particles.size(), false);
	const unsigned int boundary_id = nsearch.add_point_set(
		boundary_particles[0].data(), boundary_particles.size(), false);
	nsearch.find_neighbors();

	const Eigen::Vector3d displacement =
		fluid_particles[0] - boundary_particles[0];
	const Eigen::Vector3d grad_w = learnSPH::kernel::cubic_grad_spline(
		fluid_particles[0], boundary_particles[0], h);
	const double denominator = displacement.squaredNorm() + 0.01 * h * h;
	const double directional_term = displacement.dot(grad_w) / denominator;

	std::vector<Eigen::Vector3d> accelerations(
		fluid_particles.size(), Eigen::Vector3d::Zero());
	const std::vector<Eigen::Vector3d> stationary_boundary_velocities;
	learnSPH::sph_math::compute_friction_acceleration(
		fluid_particles, boundary_particles, fluid_velocities,
		stationary_boundary_velocities, densities, boundary_volumes,
		accelerations, nsearch, fluid_id, boundary_id,
		particle_mass, 0.0, boundary_viscosity, h);

	const Eigen::Vector3d stationary_expected = 2.0 * boundary_viscosity
		* boundary_volume * directional_term * fluid_velocities[0];
	check_vector(accelerations[0], stationary_expected);

	const std::vector<Eigen::Vector3d> moving_boundary_velocities = {
		Eigen::Vector3d(0.0, 1.0, 0.0)
	};
	accelerations.assign(fluid_particles.size(), Eigen::Vector3d::Zero());
	learnSPH::sph_math::compute_friction_acceleration(
		fluid_particles, boundary_particles, fluid_velocities,
		moving_boundary_velocities, densities, boundary_volumes,
		accelerations, nsearch, fluid_id, boundary_id,
		particle_mass, 0.0, boundary_viscosity, h);

	const Eigen::Vector3d moving_expected = 2.0 * boundary_viscosity
		* boundary_volume * directional_term
		* (fluid_velocities[0] - moving_boundary_velocities[0]);
	check_vector(accelerations[0], moving_expected);
}

TEST_CASE("body acceleration adds gravity to all particles", "[forces][body]")
{
	std::vector<Eigen::Vector3d> accelerations = {
		Eigen::Vector3d(1.0, 2.0, 3.0),
		Eigen::Vector3d(-1.0, 0.0, 1.0)
	};
	const Eigen::Vector3d gravity(0.0, 0.0, -9.81);

	learnSPH::sph_math::compute_body_acceleration(accelerations, gravity);

	check_vector(accelerations[0], Eigen::Vector3d(1.0, 2.0, 3.0) + gravity);
	check_vector(accelerations[1], Eigen::Vector3d(-1.0, 0.0, 1.0) + gravity);
}

TEST_CASE("force functions work with a sampled OBJ boundary", "[forces][obj][integration]")
{
	const std::string box_path =
		std::string(LEARNING_SPH_SOURCE_DIR) + "/res/box.obj";
	const std::vector<learnSPH::TriMesh> meshes =
		learnSPH::read_tri_meshes_from_obj(box_path);
	REQUIRE(!meshes.empty());

	const learnSPH::TriMesh &box = meshes[0];
	const double boundary_sampling_distance = 0.4;
	const double h = 0.6;
	const double compact_support = 2.0 * h;
	const double boundary_density = 1000.0;
	const double rest_density = 1000.0;
	const double fluid_density = 1100.0;
	const double fluid_particle_mass = 1000.0 * 0.02 * 0.02 * 0.02;

	std::vector<Eigen::Vector3d> boundary_particles;
	learnSPH::sampling::triangle_mesh(
		boundary_particles, box.vertices, box.triangles,
		boundary_sampling_distance);
	REQUIRE(!boundary_particles.empty());

	std::vector<Eigen::Vector3d> fluid_particles = {
		Eigen::Vector3d(0.0, 0.0, 1.05)
	};

	CompactNSearch::NeighborhoodSearch nsearch(compact_support);
	const unsigned int boundary_id = nsearch.add_point_set(
		boundary_particles[0].data(), boundary_particles.size(), false);
	const unsigned int fluid_id = nsearch.add_point_set(
		fluid_particles[0].data(), fluid_particles.size(), false);
	nsearch.find_neighbors();

	std::vector<double> boundary_volumes;
	std::vector<double> boundary_masses;
	learnSPH::sph_math::compute_boundary_volumes(
		boundary_particles, boundary_volumes, nsearch, boundary_id, h);
	learnSPH::sph_math::compute_boundary_masses(
		boundary_volumes, boundary_masses, boundary_density);

	REQUIRE(boundary_volumes.size() == boundary_particles.size());
	REQUIRE(boundary_masses.size() == boundary_particles.size());
	for (double volume : boundary_volumes)
	{
		CHECK(volume > 0.0);
		CHECK(std::isfinite(volume));
	}

	const std::vector<double> fluid_densities = { fluid_density };
	std::vector<double> fluid_pressures;
	learnSPH::sph_math::compute_pressure(
		fluid_densities, fluid_pressures, rest_density, 1000.0);
	REQUIRE(fluid_pressures.size() == 1);
	CHECK(fluid_pressures[0] > 0.0);

	const auto &ps = nsearch.point_set(fluid_id);
	REQUIRE(ps.n_neighbors(boundary_id, 0) > 0);

	std::vector<Eigen::Vector3d> accelerations(
		fluid_particles.size(), Eigen::Vector3d::Zero());
	learnSPH::sph_math::compute_pressure_acceleration(
		fluid_particles, boundary_particles, fluid_densities, fluid_pressures,
		boundary_volumes, accelerations, nsearch, fluid_id, boundary_id,
		fluid_particle_mass, rest_density, h);
	CHECK(accelerations[0].allFinite());
	CHECK(accelerations[0].norm() > 0.0);

	std::vector<Eigen::Vector3d> fluid_velocities = {
		Eigen::Vector3d(0.5, 0.0, 0.0)
	};
	const std::vector<Eigen::Vector3d> boundary_velocities;
	accelerations.assign(fluid_particles.size(), Eigen::Vector3d::Zero());
	learnSPH::sph_math::compute_friction_acceleration(
		fluid_particles, boundary_particles, fluid_velocities,
		boundary_velocities, fluid_densities, boundary_volumes,
		accelerations, nsearch, fluid_id, boundary_id,
		fluid_particle_mass, 0.1, 0.5, h);
	CHECK(accelerations[0].allFinite());
	CHECK(accelerations[0].norm() > 0.0);
}
