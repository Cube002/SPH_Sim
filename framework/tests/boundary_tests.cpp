#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <vector>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

#include "../learnSPH/kernel.h"
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
}

TEST_CASE("boundary volume and mass follow the compensation equations", "[boundary]")
{
	const CourseParticleSizes sizes(0.0025);
	const double dx = sizes.boundary_sampling_distance;
	const double boundary_density = 1000.0;
	const double base_mass = boundary_density * dx * dx * dx;
	const double h = sizes.smoothing_length;

	std::vector<Eigen::Vector3d> boundary_particles = {
		Eigen::Vector3d(0.0, 0.0, 0.0),
		Eigen::Vector3d(dx, 0.0, 0.0)
	};

	CompactNSearch::NeighborhoodSearch nsearch(sizes.compact_support);
	const unsigned int boundary_id = nsearch.add_point_set(
		boundary_particles[0].data(), boundary_particles.size(), false);
	nsearch.find_neighbors();

	std::vector<double> volumes;
	std::vector<double> masses;
	learnSPH::sph_math::compute_boundary_volumes(
		boundary_particles, volumes, nsearch, boundary_id, h);
	learnSPH::sph_math::compute_boundary_masses(
		volumes, masses, boundary_density);

	const double self = learnSPH::kernel::cubic_spline(
		boundary_particles[0], boundary_particles[0], h);
	const double neighbor = learnSPH::kernel::cubic_spline(
		boundary_particles[0], boundary_particles[1], h);
	const double expected_volume = 1.0 / (self + neighbor);
	const double expected_mass = expected_volume * boundary_density;

	REQUIRE(volumes.size() == boundary_particles.size());
	REQUIRE(masses.size() == boundary_particles.size());
	CHECK(volumes[0] == Catch::Approx(expected_volume).epsilon(1e-12));
	CHECK(volumes[1] == Catch::Approx(expected_volume).epsilon(1e-12));
	CHECK(masses[0] == Catch::Approx(expected_mass).epsilon(1e-12));
	CHECK(masses[1] == Catch::Approx(expected_mass).epsilon(1e-12));

	std::vector<double> densities;
	std::vector<double> compensated_masses;
	learnSPH::sph_math::compute_density(
		boundary_particles, densities, nsearch, boundary_id, base_mass, h);
	learnSPH::sph_math::compute_compensated_mass(
		densities, compensated_masses, base_mass, boundary_density);

	CHECK(compensated_masses[0] == Catch::Approx(expected_mass).epsilon(1e-12));
	CHECK(compensated_masses[1] == Catch::Approx(expected_mass).epsilon(1e-12));
}

TEST_CASE("fluid density includes the boundary volume contribution", "[boundary][density]")
{
	const CourseParticleSizes sizes(0.0025);
	const double dx = sizes.fluid_sampling_distance;
	const double fluid_density = 1000.0;
	const double fluid_mass = fluid_density * dx * dx * dx;
	const double h = sizes.smoothing_length;
	const double boundary_volume = 0.25 * dx * dx * dx;

	std::vector<Eigen::Vector3d> fluid_particles(1, Eigen::Vector3d::Zero());
	std::vector<Eigen::Vector3d> boundary_particles = {
		Eigen::Vector3d(dx, 0.0, 0.0)
	};
	const std::vector<double> boundary_volumes = { boundary_volume };

	CompactNSearch::NeighborhoodSearch nsearch(sizes.compact_support);
	const unsigned int fluid_id = nsearch.add_point_set(
		fluid_particles[0].data(), fluid_particles.size(), false);
	const unsigned int boundary_id = nsearch.add_point_set(
		boundary_particles[0].data(), boundary_particles.size(), false);
	nsearch.find_neighbors();

	const auto &ps = nsearch.point_set(fluid_id);
	REQUIRE(ps.n_neighbors(boundary_id, 0) == 1);

	std::vector<double> densities;
	learnSPH::sph_math::compute_density(
		fluid_particles, boundary_particles, densities,
		nsearch, fluid_id, boundary_id, fluid_mass,
		boundary_volumes, fluid_density, h);

	const double self = learnSPH::kernel::cubic_spline(
		fluid_particles[0], fluid_particles[0], h);
	const double boundary_contribution = fluid_density * boundary_volume
		* learnSPH::kernel::cubic_spline(
			fluid_particles[0], boundary_particles[0], h);
	const double expected = fluid_mass * self + boundary_contribution;

	REQUIRE(densities.size() == fluid_particles.size());
	CHECK(densities[0] == Catch::Approx(expected).epsilon(1e-12));
	CHECK(densities[0] > fluid_mass * self);
}
