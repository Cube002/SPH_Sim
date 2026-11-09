#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <vector>

#include <Eigen/Dense>
#include <CompactNSearch/CompactNSearch>

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
}

TEST_CASE("compute_density includes the self contribution", "[density]")
{
	const CourseParticleSizes sizes(0.0025);
	const double rho_f = 1000.0;
	const double particle_mass = rho_f
		* sizes.fluid_sampling_distance
		* sizes.fluid_sampling_distance
		* sizes.fluid_sampling_distance;

	std::vector<Eigen::Vector3d> particles(1, Eigen::Vector3d::Zero());
	CompactNSearch::NeighborhoodSearch nsearch(sizes.compact_support);
	const unsigned int fluid_id = nsearch.add_point_set(
		particles[0].data(), particles.size(), false);
	nsearch.find_neighbors();

	std::vector<double> densities;
	learnSPH::sph_math::compute_density(
		particles, densities, nsearch, fluid_id,
		particle_mass, sizes.smoothing_length);

	const double expected = particle_mass * learnSPH::kernel::cubic_spline(
		particles[0], particles[0], sizes.smoothing_length);

	REQUIRE(densities.size() == 1);
	CHECK(densities[0] == Catch::Approx(expected).epsilon(1e-12));
}

TEST_CASE("compute_density includes each neighbor exactly once", "[density]")
{
	const CourseParticleSizes sizes(0.0025);
	const double rho_f = 1000.0;
	const double dx = sizes.fluid_sampling_distance;
	const double particle_mass = rho_f * dx * dx * dx;

	std::vector<Eigen::Vector3d> particles = {
		Eigen::Vector3d(0.0, 0.0, 0.0),
		Eigen::Vector3d(dx, 0.0, 0.0)
	};

	CompactNSearch::NeighborhoodSearch nsearch(sizes.compact_support);
	const unsigned int fluid_id = nsearch.add_point_set(
		particles[0].data(), particles.size(), false);
	nsearch.find_neighbors();

	const auto &ps = nsearch.point_set(fluid_id);
	REQUIRE(ps.n_neighbors(fluid_id, 0) == 1);
	REQUIRE(ps.n_neighbors(fluid_id, 1) == 1);

	std::vector<double> densities;
	learnSPH::sph_math::compute_density(
		particles, densities, nsearch, fluid_id,
		particle_mass, sizes.smoothing_length);

	const double self = learnSPH::kernel::cubic_spline(
		particles[0], particles[0], sizes.smoothing_length);
	const double neighbor = learnSPH::kernel::cubic_spline(
		particles[0], particles[1], sizes.smoothing_length);
	const double expected = particle_mass * (self + neighbor);

	REQUIRE(densities.size() == particles.size());
	CHECK(densities[0] == Catch::Approx(expected).epsilon(1e-12));
	CHECK(densities[1] == Catch::Approx(expected).epsilon(1e-12));
}

TEST_CASE("uniform fluid box interior density matches the reference density", "[density][box]")
{
	const CourseParticleSizes sizes(0.0025);
	const double rho_f = 1000.0;
	const double dx = sizes.fluid_sampling_distance;
	const double particle_mass = rho_f * dx * dx * dx;
	const int n_per_axis = 10;

	std::vector<Eigen::Vector3d> particles;
	learnSPH::sampling::fluid_box(
		particles,
		Eigen::Vector3d::Zero(),
		Eigen::Vector3d::Constant(n_per_axis * dx),
		dx);

	REQUIRE(particles.size() == static_cast<std::size_t>(n_per_axis * n_per_axis * n_per_axis));
	CHECK(particles.front()[0] == Catch::Approx(0.5 * dx).epsilon(1e-12));
	CHECK(particles.front()[1] == Catch::Approx(0.5 * dx).epsilon(1e-12));
	CHECK(particles.front()[2] == Catch::Approx(0.5 * dx).epsilon(1e-12));
	CHECK(particles.back()[0] == Catch::Approx((n_per_axis - 0.5) * dx).epsilon(1e-12));
	CHECK(particles.back()[1] == Catch::Approx((n_per_axis - 0.5) * dx).epsilon(1e-12));
	CHECK(particles.back()[2] == Catch::Approx((n_per_axis - 0.5) * dx).epsilon(1e-12));

	CompactNSearch::NeighborhoodSearch nsearch(sizes.compact_support);
	const unsigned int fluid_id = nsearch.add_point_set(
		particles[0].data(), particles.size(), false);
	nsearch.find_neighbors();

	std::vector<double> densities;
	learnSPH::sph_math::compute_density(
		particles, densities, nsearch, fluid_id,
		particle_mass, sizes.smoothing_length);

	double interior_sum = 0.0;
	std::size_t interior_count = 0;
	const double lower = sizes.compact_support;
	const double upper = (n_per_axis - 1) * dx - sizes.compact_support;

	for (std::size_t i = 0; i < particles.size(); ++i)
	{
		const Eigen::Vector3d &p = particles[i];
		if ((p.array() >= lower).all() && (p.array() <= upper).all())
		{
			interior_sum += densities[i];
			++interior_count;
		}
	}

	REQUIRE(interior_count > 0);
	const double interior_mean = interior_sum / static_cast<double>(interior_count);
	CHECK(interior_mean == Catch::Approx(rho_f).epsilon(0.02));
}
