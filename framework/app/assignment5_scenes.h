#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "../learnSPH/sampling.h"
#include "simulation_runner.h"

namespace learnsph_app
{
	inline void append_fluid_box(
		SimulationScene &scene,
		const Eigen::Vector3d &bottom,
		const Eigen::Vector3d &top,
		double sampling_distance,
		const Eigen::Vector3d &velocity = Eigen::Vector3d::Zero())
	{
		std::vector<Eigen::Vector3d> particles;
		learnSPH::sampling::fluid_box(
			particles, bottom, top, sampling_distance);

		scene.fluid_particles.insert(
			scene.fluid_particles.end(), particles.begin(), particles.end());
		scene.fluid_velocities.insert(
			scene.fluid_velocities.end(), particles.size(), velocity);
	}

	inline SimulationScene single_fluid_box_scene(
		const SimulationSizes &sizes,
		const Eigen::Vector3d &bottom = Eigen::Vector3d::Zero(),
		const Eigen::Vector3d &top = Eigen::Vector3d::Ones())
	{
		SimulationScene scene;
		append_fluid_box(
			scene, bottom, top, sizes.fluid_sampling_distance);
		return scene;
	}

	inline SimulationScene two_fluid_blocks_scene(
		const SimulationSizes &sizes,
		double relative_velocity = 2.0)
	{
		// Two 1 m^3 blocks separated by 10 cm. The left block moves to +x and
		// the right block moves to -x, giving the requested relative speed.
		const double half_gap = 0.05;
		const double half_velocity = 0.5 * relative_velocity;

		SimulationScene scene;
		append_fluid_box(
			scene,
			Eigen::Vector3d(-1.0 - half_gap, 0.0, 0.0),
			Eigen::Vector3d(-half_gap, 1.0, 1.0),
			sizes.fluid_sampling_distance,
			Eigen::Vector3d(half_velocity, 0.0, 0.0));
		append_fluid_box(
			scene,
			Eigen::Vector3d(half_gap, 0.0, 0.0),
			Eigen::Vector3d(1.0 + half_gap, 1.0, 1.0),
			sizes.fluid_sampling_distance,
			Eigen::Vector3d(-half_velocity, 0.0, 0.0));
		return scene;
	}

	inline std::vector<Eigen::Vector3d> grid_plane_particles(
		double x_min, double x_max,
		double y_min, double y_max,
		double z,
		double spacing)
	{
		const int nx = std::max(
			1, static_cast<int>(std::round((x_max - x_min) / spacing)) + 1);
		const int ny = std::max(
			1, static_cast<int>(std::round((y_max - y_min) / spacing)) + 1);

		std::vector<Eigen::Vector3d> particles;
		particles.reserve(static_cast<size_t>(nx) * ny);
		for (int i = 0; i < nx; ++i)
		{
			const double x = (nx == 1)
				? x_min
				: x_min + i * (x_max - x_min) / (nx - 1);
			for (int j = 0; j < ny; ++j)
			{
				const double y = (ny == 1)
					? y_min
					: y_min + j * (y_max - y_min) / (ny - 1);
				particles.push_back(Eigen::Vector3d(x, y, z));
			}
		}
		return particles;
	}

	inline SimulationScene floor_scene(const SimulationSizes &sizes)
	{
		// Keep a small gap between the fluid and floor at t = 0. The floor is
		// sampled as a single layer of boundary particles.
		SimulationScene scene;
		append_fluid_box(
			scene,
			Eigen::Vector3d(0.0, 0.0, 0.05),
			Eigen::Vector3d(1.0, 1.0, 1.05),
			sizes.fluid_sampling_distance);
		scene.boundary_particles = grid_plane_particles(
			-1.0, 2.0, -1.0, 2.0, 0.0,
			sizes.boundary_sampling_distance);
		return scene;
	}

	inline std::vector<Eigen::Vector3d> rectangular_container_particles(
		const Eigen::Vector3d &box_min,
		const Eigen::Vector3d &box_max,
		double spacing)
	{
		const int nx = std::max(
			1, static_cast<int>(std::round((box_max[0] - box_min[0]) / spacing)) + 1);
		const int ny = std::max(
			1, static_cast<int>(std::round((box_max[1] - box_min[1]) / spacing)) + 1);
		const int nz = std::max(
			1, static_cast<int>(std::round((box_max[2] - box_min[2]) / spacing)) + 1);

		std::vector<Eigen::Vector3d> particles;
		for (int i = 0; i < nx; ++i)
		{
			const double x = (nx == 1)
				? box_min[0]
				: box_min[0] + i * (box_max[0] - box_min[0]) / (nx - 1);
			for (int j = 0; j < ny; ++j)
			{
				const double y = (ny == 1)
					? box_min[1]
					: box_min[1] + j * (box_max[1] - box_min[1]) / (ny - 1);
				for (int k = 0; k < nz; ++k)
				{
					const double z = (nz == 1)
						? box_min[2]
						: box_min[2] + k * (box_max[2] - box_min[2]) / (nz - 1);

					const bool on_floor = (k == 0);
					const bool on_left_wall = (i == 0);
					const bool on_right_wall = (i == nx - 1);
					const bool on_front_wall = (j == 0);
					const bool on_back_wall = (j == ny - 1);
					if (on_floor || on_left_wall || on_right_wall
						|| on_front_wall || on_back_wall)
					{
						particles.push_back(Eigen::Vector3d(x, y, z));
					}
				}
			}
		}
		return particles;
	}

	inline SimulationScene dam_break_scene(const SimulationSizes &sizes)
	{
		// Classic dam break, z up: x = length, y = depth, z = height.
		// Fluid: 0.25 x 0.15 x 0.5 m. Tank: 1.0 x 0.15 x 0.8 m.
		// The tank is enlarged by one particle diameter so boundary particles do
		// not overlap the initial fluid block.
		const double fluid_length = 0.25;
		const double fluid_depth = 0.15;
		const double fluid_height = 0.5;
		const double tank_length = 1.0;
		const double tank_depth = 0.15;
		const double tank_height = 0.8;
		const double enlargement = 2.0 * sizes.particle_radius;

		SimulationScene scene;
		append_fluid_box(
			scene,
			Eigen::Vector3d(0.0, 0.0, 0.0),
			Eigen::Vector3d(fluid_length, fluid_depth, fluid_height),
			sizes.fluid_sampling_distance);
		scene.boundary_particles = rectangular_container_particles(
			Eigen::Vector3d(-enlargement, -enlargement, -enlargement),
			Eigen::Vector3d(
				tank_length + enlargement,
				tank_depth + enlargement,
				tank_height + enlargement),
			sizes.boundary_sampling_distance);
		return scene;
	}
}
