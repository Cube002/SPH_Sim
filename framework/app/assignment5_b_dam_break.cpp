#include <Eigen/Dense>

#include "assignment5_scenes.h"
#include "simulation_runner.h"

int main()
{
	using namespace learnsph_app;

	// Assignment 5 b: dam break in a rectangular tank.
	const SimulationSizes sizes = SimulationSizes::from_radius(0.005);

	SimulationConfig config;
	config.name = "Assignment 5 b: dam break";
	config.output_prefix = "assignment5_b_dam_break";
	config.stiffness = 1000.0;
	config.fluid_viscosity = 0.0025;
	config.boundary_viscosity = 0.0;
	config.gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
	config.max_time_step = 0.00025;
	config.simulation_time = 0.5;
	config.frame_interval = 0.01;
	config.fluid_mass_scale = 1.0;

	run_wcsph_simulation(sizes, dam_break_scene(sizes), config);
	return 0;
}
