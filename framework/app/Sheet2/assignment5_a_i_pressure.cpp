#include <Eigen/Dense>

#include "assignment5_scenes.h"
#include "simulation_runner.h"

int main()
{
	using namespace learnsph_app;

	// Assignment 5 a(i): single 1 m^3 block, pressure enabled, no gravity.
	const SimulationSizes sizes = SimulationSizes::from_radius(0.03);

	SimulationConfig config;
	config.name = "Assignment 5 a(i): pressure stiffness";
	config.output_prefix = "assignment5_a_i_pressure";
	config.stiffness = 1000.0;
	config.fluid_viscosity = 0.0;
	config.boundary_viscosity = 0.0;
	config.max_time_step = 0.001;
	config.simulation_time = 0.1;
	config.frame_interval = 0.01;
	config.fluid_mass_scale = 0.9;

	run_wcsph_simulation(
		sizes, single_fluid_box_scene(sizes), config);
	return 0;
}
