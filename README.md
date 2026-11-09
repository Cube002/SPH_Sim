# learnSPH / SPH Simulation

This repository contains the incremental learnSPH coursework implementation:
particle sampling, neighborhood search, SPH kernels, boundary handling, density
estimation, forces, and the final WCSPH simulation.

The project is organized around exercise sheets. The documentation below
describes the current state and keeps a permanent recap of previous versions.

## Current Version

- **Version:** `v0.2.4`
- **Exercise sheet:** Exercise Sheet 2
- **Status:** Complete for Exercise Sheet 2, tasks 1-5, including 5b, 5c, boundary-mass export, corrected fluid-box sampling, and the enlarged 5c scene with a smaller obstacle

## Exercise Sheet Mapping

Exercise sheet task numbers are the official task numbers. File and executable
names are implementation names and do not necessarily match the task number.

### Exercise Sheet 1

Exercise Sheet 1 contains the introductory neighborhood-search and gravity
experiments:

| Executable | Purpose |
| --- | --- |
| `assignment1_search_speed` | Compares brute-force neighborhood search with `CompactNSearch`. |
| `assignment1_gravity_sim` | Basic particle simulation with gravity and neighborhood interaction. |

The related code is in:

- `framework/app/assignment1_search_speed.cpp`
- `framework/app/assignment1_gravity_sim.cpp`
- `framework/learnSPH/bruteforcesearch.h`
- `framework/learnSPH/bruteforcesearch.cpp`

### Exercise Sheet 2

Exercise Sheet 2 contains tasks 1-5:

| Task | Topic | Main implementation |
| --- | --- | --- |
| 1 | Cubic spline kernel and gradient | `framework/learnSPH/kernel.*`, `framework/tests/kernel_tests.cpp` |
| 2 | Sampling and boundary volume/mass correction | `framework/learnSPH/sampling.*`, boundary functions in `framework/learnSPH/sph_math.*`, `framework/tests/boundary_tests.cpp` |
| 3 | Fluid density with boundary contribution | `framework/app/assignment3_boundary_density.cpp`, `framework/tests/density_tests.cpp` |
| 4 | Pressure, friction, and body accelerations | `framework/app/assignment4_forces_demo.cpp`, `framework/tests/force_tests.cpp` |
| 5 | Putting it together: WCSPH simulation | `framework/app/assignment5_*`, `framework/app/simulation_runner.h`, `framework/app/assignment5_scenes.h` |

Exercise Sheet 2, task 5 contains these incremental test cases:

| Executable | Scene |
| --- | --- |
| `assignment5_test1_no_forces` | Baseline 1 m³ fluid block with all forces disabled. |
| `assignment5_a_i_pressure` | Add pressure stiffness `B = 1000`. |
| `assignment5_a_ii_viscosity` | Add fluid viscosity `nu_f = 0.1`. |
| `assignment5_a_iii_two_blocks` | Two 1 m³ blocks collide at 2 m/s relative velocity. |
| `assignment5_a_iv_gravity` | Remove the second block and enable gravity. |
| `assignment5_a_v_floor` | Add a boundary-particle floor. |
| `assignment5_a_vi_boundary_viscosity` | Add boundary viscosity `nu_b = 1.0`. |
| `assignment5_b_dam_break` | Official 5b dam-break scene: 0.15 x 0.25 x 0.5 m fluid in a 0.15 x 0.8 x 1.0 m tank. |
| `assignment5_c_obstacle_dam_break` | Official 5c showcase: 36,750 fluid particles, 20% smaller obstacle, 2 s simulation. |

The 5b scene is mapped to the z-up coordinate system as width x length x height.
The tank is enlarged by one particle diameter. In both 5b and 5c, particles that
leave the tank are removed from the simulation.

### Future Exercise Sheets

Additional exercise sheets are expected. When a new sheet is completed:

1. Add its executables and source files to this document.
2. Add a new version entry to the version history.
3. Keep all older history entries unchanged.

## Requirements

- CMake 3.10 or newer
- C++17 compiler
- OpenMP-capable compiler recommended (GCC/Clang)
- Linux shell commands are shown below

The project vendors Eigen and CompactNSearch under `framework/extern/`, and
Catch2 under `framework/tests/Catch2/`.

## Build

Configure and build from the repository root:

```bash
cmake -S framework -B build
cmake --build build -j"$(nproc)"
```

For simulation performance, use a Release build:

```bash
cmake -S framework -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release -j"$(nproc)"
```

The default build type in the existing `build/` directory may be `Debug`. The
dam-break simulation is significantly slower in Debug mode.

## Running the Executables

Run executables from the `app` directory inside the build folder. The programs
use relative paths such as `../res` and `../obj`.

```bash
cd build/app
./assignment1_search_speed
./assignment1_gravity_sim
./assignment3_boundary_density
./assignment4_forces_demo
./assignment5_test1_no_forces
./assignment5_a_i_pressure
./assignment5_a_ii_viscosity
./assignment5_a_iii_two_blocks
./assignment5_a_iv_gravity
./assignment5_a_v_floor
./assignment5_a_vi_boundary_viscosity
./assignment5_b_dam_break
./assignment5_c_obstacle_dam_break
```

The two large simulations accept optional positional overrides and flags:

```bash
./assignment5_b_dam_break [simulation_time] [threads] [frame_interval] [--export-boundary]
./assignment5_c_obstacle_dam_break [simulation_time] [threads] [frame_interval] [--export-boundary]
```

For example:

```bash
./assignment5_c_obstacle_dam_break 2.0 8 0.01 --export-boundary
```

`--export-boundary` writes an additional boundary VTK file for every frame,
using the corrected boundary mass as the scalar field. Use
`--no-export-boundary` to disable it again when a default enables it.

The CMake build copies `framework/res` and `framework/obj` into the build
directory, so `build/res` and `build/obj` are available to the executables.

### Executable Summary

| Executable | Description |
| --- | --- |
| `assignment1_search_speed` | Times brute-force search and `CompactNSearch` on sampled particle grids. |
| `assignment1_gravity_sim` | Demonstrates particle movement under gravity with neighborhood interaction. |
| `assignment3_boundary_density` | Samples an OBJ boundary, computes corrected boundary volumes/masses, and computes fluid density including boundary particles. |
| `assignment4_forces_demo` | Computes density, pressure, pressure acceleration, friction, and gravity for a fluid block above an OBJ sphere. |
| `assignment5_test1_no_forces` | First Assignment 5 baseline: no pressure, no viscosity, no gravity, no boundary. |
| `assignment5_a_i_pressure` | Assignment 5a(i): pressure enabled, no gravity. |
| `assignment5_a_ii_viscosity` | Assignment 5a(ii): pressure and fluid viscosity enabled. |
| `assignment5_a_iii_two_blocks` | Assignment 5a(iii): two fluid blocks collide at 2 m/s relative speed. |
| `assignment5_a_iv_gravity` | Assignment 5a(iv): pressure, fluid viscosity, and gravity enabled. |
| `assignment5_a_v_floor` | Assignment 5a(v): boundary-particle floor added. |
| `assignment5_a_vi_boundary_viscosity` | Assignment 5a(vi): large boundary viscosity enabled. |
| `assignment5_b_dam_break` | Assignment 5b: official dam-break simulation with `B = 1000`, `nu_f = 0.0025`, `nu_b = 0`, and `dt_max = 0.25 ms`. |
| `assignment5_c_obstacle_dam_break` | Assignment 5c: 36,750-particle obstacle dam break with a 20% smaller obstacle, 2 s simulation time, 201 exported frames. |

## Output and Visualization

Simulation and demonstration output is written to the build folder's `res`
directory:

```text
build/res/
build-release/res/
```

The main output format is binary VTK:

- `assignment3_boundary_density` writes `boundary_volume.vtk`,
  `boundary_mass.vtk`, and `fluid_density_boundary.vtk`.
- `assignment4_forces_demo` writes `assignment4_density.vtk`,
  `assignment4_pressure.vtk`, `assignment4_acceleration.vtk`, and
  `assignment4_velocity.vtk`.
- Assignment 5 executables write numbered fluid frames such as
  `assignment5_b_dam_break_0.vtk`, `assignment5_b_dam_break_1.vtk`, etc.
- With boundary export enabled, each frame also writes
  `assignment5_b_dam_break_boundary_0.vtk`,
  `assignment5_b_dam_break_boundary_1.vtk`, etc. These files contain the
  boundary particle positions and their corrected mass as a scalar field named
  `mass`.

Each fluid frame contains particle positions, density as a scalar field, and
velocity as a vector field. Use ParaView or Blender to inspect the frames.
The Blender importer ignores `_boundary_` frames; use ParaView to inspect the
boundary-mass files. The VTK files are ignored by Git.

The validated 2 s runs produced:

| Scene | Fluid particles | Frames | Final max speed | Culled particles | Loop time |
| --- | ---: | ---: | ---: | ---: | ---: |
| `assignment5_b_dam_break` | 18,750 | 201 | 1.250 m/s | 480 | 262 s |
| `assignment5_c_obstacle_dam_break` | 36,750 | 201 | 2.162 m/s | 611 | 502 s |

Both runs used 8 OpenMP threads and `frame_interval = 0.01`. All exported
positions, densities, and velocities remained finite, with density maxima around
`1110-1120 kg/m^3`.

### Blender Visualization

`tools/blender_vtk_sequence.py` imports a binary VTK particle sequence into
Blender and colors particles by velocity with Geometry Nodes. Run it from the
repository root with Blender's background mode:

```bash
blender --python tools/blender_vtk_sequence.py -- \
  --pattern 'build-release/res/assignment5_c_obstacle_dam_break_*.vtk' \
  --radius 0.006 --speed-max 8.0 --fps 20
```

To render the full imported sequence to PNG frames:

```bash
blender --python tools/blender_vtk_sequence.py -- \
  --pattern 'build-release/res/assignment5_c_obstacle_dam_break_*.vtk' \
  --radius 0.006 --speed-max 8.0 --fps 20 --render \
  --output '//sph_animation/'
```

The script creates a 1920x1080 scene, a camera, lights, a tank floor, and an
animated point cloud. It uses a fixed-size mesh and a per-particle radius
attribute because particle counts change when escaped particles are culled. A
velocity-colored material maps `0 .. speed-max` from blue through cyan and
yellow to red. Rendering the 201 frames at 20 FPS gives an approximately 10 s
video; for a slower showcase, lower the FPS or duplicate frames.

Visualization challenges addressed by the script:

- Blender has no built-in importer for the binary, big-endian VTK files written
  by `learnSPH::write_particles_to_vtk`.
- Creating 20,000+ separate Blender objects is impractical, so the script uses
  one mesh plus Geometry Nodes instancing.
- Culling changes the particle count between frames, so the script keeps a
  fixed maximum-size mesh and sets the radius of missing particles to zero.
- Velocity is stored as a point attribute and converted to a color ramp in the
  shader, rather than assigning materials per object.
- Geometry Nodes and shader attribute APIs differ between Blender versions; the
  script targets Blender 4.x.

## Simulation Configuration

The Assignment 5a, 5b, and 5c executables based on
`run_wcsph_simulation()` use `learnsph_app::SimulationConfig` from
`framework/app/simulation_runner.h`. Change values in the corresponding
executable, then rebuild. `assignment5_test1_no_forces` keeps its original
standalone fixed-step setup.

| Field | Meaning |
| --- | --- |
| `name` | Display name printed at startup. |
| `output_prefix` | Prefix for exported VTK files. |
| `rest_density` | Rest density `rho_0` in kg/m³. |
| `stiffness` | WCSPH pressure stiffness `B`. |
| `fluid_viscosity` | Fluid viscosity `nu_f`. |
| `boundary_viscosity` | Boundary viscosity `nu_b`. |
| `gravity` | Gravity/body acceleration vector. |
| `max_time_step` | Maximum/default time step. |
| `simulation_time` | Total simulated time. |
| `frame_interval` | Time between exported VTK frames. |
| `fluid_mass_scale` | Multiplier applied to the sampled fluid particle mass. |
| `cfl_factor` | CFL safety factor used by `dt = min(dt_CFL, dt_default)`. |
| `use_cfl` | Enables or disables the CFL time-step limit. |
| `use_openmp` | Enables or disables OpenMP execution. |
| `num_threads` | OpenMP thread count. `0` means automatic, capped at 8. |
| `cull_outside_bounds` | Removes escaped or non-finite fluid particles before each step. |
| `bounds_min` | Minimum corner of the allowed fluid region when culling is enabled. |
| `bounds_max` | Maximum corner of the allowed fluid region when culling is enabled. |
| `max_wall_clock_seconds` | Graceful wall-clock limit for very long simulations; default 3600 s. |
| `export_boundary_particles` | Writes a separate boundary VTK frame with corrected boundary mass as scalar data. |

Example:

```cpp
SimulationConfig config;
config.stiffness = 1000.0;
config.fluid_viscosity = 0.1;
config.gravity = Eigen::Vector3d(0.0, 0.0, -9.81);
config.use_openmp = true;
config.num_threads = 4;
config.cull_outside_bounds = true;
config.bounds_min = Eigen::Vector3d(-0.01, -0.01, -0.01);
config.bounds_max = Eigen::Vector3d(1.21, 0.36, 0.71);
config.export_boundary_particles = true;
```

The `OMP_NUM_THREADS` environment variable is also respected when
`num_threads = 0`.

## Timer

Every executable based on `run_wcsph_simulation()` prints:

- setup time
- simulation-loop time
- total wall-clock time

Example:

```text
  OpenMP: enabled, threads = 8
  finished at t = 0.1 s after 110 steps
  setup time: 0.00861902 s
  simulation loop time: 0.516039 s
  total wall-clock time: 0.524658 s
```

## Project Structure

```text
SPH_Sim/
├── README.md
├── KERNEL_DOCUMENTATION.md
├── NOTES.md
├── NotesForPresentation.md
├── tools/
│   └── blender_vtk_sequence.py    # Binary VTK sequence importer for Blender
├── build/                         # Local build output, ignored by Git
└── framework/
    ├── CMakeLists.txt             # Top-level build, OpenMP, dependency setup
    ├── app/
    │   ├── CMakeLists.txt
    │   ├── assignment1_search_speed.cpp
    │   ├── assignment1_gravity_sim.cpp
    │   ├── assignment3_boundary_density.cpp
    │   ├── assignment4_forces_demo.cpp
    │   ├── assignment5_test1_no_forces.cpp
    │   ├── assignment5_a_i_pressure.cpp
    │   ├── assignment5_a_ii_viscosity.cpp
    │   ├── assignment5_a_iii_two_blocks.cpp
    │   ├── assignment5_a_iv_gravity.cpp
    │   ├── assignment5_a_v_floor.cpp
    │   ├── assignment5_a_vi_boundary_viscosity.cpp
    │   ├── assignment5_b_dam_break.cpp
    │   ├── assignment5_c_obstacle_dam_break.cpp
    │   ├── assignment5_scenes.h     # Reusable Assignment 5 scene builders
    │   └── simulation_runner.h      # WCSPH loop, CFL, culling, timer, export, OpenMP
    ├── learnSPH/
    │   ├── CMakeLists.txt
    │   ├── bruteforcesearch.h/.cpp
    │   ├── io.h/.cpp
    │   ├── kernel.h/.cpp
    │   ├── parallel.h               # OpenMP runtime configuration
    │   ├── sampling.h/.cpp
    │   └── sph_math.h/.cpp          # Density, boundary, pressure, friction, body forces
    ├── tests/
    │   ├── CMakeLists.txt
    │   ├── kernel_tests.cpp
    │   ├── sample_tests.cpp
    │   ├── boundary_tests.cpp
    │   ├── density_tests.cpp
    │   └── force_tests.cpp
    ├── obj/
    │   └── IcoSphere.obj
    ├── res/
    │   └── box.obj
    └── extern/
        ├── Eigen/
        ├── CompactNSearch/
        └── tests/Catch2/
```

Generated simulation frames are written to `build/res` and are not source
files.

## Core Modules

- `learnSPH::kernel`: cubic spline kernel and gradient.
- `learnSPH::sampling`: fluid-box and triangle-mesh sampling.
- `learnSPH::sph_math`: boundary correction, density, pressure, pressure
  acceleration, friction acceleration, and body acceleration.
- `learnSPH::parallel`: runtime OpenMP enable/disable and thread count.
- `learnsph_app::SimulationSizes`: particle radius and derived sampling lengths.
- `learnsph_app::SimulationConfig`: scene-independent simulation parameters.
- `learnsph_app::SimulationScene`: fluid/boundary particles and velocities.
- `learnsph_app::run_wcsph_simulation`: Algorithm 1 loop, time integration,
  VTK export, and timing.

## Algorithm Flow

The simulation runner follows the assignment outline:

1. Build the scene (fluid and optional boundary particles).
2. Compute corrected boundary volumes and masses once.
3. For each time step:
   1. Remove particles outside the configured bounds or with non-finite state.
   2. Update `CompactNSearch` neighbors.
   3. Compute fluid densities.
   4. Compute pressures.
   5. Compute pressure, friction, and body accelerations.
   6. Compute `dt = min(dt_CFL, dt_default)`.
   7. Perform semi-implicit Euler integration.
   8. Export a VTK frame when due.
   9. Stop gracefully if the wall-clock limit is reached.

## Tests

Build and run the test binary directly:

```bash
cmake --build build -j"$(nproc)"
./build/tests/learnSPH_tests
```

Current test result: **934 assertions in 19 test cases**.

The tests cover:

- cubic spline values, derivatives, and compact support
- box and triangle-mesh sampling
- boundary volume/mass correction
- density with self, fluid, and boundary contributions
- pressure equation of state and pressure/friction/body accelerations
- OBJ-based fluid-boundary force integration

CTest is not currently configured, so use the test binary directly.

## OpenMP Parallelization

OpenMP is enabled in CMake when available. The main SPH loops use
`#pragma omp parallel for`, including:

- density computation
- boundary volume/mass computation
- pressure computation
- pressure, friction, and body acceleration
- semi-implicit Euler integration
- CFL maximum-speed reduction

Small workloads below 128 particles run serially to avoid OpenMP overhead.

Runtime control is available through `SimulationConfig`:

```cpp
config.use_openmp = false;
config.num_threads = 1;
```

By default, automatic mode uses up to 8 threads. The built-in timer makes it
easy to compare serial and parallel runs.

## Known Limitations

- Moving boundaries are prepared at the API level through
  `boundary_velocities`, but scene updates and boundary re-sampling are not
  implemented yet.
- The floor scene uses a large single-layer boundary particle plane.
- The 5b and 5c executables run to `2.0 s` by default. Use the positional
  simulation-time override for shorter tests or longer showcase runs.
- `CompactNSearch::find_neighbors()` and VTK export are still serial, so
  parallel speedup is limited by those stages.
- Culling changes the fluid particle count, so later frames contain fewer
  particles than the initial frame.
- The Blender importer targets Blender 4.x and was not run in this environment;
  it requires a local Blender installation.
- CTest integration is not configured; tests are run as a single binary.

## Versioning Policy

The repository uses exercise-sheet-based versions:

- `v0.1.0`: Exercise Sheet 1
- `v0.2.0`: Exercise Sheet 2 baseline
- `v0.2.1`: Completion of Exercise Sheet 2 tasks 5b and 5c
- `v0.2.2`: Optional boundary-particle mass export
- `v0.2.3`: Corrected fluid-box sampling near boundaries
- `v0.2.4`: Larger 5c fluid reservoir and smaller obstacle
- `v0.3.0`: Exercise Sheet 3 (future)
- Patch versions are used for additional work or fixes within the same exercise
  sheet.

When creating a new version:

1. Add a row to the version history table below.
2. Summarize the new exercise sheet and its tasks.
3. Keep all older rows; never replace history.
4. Update `Current Version` at the top of this document.
5. Keep this README in the repository so future versions continue the recap.

### Version History

| Version | Exercise Sheet | Summary |
| --- | --- | --- |
| `v0.1.0` | Exercise Sheet 1 | Initial framework usage, brute-force neighborhood search, `CompactNSearch` comparison, and gravity particle demonstration. |
| `v0.2.0` | Exercise Sheet 2 | Tasks 1-5 baseline: cubic spline kernel and tests, sampling and boundary volume/mass correction, density with boundary contribution, pressure/friction/body forces, incremental Assignment 5 scenes, dam-break simulation, VTK velocity export, wall-clock timer, OpenMP parallelization, and project documentation. |
| `v0.2.1` | Exercise Sheet 2 | Completed 5b with the official 0.15 x 0.25 x 0.5 m fluid and 0.15 x 0.8 x 1.0 m tank configuration; added the 5c 23,660-particle obstacle dam-break showcase; added escaped-particle culling, non-finite particle removal, command-line simulation overrides, graceful 1 h wall-clock stopping, 2 s validation runs, and the Blender binary-VTK sequence importer with velocity coloring. |
| `v0.2.2` | Exercise Sheet 2 | Added optional boundary-particle VTK export with corrected boundary mass as scalar data, including the `--export-boundary` and `--no-export-boundary` command-line flags. |
| `v0.2.3` | Exercise Sheet 2 | Fixed fluid-box sampling so particles are consistently placed half a spacing inside every box face. This removes the initial over-density at the floor/walls and the resulting startup bounce. Updated 5b to 18,750 fluid particles and 5c to 21,875 fluid particles. |
| `v0.2.4` | Exercise Sheet 2 | Enlarged the 5c fluid reservoir to 0.35 x 0.35 x 0.30 m and reduced the obstacle to 80% of its original width and height. Release-mode 2 s run produced 36,750 initial fluid particles, 201 frames, final max speed 2.162 m/s, and 611 culled particles. |

### Adding the Next Version

Use this template when Exercise Sheet 3 is completed:

```markdown
| `v0.3.0` | Exercise Sheet 3 | Brief summary of tasks and new modules/executables. |
```

Then add a corresponding short section above the version table describing the
new sheet, its tasks, and the executables that were added.
