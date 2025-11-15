MolSim
===
_Group C_

The Molecular Dynamics teaching code.

#### Build instructions
To build with GCC, use:
```
cmake --workflow release
```
The `MolSim` binary will be located under `build/release/MolSim`.

To build with Clang, use
```
cmake --workflow release-clang
```
The `MolSim` binary will be located under `build/release-gcc/MolSim`.

The classic CMake switches of the type `-DOPTION=ON` can be used to further customize the project.
Alternatively, one can edit the configuration presets within `CMakePresets.json` to set cache variables
for every run from then on.

#### Running the simulation
```
Usage: MolSim [--help] [--version] [--delta-t VAR] [--end-time VAR] [--output VAR] [--log-level VAR] files...

Positional arguments:
  files           input file(s) to read initial condition of particles from [nargs: 1 or more]

Optional arguments:
  -h, --help      shows help message and exits
  -v, --version   prints version information and exits
  -d, --delta-t   the tick length used for the simulation [nargs=0..1] [default: 0.014]
  -e, --end-time  point in time to simulate until [nargs=0..1] [default: 1000]
  -o, --output    output directory for resulting files [nargs=0..1] [default: "."]
  --log-level     verbosity of log output [nargs=0..1] [default: "info"]
```

#### Running tests
To run the tests after building with GCC, use:
```
cmake --workflow test
```
`test-clang` is also available.

#### Building documentation
To build the project documentation via Doxygen, use:
```
cmake --workflow documentation
```
The documentation will be avaliable under `build/release/documentation`.

Use a browser to open `documentation/html/index.html` to view the documentation.

> [!IMPORTANT]
> In order to properly parse the LaTeX within the documentation, a LaTeX engine (that provides the binary `latex`) and Ghostscript (that provides the binary `gs`) is required.

#### Nix and direnv
This project makes use of the nix package manager. If you have flakes enabled, you can use `nix develop` to get an identical development environment to the one that's used by the developers.

Or, running `direnv allow` will automatically put you into the environment whenever you change directory into
here, if you have `direnv` support enabled in your shell.

<details>
<summary>Building the documentation manually via `make`</summary>
After navigating to the binary folder created using the build instructions (`build/<PRESET_USED>`), run:
```
make doc_doxygen
```
The documentation will be avaliable under `documentation`.
</details>

<details>
<summary>Instructions for older versions of CMake than 3.31</summary>

- Below CMake 3.31:

Use `cmake --workflow --preset <WORKFLOW_NAME>` to run workflows, see [docs](https://cmake.org/cmake/help/latest/manual/cmake.1.html#cmdoption-cmake-workflow-preset).

- Below CMake 3.19:

Use `cmake -S . -B build -DENABLE_VTK_OUTPUT=ON && cmake --build build`, see [docs](https://cmake.org/cmake/help/latest/manual/cmake-presets.7.html#introduction).
</details>
