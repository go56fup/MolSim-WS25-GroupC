MolSim
===
_Group C_

The Molecular Dynamics teaching code.

#### Build instructions
To build with Clang, use:
```
cmake --workflow release
```
The `MolSim` binary will be located under `build/release/MolSim`.

To build with GCC, use
```
cmake --workflow release-gcc
```
The `MolSim` binary will be located under `build/release-gcc/MolSim`.

The classic CMake switches of the type `-DOPTION=ON` can be used to further customize the project.
Alternatively, one can edit the configuration presets within `CMakePresets.json` to set cache variables
for every run from then on.

#### Running tests
To run the tests after building with Clang, use:
```
cmake --workflow test
```
`test-gcc` is also available.

#### Building documentation
To build the project documentation via Doxygen, use:
```
cmake --workflow documentation
```
The documentation will be avaliable under `build/release/documentation`.

Use a browser to open `documentation/html/index.html` to view the documentation.

> [!IMPORTANT]
> In order to properly parse the LaTeX within the documentation, a LaTeX engine (that provides the binary `latex`) and Ghostscript (that provides the binary `gs`) is required.

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
