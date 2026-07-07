# FEM Solver

A 2D truss Finite Element Method (FEM) solver written in C++.
Reads a structure from a JSON file, assembles the global stiffness matrix,
solves for node displacements, recovers element stresses, and visualizes
the original and deformed structure using SFML.

## What it does

- Reads any 2D truss structure from `structure.json`
- Builds each bar element's 4×4 stiffness matrix using coordinate transforms
- Assembles the global stiffness matrix via the standard FEM assembly loop
- Applies boundary conditions (pinned supports, rollers) using the penalty method
- Solves the linear system K·u = F for node displacements (Eigen LDLT)
- Recovers and prints axial stress in each element
- Renders the original structure (gray) and deformed shape (red) in an SFML window

## Dependencies

- [Eigen](https://eigen.tuxfamily.org/) — linear algebra and sparse solve
- [nlohmann/json](https://github.com/nlohmann/json) — JSON input parsing
- [SFML 2.6](https://www.sfml-dev.org/) — 2D rendering and window management
- CMake 3.16+
- g++ with C++17 support

## Build

```bash
cmake -S . -B build
cmake --build build
```

## Run

```bash
./build/fem_solver
```

## Defining a structure (structure.json)

All geometry, connectivity, material properties, and loads are defined
in `structure.json`. No recompilation needed to change the structure.

```json
{
    "nodes": [
        {"x": 0.0, "y": 0.0},
        {"x": 2.0, "y": 0.0},
        {"x": 1.0, "y": 1.0}
    ],
    "elements": [
        {"node1": 0, "node2": 1, "E": 200000, "A": 0.01},
        {"node1": 1, "node2": 2, "E": 200000, "A": 0.01},
        {"node1": 0, "node2": 2, "E": 200000, "A": 0.01}
    ],
    "force": -1000.0,
    "forceDOF": 5
}
```

### Node fields

| Field | Description | Units |
|---|---|---|
| `x` | Horizontal position | m |
| `y` | Vertical position | m |

### Element fields

| Field | Description | Units |
|---|---|---|
| `node1` | Index of first node | — |
| `node2` | Index of second node | — |
| `E` | Young's modulus | Pa |
| `A` | Cross-sectional area | m² |

### Load fields

| Field | Description |
|---|---|
| `force` | Force magnitude (negative = downward) |
| `forceDOF` | Global DOF index where force is applied |

### Boundary conditions

Boundary conditions are currently applied by pinning specific DOFs
via the penalty method in `main.cpp`. Node 0 is fully pinned (DOFs 0
and 1), node 1 is a roller (DOF 3 fixed).

## Physics

Each bar element has axial stiffness `k = EA/L`. The local 2×2 stiffness
matrix is expanded to a 4×4 global-frame matrix using a coordinate
transform based on the bar's angle: `c = dx/L`, `s = dy/L`. Assembly
adds each element's 4×4 contribution into the global matrix at the DOF
slots corresponding to its two nodes. Stress is recovered after solving:
`σ = E · elongation / L`, where elongation projects the nodal
displacements onto the bar's axis.

## Roadmap

- [x] 1D axial bar solver
- [x] 2D truss solver with coordinate transforms
- [x] JSON-driven structure input
- [x] SFML visualization of original and deformed shape
- [ ] Dear ImGui interactive control panel
- [ ] Boundary conditions defined in JSON (not hardcoded)
- [ ] User-defined load direction and location via JSON