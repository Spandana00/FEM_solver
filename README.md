# FEM Solver

A 1D axial bar Finite Element Method (FEM) solver written in C++.
Computes node displacements and element stresses for a chain of bars
in series under an axial load.

## What it does

- Reads structure definition from `structure.json`
- Assembles the global stiffness matrix from individual bar elements
- Applies a fixed boundary condition at node 0 (wall)
- Solves the linear system K·u = F for node displacements
- Recovers and prints axial stress in each element

## Dependencies

- [Eigen](https://eigen.tuxfamily.org/) — linear algebra
- [nlohmann/json](https://github.com/nlohmann/json) — JSON parsing
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

## Structure definition (structure.json)

```json
{
    "numElements": 5,
    "barLength": 1.0,
    "E": 200.0,
    "A": 100.0,
    "force": 1000.0
}
```

| Field | Description | Units |
|---|---|---|
| `numElements` | Number of bar elements | — |
| `barLength` | Length of each bar | m |
| `E` | Young's modulus | Pa |
| `A` | Cross-sectional area | m² |
| `force` | Applied axial force at free end | N |

## Physics

Each bar element has stiffness `k = EA/L`. The global stiffness matrix
is assembled by summing element contributions at shared nodes. Stress
is recovered after solving: `σ = E · (u₂ - u₁) / L`.
