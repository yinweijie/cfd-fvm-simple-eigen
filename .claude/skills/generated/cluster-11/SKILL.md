---
name: cluster-11
description: "Skill for the Cluster_11 area of cfd-fvm-simple-eigen. 6 symbols across 1 files."
---

# Cluster_11

6 symbols | 1 files | Cohesion: 100%

## When to Use

- Working with code in `src/`
- Understanding how pressure_index, u_index, v_index work
- Modifying cluster_11-related functionality

## Key Files

| File | Symbols |
|------|---------|
| `src/mesh.cpp` | pressure_index, u_index, v_index, cell_center_x, cell_center_y (+1) |

## Entry Points

Start here when exploring this area:

- **`pressure_index`** (Method) — `src/mesh.cpp:35`
- **`u_index`** (Method) — `src/mesh.cpp:42`
- **`v_index`** (Method) — `src/mesh.cpp:47`
- **`cell_center_x`** (Method) — `src/mesh.cpp:57`
- **`cell_center_y`** (Method) — `src/mesh.cpp:63`

## Key Symbols

| Symbol | Type | File | Line |
|--------|------|------|------|
| `pressure_index` | Method | `src/mesh.cpp` | 35 |
| `u_index` | Method | `src/mesh.cpp` | 42 |
| `v_index` | Method | `src/mesh.cpp` | 47 |
| `cell_center_x` | Method | `src/mesh.cpp` | 57 |
| `cell_center_y` | Method | `src/mesh.cpp` | 63 |
| `validate_pressure_ij` | Method | `src/mesh.cpp` | 69 |

## How to Explore

1. `gitnexus_context({name: "pressure_index"})` — see callers and callees
2. `gitnexus_query({query: "cluster_11"})` — find related execution flows
3. Read key files listed above for implementation details
