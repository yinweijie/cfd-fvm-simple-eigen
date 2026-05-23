---
name: tests
description: "Skill for the Tests area of cfd-fvm-simple-eigen. 7 symbols across 3 files."
---

# Tests

7 symbols | 3 files | Cohesion: 100%

## When to Use

- Working with code in `tests/`
- Understanding how main, solver, solver work
- Modifying tests-related functionality

## Key Files

| File | Symbols |
|------|---------|
| `tests/test_cli_smoke.cpp` | has_solver_marker, has_standard_outputs, main |
| `tests/test_projection_smoke.cpp` | expect_near, solver |
| `tests/test_simple_smoke.cpp` | expect_near, solver |

## Entry Points

Start here when exploring this area:

- **`main`** (Function) — `tests/test_cli_smoke.cpp:29`
- **`solver`** (Function) — `tests/test_projection_smoke.cpp:21`
- **`solver`** (Function) — `tests/test_simple_smoke.cpp:23`

## Key Symbols

| Symbol | Type | File | Line |
|--------|------|------|------|
| `main` | Function | `tests/test_cli_smoke.cpp` | 29 |
| `solver` | Function | `tests/test_projection_smoke.cpp` | 21 |
| `solver` | Function | `tests/test_simple_smoke.cpp` | 23 |
| `has_solver_marker` | Function | `tests/test_cli_smoke.cpp` | 8 |
| `has_standard_outputs` | Function | `tests/test_cli_smoke.cpp` | 19 |
| `expect_near` | Function | `tests/test_projection_smoke.cpp` | 7 |
| `expect_near` | Function | `tests/test_simple_smoke.cpp` | 7 |

## How to Explore

1. `gitnexus_context({name: "main"})` — see callers and callees
2. `gitnexus_query({query: "tests"})` — find related execution flows
3. Read key files listed above for implementation details
