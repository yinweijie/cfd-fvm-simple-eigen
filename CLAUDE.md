<!-- gitnexus:start -->
# GitNexus — Code Intelligence

This project is indexed by GitNexus as **cfd-fvm-simple-eigen** (741 symbols, 1604 relationships, 63 execution flows). Use the GitNexus MCP tools to understand code, assess impact, and navigate safely.

> If any GitNexus tool warns the index is stale, run `npx gitnexus analyze` in terminal first.

## Always Do

- **MUST run impact analysis before editing any symbol.** Before modifying a function, class, or method, run `gitnexus_impact({target: "symbolName", direction: "upstream"})` and report the blast radius (direct callers, affected processes, risk level) to the user.
- **MUST run `gitnexus_detect_changes()` before committing** to verify your changes only affect expected symbols and execution flows.
- **MUST warn the user** if impact analysis returns HIGH or CRITICAL risk before proceeding with edits.
- When exploring unfamiliar code, use `gitnexus_query({query: "concept"})` to find execution flows instead of grepping. It returns process-grouped results ranked by relevance.
- When you need full context on a specific symbol — callers, callees, which execution flows it participates in — use `gitnexus_context({name: "symbolName"})`.

## Never Do

- NEVER edit a function, class, or method without first running `gitnexus_impact` on it.
- NEVER ignore HIGH or CRITICAL risk warnings from impact analysis.
- NEVER rename symbols with find-and-replace — use `gitnexus_rename` which understands the call graph.
- NEVER commit changes without running `gitnexus_detect_changes()` to check affected scope.

## Resources

| Resource | Use for |
|----------|---------|
| `gitnexus://repo/cfd-fvm-simple-eigen/context` | Codebase overview, check index freshness |
| `gitnexus://repo/cfd-fvm-simple-eigen/clusters` | All functional areas |
| `gitnexus://repo/cfd-fvm-simple-eigen/processes` | All execution flows |
| `gitnexus://repo/cfd-fvm-simple-eigen/process/{name}` | Step-by-step execution trace |

## CLI

| Task | Read this skill file |
|------|---------------------|
| Understand architecture / "How does X work?" | `.claude/skills/gitnexus/gitnexus-exploring/SKILL.md` |
| Blast radius / "What breaks if I change X?" | `.claude/skills/gitnexus/gitnexus-impact-analysis/SKILL.md` |
| Trace bugs / "Why is X failing?" | `.claude/skills/gitnexus/gitnexus-debugging/SKILL.md` |
| Rename / extract / split / refactor | `.claude/skills/gitnexus/gitnexus-refactoring/SKILL.md` |
| Tools, resources, schema reference | `.claude/skills/gitnexus/gitnexus-guide/SKILL.md` |
| Index, status, clean, wiki CLI commands | `.claude/skills/gitnexus/gitnexus-cli/SKILL.md` |
| Work in the Cfd area (51 symbols) | `.claude/skills/generated/cfd/SKILL.md` |
| Work in the Scripts area (21 symbols) | `.claude/skills/generated/scripts/SKILL.md` |
| Work in the Cluster_8 area (12 symbols) | `.claude/skills/generated/cluster-8/SKILL.md` |
| Work in the Cluster_12 area (12 symbols) | `.claude/skills/generated/cluster-12/SKILL.md` |
| Work in the Tests area (7 symbols) | `.claude/skills/generated/tests/SKILL.md` |
| Work in the Cluster_6 area (6 symbols) | `.claude/skills/generated/cluster-6/SKILL.md` |
| Work in the Cluster_9 area (6 symbols) | `.claude/skills/generated/cluster-9/SKILL.md` |
| Work in the Cluster_11 area (6 symbols) | `.claude/skills/generated/cluster-11/SKILL.md` |
| Work in the Cluster_10 area (3 symbols) | `.claude/skills/generated/cluster-10/SKILL.md` |

<!-- gitnexus:end -->
