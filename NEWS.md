# AdapSamp News

## Version 1.2.0 (2026-08-25)

### Bug Fixes

- **rARS**: Removed redundant `solve()` call (computed twice per iteration).
  Added special handling for degenerate `tangent = 0` segments (uniform
  envelope) to prevent division-by-zero NaN. Added `p(x) <= 0` guard in
  `V(x)` to avoid `log(0)` / `log(negative)` crashes.

- **rMARS**: Replaced `ifelse()` control-flow misuse with proper `if/else`
  for `x_star` assignment (#P0-6). Added fallback for `deriv2 == 0`
  (inflexion-point boundary) to prevent silent region skipping (#P0-7).
  Moved `deriv1`/`deriv2` symbolic derivative construction outside the
  while-loop to avoid redundant re-parsing (#P1-5). Added `is.finite()`
  guards on integration limits. Fixed spelling "shouble" -> "should".

- **rCCARS**: Added duplicate-point removal after sorting cross-points
  to prevent zero-length `integrate()` calls (#P1-8). Added NaN guard
  in inverse-CDF formula with uniform-sampling fallback (#P0-9). Renamed
  `rubbish1` to `merged`. Fixed concave segment replication logic.

- **rASS**: Changed acceptance from `>` to `>=` to handle floating-point
  edge cases (#P0-11). Added `max_steps` parameter to bound stepping-out
  and shrinkage loops, preventing infinite loops on heavy-tailed
  distributions (#P0-10). Fixed return length: now returns exactly `n`
  samples excluding the initial `x0` (#P2-26). Pre-allocated result
  vector. Added input validation.

### Structural Changes

- Added comprehensive input validation to all functions.
- Abstracted 8 repetitive internal functions in `rMARS` into two factory
  helpers (`make_tangent_seg`, `make_secant_seg`) with a `switch` dispatcher.
- Updated `DESCRIPTION` to use `Authors@R` field.
- Bumped version to 1.2.0.


## Version 1.1.1 (2018-03-21)

- Initial CRAN release.
