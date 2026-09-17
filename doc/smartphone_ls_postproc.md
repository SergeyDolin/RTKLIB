# Smartphone relative post-processing by batch least squares

This note defines the target architecture for a Bernese-like smartphone GNSS
post-processor in this RTKLIB tree. The goal is a separate post-processing
engine, not another tuning pass inside the existing sequential relative Kalman
filter.

## Motivation

The current relative solution is implemented in `src/rtkpos.c` by `relpos()`.
It forms zero-difference residuals, converts them to double differences, and
then updates the state with `filter()` epoch by epoch. Smartphone-specific
logic already exists there (`nf == 6` means logical L1/L5 using physical slots
0/2, C/N0-dependent variances, receiver-reported `Lstd/Pstd`, robust residual
gating, and partial ambiguity resolution), but the solution is still driven by
a sequential Kalman state.

For noisy smartphone data this is fragile: cycle slips, duty cycling, channel
dependent phase biases, multipath, and wrong early ambiguities can pollute the
state before enough redundancy has accumulated. A post-processing engine should
instead build observation arcs first, screen them, estimate all selected epochs
and ambiguities together, then fix only defensible ambiguity subsets.

## Processing mode

Initial implementation target:

- Relative base-rover processing from RINEX observation files.
- Static and stop-and-go first; kinematic can be added after the static engine
  is stable.
- Multi-GNSS and multi-frequency, using all available code and phase
  observables selected by RTKLIB signal priority.
- Uncombined double-difference model by default. Ionosphere is estimated or
  constrained, not eliminated unless the user explicitly chooses an IF
  combination.
- Float solution, robust reweighting, ambiguity subset fixing with LAMBDA, and
  fixed solution recomputation.

## Proposed new files

- `src/postls.c`: batch least-squares relative solver.
- `src/postls.h` or exported declarations in `src/rtklib.h`: public API.
- `app/consapp/sppostls/`: standalone command line program, similar to
  `rnx2rtkp`, for testing and operational use.
- `test/postls/`: small synthetic and RINEX regression cases.

## High-level pipeline

1. Read rover/base observations and navigation products through existing
   RTKLIB readers.
2. Synchronize epochs and interpolate base residuals only when the time offset
   is defensible.
3. Build per-satellite/per-frequency phase arcs. Split arcs on LLI, missing
   phase, half-cycle flags, geometry-free jumps, Doppler/TDCP jumps, large code
   inconsistencies, or clock events.
4. Create zero-difference modeled observables using existing satellite
   position, clock, antenna, tide, troposphere, and ionosphere helpers.
5. Form single and double differences using a stable reference satellite per
   constellation/frequency/arc block. Reference selection must prefer long,
   high-elevation, high-C/N0, low-residual arcs.
6. Assemble weighted normal equations for all accepted code and phase
   observations in the session/window.
7. Solve iteratively by weighted least squares. Use robust weighting between
   iterations and remove only confirmed gross outliers.
8. Resolve ambiguity subsets with LAMBDA. For smartphone data, require arc
   continuity, C/N0/standard-deviation gates, residual agreement, ratio test,
   and post-fit validation.
9. Recompute the fixed solution using fixed ambiguity constraints. Output
   coordinates, covariance, residuals, rejected observations, fixed ambiguity
   report, and quality flags.

## State vector

For a static baseline, the first state should be small and observable:

- Rover position correction, 3 parameters.
- Optional zenith wet delay correction for rover, 0 or 1 parameter.
- Optional slant ionosphere parameters per satellite arc, for longer baselines.
- Carrier phase ambiguity per valid DD arc and frequency.
- Optional smartphone phase-bias nuisance terms when integer behavior is not
  reliable for a device/signal.

For kinematic support later:

- Rover position per epoch or per knot, with optional smoothness constraints.
- Ambiguities remain arc constants between detected slips.
- TDCP/Doppler can be used as constraints, but not as a Kalman dynamic model.

## Observation equations

Default short-baseline uncombined DD model:

- Code DD residual includes geometry, troposphere, ionosphere sign by
  frequency, code bias residual, and measurement noise.
- Phase DD residual includes the same geometry terms plus wavelength times
  integer or float ambiguity.
- Satellite and receiver clocks cancel in DD, but receiver clock events still
  matter during arc construction.

Longer baselines should switch from "ionosphere nearly cancels" to one of:

- Estimated slant ionosphere per satellite arc with constraints.
- External IONEX/CODE ionosphere constraints.
- Quasi-ionosphere-free ambiguity strategy for dual-frequency data.

## Smartphone-specific reliability rules

- Prefer receiver-reported standard deviations when present, but use them as a
  floor/inflator rather than blindly trusting them.
- Weight by C/N0 and elevation. C/N0 should dominate for smartphone antennas.
- Treat L5/E5 code as potentially cleaner than L1 code, but do not assume
  phase integer behavior without validation.
- Detect and model common receiver clock jumps separately from per-satellite
  slips.
- Avoid global fix-and-hold. Post-processing should fix, validate, and
  recompute; it should not feed a wrong fix forward.
- Allow float/decimeter output when integer ambiguity is not defensible. A
  reliable float solution is better than a false centimeter solution.

## Milestones

1. Static float DD batch LS on code and phase, no ambiguity fixing.
2. Robust weighting and residual/outlier reports.
3. Arc builder with smartphone slip detection and clock-jump separation.
4. LAMBDA ambiguity subset fixing plus fixed recomputation.
5. Multi-frequency strategy: L1/L5 first, then all RTKLIB-selected frequencies.
6. Kinematic/windowed batch LS with position knots and TDCP constraints.
7. Benchmark suite against current `rnx2rtkp`, Justin/Bernese results if
   available, and known-base smartphone data.

## Accuracy target

The realistic first target is repeatable decimeter-level static/stop-and-go
relative positioning with smartphone internal antenna data in open sky, and
centimeter-level only when the antenna/device data actually preserves usable
phase arcs. The solver must report when the data support only float ambiguities
instead of forcing a misleading fixed result.
