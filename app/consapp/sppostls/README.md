# sppostls

`sppostls` is a smartphone-oriented static relative post-processor using
batch least squares over double-difference code and carrier phase observations.

## Build

```sh
cd /Users/sergeidolin/RTKLIB/app/consapp/sppostls/gcc
make sppostls
```

## Run

The rover observation file must be first, the base observation file second, and
navigation files after them.

```sh
cd /Users/sergeidolin/RTKLIB/app/consapp/sppostls/gcc

./sppostls -k phone_static.conf -o sol.pos \
  ../test/GEOP1895.26o \
  ../test/NSK100RUS_R_20261890000_01D_30S_MO.rnx \
  ../test/BRDC00IGS_R_20261890000_01D_MN.rnx \
  ../test/IGS0OPSFIN_20261890000_01D_15M_ORB.SP3
```

The default test config uses the NSK1 base position:

```text
ant2-postype = llh
ant2-pos1    = 55.0122550194
ant2-pos2    = 82.9850248083
ant2-pos3    = 141.165
```

Use `-r X Y Z` for ECEF base coordinates or `-l lat lon h` to override the
base position from the config on the command line.

The default config uses broadcast navigation because the bundled test SP3 files
are GPS-only while the phone data are multi-GNSS. Use `-precise` when a complete
precise orbit/clock set is available for the enabled systems.

For kinematic processing use `phone_kine.conf`:

```sh
./sppostls -k phone_kine.conf -o kine.pos \
  ../test/kine/phone.26o \
  ../test/kine/NSK100RUS_R_20262020000_06H_30S_MO.rnx \
  ../test/kine/BRDC00IGS_R_20262020000_01D_MN.26N \
  ../test/kine/COD0MGXFIN_20262020000_01D_05M_ORB.SP3 \
  ../test/kine/COD0MGXFIN_20262020000_01D_30S_CLK.CLK
```

Kinematic mode is batch least squares as well: each epoch gets its own rover
XYZ unknowns, while carrier-phase ambiguities are shared by continuous arcs.
The solver runs this as 60 s windows so 1 Hz rover files can be processed
without creating one very large normal equation. Rover/base epochs are paired to
the nearest base epoch within `pos2-maxage`, which allows a 1 Hz rover file to
work with a 30 s base RINEX. The first implementation is float-only and does not
yet add velocity or smoothness constraints between epochs.

ANTEX files from `file-satantfile` and `file-rcvantfile` are now loaded by
`sppostls`. Receiver antenna corrections are applied only when
`pos1-posopt2 = on`; satellite PCV corrections are applied only when
`pos1-posopt1 = on`. The default phone config keeps both off because the current
base coordinates are supplied explicitly and this test set performs better
without applying the RINEX antenna height as an extra offset.

For smartphone raw L1/L5 processing, the solver now uses the code observations
for the first geometry pass and automatically reduces their weight on later
iterations so the carrier phase dominates the final static estimate. Setting
`pos1-ionoopt = dual-freq` enables an experimental ionosphere-free L1/L5 path,
but it requires enough satellites with simultaneous L1 and L5 code/phase at both
rover and base.
