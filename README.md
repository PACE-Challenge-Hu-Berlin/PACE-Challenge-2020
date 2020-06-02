# PACE'20 Solver for Treedepth

## Building

The code can be compiled using the Meson build system:

```
mkdir build && cd build
meson --buildtype=debugoptimized ..
```

## Running

Run the solver via: `./td --solver simple-pid INSTANCE` where `INSTANCE` is a `.gr` file in PACE format.


## PACE 2020 submission

The PACE 2020 version of the solver (i.e., the binary that was uploaded to optil.io), was built using the
flags `-Dtd-static=true -Dtd-avoid-madvfree=true` for `meson`. The solver is invoked by `./td --solver simple-pid --witness /dev/stdin`.
