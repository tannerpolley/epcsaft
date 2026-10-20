# Contributing

This page is being prepared. Build and check commands are in the
[repository README](https://github.com/tannerpolley/ePC-SAFT/blob/main/README.md).

## Optional process-example tools

The documentation build uses the Engine with MUMPS. The optional upstream IDAES
Ipopt executable includes HSL; selecting MUMPS does not remove that code.
Obtain solver tools directly from upstream. Preserve their notices and keep
these external build tools outside the checkout, Engine builds, documentation
outputs and shared caches. Documentation uploads contain rendered pages,
notebooks, figures, numerical records and text logs only.

HSL, a collection of Fortran codes for large-scale scientific computation. See https://www.hsl.rl.ac.uk/.
