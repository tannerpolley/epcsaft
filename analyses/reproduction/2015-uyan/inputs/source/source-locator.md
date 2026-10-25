# Source locator

## Uyan et al. 2015 model

Local source: `/home/tnnrpolley21/Workspaces/Engineering/MEA-Thermodynamics/docs/papers/md/Uyan et al.md`

SHA-256: `3bec38e73260af8b820d460e914a748a0a58cce77659069ea89938a7692f5347`

- reactions and activity conventions: lines 234-337
- reaction constants and starting extents: Tables 1-2, lines 319-376
- molecular and ion parameters: Tables 3-4, lines 397-435
- binary parameters: Table 5 and Eqs. 26-27, lines 439-494
- finite reacting liquid and neutral volatile VLE: lines 533-558

The source is internally inconsistent about water-MDEAH+ dispersion: Section
4.2.1 says zero, Table 5 says Eq. 27, and the VLE results section says it was
inherited from water-MDEA. The executable reconstruction follows Table 5 and
the results section because those statements describe the reported VLE.

## Ermatchkov et al. 2006 observations

DOI: `10.1021/ie0604270`

Local PDF SHA-256: `54d60f67736ac44633078fc39677eee75521cd79f2eed734d5fb514d9462850a`

Local Mathpix companion SHA-256:
`f826ebcb00084f36225e1b711c1a959fbe7be970a09b9ad624f3c6fdeac1cd5a`

Exact observations are materialized from Tables 1-3 at companion lines
47-95, 100-135, and 148-185. The retained CSV contains all 101 rows and the
reported composition uncertainty, pressure uncertainty, and repeatability
standard deviation. The Uyan Figure 6 comparison uses Ermatchkov data in its
approximately 0.19 MDEA mass-fraction lane; the executable target subset is
therefore the 30 rows whose CO2-free solvent MDEA mass fraction is between
0.18 and 0.20.

## Hsieh et al. 2007 permittivity support

DOI: `10.1021/je600515j`

The six pure-MDEA relative-permittivity rows are copied from the official NIST
ThermoML record. Their least-squares line is
`epsilon_MDEA = 52.905597142857 - 0.102457142857 T/K` with a maximum source-row
deviation of 0.0112. Uyan applies the Hsieh correlation above the 303.15-328.15
K measurement range; the 353.15 and 393.15 K campaign values retain that
source extrapolation explicitly.

Uyan cites the pressure-dependent Floriano--Nascimento water correlation but
does not publish the evaluated correlation or values used in its calculations.
The executable replay therefore uses and explicitly labels the Archer--Wang
polynomial approximation already materialized in the Wangler transfer
campaign. It must not be described as an exact implementation of the cited
Floriano--Nascimento correlation.
