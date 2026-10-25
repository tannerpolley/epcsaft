# PC-SAFT superancillary coefficients

Bell and Deiters, *Industrial & Engineering Chemistry Research* 62(4), 1958–1967 (2023),
doi:10.1021/acs.iecr.2c02916 (Zotero `BEC7SFNS`). Coefficients are from
https://github.com/usnistgov/SAFTsuperanc, commit
`a3aaf0c39f78b053ef12a97339449140f2c7b0f8` (version 0.0.7).

The retained JSON is a lossless conversion of these upstream files:

| File | SHA-256 |
|---|---|
| `implementation/domain.cpp` | `54b192e11712460760862c87527db46b765542ec9e8d726e8d9d332d5d764907` |
| `implementation/PCSAFTsuperancillary_crit.hpp` | `9070912c382d09b1d1dc7ddf694c63429f40b9b22881d73bac35be1983b56788` |
| `LICENSE.md` | `8788d7bda501b931d0e138e46b4fa545afeff3a9428132ad75d5620094c3219c` |
| `pcsaft.json` | `dba6e908ec067008f9bffa147ea02cc0bbf81543ae615cf5c2a0a613844b0384` |

Regenerate with the upstream checkout at that commit:

```sh
python3 engine/data/superancillary/convert.py /path/to/SAFTsuperanc engine/data/superancillary/pcsaft.json a3aaf0c39f78b053ef12a97339449140f2c7b0f8
```

Every coefficient passes through a binary64 value and its round-trip decimal
representation. The 166,526 upstream numeric values and retained values agree
bit for bit. Meson generates hexadecimal C++ literals in its build directory;
no calculation reads a file at runtime. The generator checks interval and
coefficient dimensions and contiguous temperature segments.

The table uses $w=1/m$, $\tilde T=T/(\varepsilon/k)$ and
$\tilde\rho=\rho N_A\sigma^3$ (molecular density, not segment density).
The two lower-temperature coefficients define
$\tilde T_{min}=\exp(a)m^b\tilde T_c$.
The algorithm and its narrower start domain are owned by ALG-EQUILIBRIUM.

`LICENSE-NIST.md` is the verbatim upstream notice. NIST employee work is in the
United States public domain; use elsewhere is permitted when the notice and
warranty disclaimer accompany copies. No paper text or figures are copied.
The table is source data, outside the runtime line budget. The wheel is not
distributed; distribution would require including the notice with its licence files.
