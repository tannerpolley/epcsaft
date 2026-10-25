| Numerical verification | Retained evidence |
| --- | --- |
| Base-wheel reduction | Omitted, empty and full default maps against the merged-base wheel (main 1218dc78): 12 declared cases with first T/P/EOS actions, bitwise equality |
| Native column reduction | 16 raw/common-molality cases against the merged base; values, all lower and higher columns, both terminal changes and availability, bitwise equality |
| Convention switch | Density target maximum error 7.99360578e-15; mole-fraction maximum error 1.80411242e-15, each against 1e-10 |
| CO₂ pure-vapor identities | Three source temperatures: value error ≤1e-10; scaled T/P slope error ≤1e-9 |
| Ideal-gas limit at 0.01 Pa | 4.37889724e-10 / 1.30016198e-09 / 4.37888725e-10 for L, T and P limits; each ≤1e-7 |
| Refined first/second derivative maximum scaled error | 2.82875635e-07 against 1e-6 |
| Refined third/fourth derivative and caloric maximum scaled error | 6.42057101e-06 against 1e-5 |
| Fixed-source system-pressure columns | Exact zeros |
| Vapor branch at ethane saturation | 180 K; pure-gas source slope matches vapor volume/(RT), distinct from liquid volume/(RT) |
| Declared branch admission | Ethane at 180 K and 3 MPa declared pure gas: single liquid-side root, refused with species 1 and cause; accepted before the fix. Above the EOS critical temperature the single root serves either branch |
| Caloric pure override | 260 K, 3 MPa system and 5 MPa source: s1 declared pure liquid, the only native root; order-3/4 rows unchanged because the earlier pure-gas declaration had silently used this root |
| Failure meaning | Charged pure selection and pure density-domain failure retain their cause; a refused T slope keeps the value and feed direction |
| Mixed common molality | Native refusal: mixed species references require raw mole fractions; common molality is not qualified |
| Key solves | case-6/figure-2-point-7/cold: (16, 53, 17), work identical, state difference 0; case-7/flash-5918.5286-Pa/cold: (43, 57, 44), work identical, state difference 0 |
