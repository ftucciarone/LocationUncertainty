# Anisotropic diffusion, div(a grad) form

Diffusion of tracers and momentum by the six-component variance tensor `var_ten`. Not built or tested here. It targets the NEMO interface of 2022 (`un`, `ua`, `tsb`, `tsa`), so it needs porting before use.

**Origin:** old `develop` branch, commit `9cc82cc` (2022-03-09), folder `srcNEMO`. Files are copied unchanged. The whole tree is in tag `archive-develop-2022`.

| Term | Tracer (`tracer/tlutradiff.F90`) | Momentum (`momentum/`) |
|---|---|---|
| Horizontal: a11, a12, a22 | `tlu_trahhdiff` | `tlu_hhdiff` in `tluhdyn.F90` |
| Cross: a13, a23 | `tlu_trahzdiff`, **commented out at the call** | `tlu_hzdiff` in `tluhzdyn.F90` |
| Vertical: a33 | `tlu_trazzdiff` | `tlu_zzdiff` in `tluzdyn.F90` |

The momentum files also hold the matching drift routines (`tlu_hhdrift`, `tlu_hzdrift`, `tlu_zzdrift`).

## Support files
- `support/tlu_step.F90`: calling order, in `tlu_bcdyn` and `tlu_tradyn`.
- `support/tlu.F90`: declares `var_ten` and the index map `ndiffidx`.
- `support/tlu_noi_crt.F90`: module `tlusvd`, which builds `var_ten`.
- `support/tlutrasbc.F90`: `tlu_tradflux`, the surface flux used by `tlu_trazzdiff`.

`tlu_step.F90` also uses `tlubias`, `tlucmp`, `tlusbc`, `tlustp`, `tlutraadv` and, through `tlu_noi_crt.F90`, `tludgns`. These are not diffusion and are only in the archive tag.

## Check before reuse
- The routines carry the warning "This code has some strange index in variance tensor". It was never resolved.
- `tlu_step.F90:269` comments out `tlu_trahzdiff`, so tracers ran without the a13 and a23 terms.
- A commit of 2022-03-01 says "Bottom boundary fluxes killed for stability".
- `tlu_noi_crt.F90` has the TODO "Make sure var_ten is defined on the T grid".
