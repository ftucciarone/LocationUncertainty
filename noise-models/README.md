# Noise models

The noise models under test, one file per model. A model is added here only when it is chosen for
testing; the older pipelines in `py-kits/noisekit/pipelines` that never get here will be dropped. No
model is here yet.

## Rules for a model
- [ ] Reads every input through `datakit`, with no file paths or dataset-specific names in the model.
- [ ] Takes preprocessed fields (`preprocessing/`) as a source. It does not compute them itself.
- [ ] Runs one source, one calculation and one target grid.
- [ ] Takes data locations from `py-kits/paths.toml`.
- [ ] Flips the sign of velocity means at the north fold, but not of variances or masks.
- [ ] Writes provenance into its output (`datakit.provenance`).
- [ ] Writes the file NEMO reads under the names the Fortran expects.
- [ ] Is documented: idea, formulas, code, discretisation, run, output, status.
- [ ] Is validated: its results look good to F. Tucciarone.

Fields it computes on the way are written with `datakit.write_field`, so they read back as a source.

## Models

| Model | Source | Calculation | Target | Validated |
|---|---|---|---|---|
| | | | | |
