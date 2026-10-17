# Reproducing the cortical Poisson workflow figure

The README contains one six-panel figure generated in MATLAB from the bundled left-hemisphere pial and flat surfaces. All scientific computations use the repository's existing functions. The example is a visualization of the supplied geometry, not a patient-control analysis.

## Run

From the repository root, prepare lossless ASCII GIfTI copies with Python 3:

```sh
python3 assets/readme/prepare_demo_surfaces.py
```

Then in MATLAB:

```matlab
addpath(fullfile(pwd, 'assets', 'readme'));
generate_readme_figure;
```

The converter uses Python's standard library and preserves the numeric arrays and GIfTI indexing conventions. It writes copies under `build/readme/` so the bundled reader can load them without Java. The source surface files remain unchanged. Curvature estimation, finite-element assembly, Poisson solves, and figure rendering are performed in MATLAB.

## Figure construction

- **A:** left pial geometry, with 32,492 vertices and 64,980 triangular faces.
- **B:** curvature from `surf2quadratic` and `surf2curvature`, centered using the FEM mass matrix.
- **C:** potential at lambda = 0.1 and face-centroid arrows from `flux_to_display` in face mode. Every 55th face is sampled; arrows are normalized and scaled to 2.4 coordinate units for display.
- **D-E:** potentials at lambda = 1 and 0.01, illustrating changes in spatial detail.
- **F:** the lambda = 0.1 scalar potential displayed by vertex index on the supplied flat mesh, which has the same 32,492 vertices and 59,013 faces. No pial vectors are transferred to the flat mesh.

All fields are computed on the pial mesh using double-precision coordinates. No anatomical region mask is applied. The solver matches each potential's mass-weighted RMS to the centered curvature source. Scalar colors use the common interval [-0.5, 0.5]; values outside the interval saturate. The fractions of vertices outside that interval are 7.50% for curvature and 7.03%, 3.37%, and 2.30% for lambda = 1, 0.1, and 0.01. Lighting provides shape cues on the pial panels; the flat panel is unlit.

## Verified implementation conventions

The underlying linear solve is `(C + lambda*A)*w = A*hc`. The returned `u` includes RMS rescaling and mean-drift removal, so residual checks account for that rescaling. Lambda = 0 is also tested with the mass-weighted zero-mean constraint, although it is not plotted.

At source commit `d144e3a720adb4c565b2c238c287e075f92ad0c3`, the array `Jf` returned by `FEM_poisson_regualarized` has the uphill gradient orientation. This is verified independently from its dot products with triangle edges and the corresponding differences in `u`. The helper `flux_to_display` negates `Jf`, so its face-mode arrows point downhill. Use `-Jf` for the numerical downhill field. The figure does not modify either function.

Results in MATLAB R2026b (26.2.0.3386108):

| Check | Result |
|---|---:|
| Maximum relative scaled-system residual over lambda = 1, 0.1, 0.01, 0 | 6.34e-14 |
| Maximum absolute mass-weighted mean of returned potentials | 4.85e-18 |
| Potential/source RMS ratios | 1 to numerical precision |
| Returned vector versus triangle-edge potential differences | 5.15e-16 relative residual |
| Face-vector tangency | 3.43e-16 relative residual |
| Displayed arrows versus normalized negative returned vectors | 9.23e-17 relative difference |

Exact values are saved in [`numerical-checks.json`](numerical-checks.json). Computed arrays are saved locally in `build/readme/figure-data.mat`. These checks validate this example and its documented conventions; they do not constitute a full audit of all tutorial exercises or the research studies.
