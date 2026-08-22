# cad/ — foot-model geometry

Source geometry for the three jerboa foot models. Publication schematics are
**rendered from these STLs** by `scripts/fig_foot_schematic.m`, never drawn by
hand and never AI-generated, so every figure is reproducible from the CAD.

## Expected files

| file | model |
|---|---|
| `jerboa_foot_model_tighter_toe_spacing.stl` | Tight |
| `jerboa_foot_model_rectangularbeam.stl` | Default |
| `jerboa_foot_model_wider_toe_spacing.stl` | Wide |

The three differ **only** in toe splay. Beam, linkage and mount are identical,
which is why the bare (projected) intruding area is the same for all three and
only the convex hull grows with splay.

`scripts/fig_foot_schematic.m` errors naming the missing path if any file is
absent. It never substitutes, approximates or synthesizes geometry — a figure
that silently drew something other than the real model would be worse than no
figure.

> **Note:** a copy of `jerboa_foot_model_rectangularbeam.stl` already exists in
> `src/`, left over from earlier work. It is the Default model. Move it here
> rather than exporting a fresh one if it is the same revision; the scripts read
> only from `cad/`.

## Provenance

Exported from the OpenSCAD models. Re-export rather than editing an STL: the
`.scad` sources are the definition, and these are a build product.

## STL coordinate convention

Established by the OpenSCAD models and relied on by the schematic script:

| axis | meaning |
|---|---|
| **−Y** | penetration / drop direction. The foot travels toes-first along −Y. |
| **XZ** | the toe-splay plane. Splay differences between models appear here. |
| toe tips | `y ≈ −61.2` |
| top of the rectangular beam / foot-bar junction | `y = −50` (the highest point of the toe structure is `y = −50.17`; the inclined bar, marker post and mount start here) |

### Depth landmarks

Re-derived from each STL by `local_landmarks` in `scripts/fig_foot_schematic.m`,
which asserts them against these values to 0.05 mm on every run, so a label can
never drift from the CAD. `z` is toe-tip depth below the bed surface,
`z = y − (−61.192)`.

| landmark | `y` (mm) | `z` (cm) |
|---|---|---|
| toe tips | −61.19 | 0.00 |
| rectangular beam: underside | −53.20 | 0.80 |
| rectangular beam: top / cut | −50.00 | 1.12 |
| marker post: top | −21.42 | 3.98 |
| mount: bottom | −2.50 | 5.87 |

The beam top is the default cut plane (`'CutY'`, default `-50`). Everything at
`y <= CutY` is what enters the bed, so the projected area below that plane is
the intruding footprint.

Two consequences worth stating, since both are easy to get wrong:

- Projection for the **footprint** is onto **XZ** — normal to the drop axis.
  Projecting onto XY instead would give the silhouette, not the intruding area.
- Projection for the **assembly side view** is along **Z** onto **XY**, which is
  the view that shows the drop direction. That figure uses the *full* mesh, with
  no cut.

## Reference areas

Computed independently (trimesh/shapely, 2026-08-20) at `CutY = -50`, and
asserted by the script at 2% relative tolerance on every run:

| model | `A_bare` (cm²) | `A_hull` (cm²) |
|---|---|---|
| Tight | 2.122 | 2.607 |
| Default | 2.122 | 3.495 |
| Wide | 2.122 | 4.052 |

`A_bare` is identical across models by construction. If an assertion fails, the
STLs and these reference values disagree — check whether the CAD was revised
before trusting either number. See `docs/QUANTITIES.md` for the definitions.

## Reference depth sweep

`foot_area_vs_depth_reference.csv` is an independent sweep of the same three
models over toe-tip depth (trimesh/shapely, 2026-08-22): `A_bare`, `A_hull` and
the local cross-section `a_local`, in cm², at each cut plane. It is a
**verification file only** — the MATLAB script never reads it at runtime. The
reference values it contains are transcribed into `scripts/fig_foot_schematic.m`,
which asserts its own depth sweep against them at `z = 6, 10, 20, 29 mm` with a
2% relative tolerance. The script's own sweep is written to
`foot_area_vs_depth.csv` in the figures folder and plotted in
`fig_foot_area_vs_depth`.
