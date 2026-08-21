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
| beam top | `y = −50` |

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
