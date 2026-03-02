# TRUDI 2026 -- TrixiShallowWater tutorial

This provides a copy of the slides regarding the TrixiShallowWater features
as well as the files needed to setup an interesting bottom topography and flooding run.

TODO: Add details and such

## What a bottom topography!

This will be the flooding of Cologne downtown, but we do not want to spoil the surprise!
So keep these descriptions vague about files to download and instructions describing the steps.

## What a Cartesian run!

The chuck of downtown Cologne is on a box, so we can initially setup the simulation using
`TreeMesh`. It would be interesting to use a wave maker style BC along the Rhine to flood
the region near Koeln Hbf.

## What a curvilinear run!

A straight-sided mesh is a bit boring! Have them download the boundary data and make a mesh
using HOHQMesh. This will be the Trixi witch outline (resized to encompass most of downtown Koeln).
I designed and adjusted the scaling do that it looks like the cathedral is the witch's nose.
Rerun the simulation just for fun. This isn't a *real* test case after all.
