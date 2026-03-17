
# TODO: (1) Add file with endpoints for the joining of nozzles and end point lines
#       (2) Put this and the nozzle.control into the repo
#       (3) Add the manifest to the repository

using HOHQMesh, GLMakie

# make a new project
p = newProject("nozzle", "out")

# Reset params
setPolynomialOrder!(p, 4)
setPlotFileFormat!(p, "sem")

# Set the background
addBackgroundGrid!(p, [0.5, 0.5, 0.0])

# Outer boundary for this example
nozzTop = new("nozzTop_spline", joinpath(@__DIR__, "Data", "NozzleTop.txt"))
add!(p, nozzTop)

inlet = newEndPointsLineCurve("inlet", [-1.2610563, 1.6411425, 0.0], [-1.4184389, 1.0537827, 0.0])
add!(p, inlet)

nozzBot = newSplineCurve("nozzBot_spline", joinpath(@__DIR__, "Data", "NozzleBottom.txt"))
addCurveToOuterBoundary!(p, nozzBot)

ext1 = newEndPointsLineCurve("ext1", [0.0, 0.72965925, 0.0], [1.0, 0.5, 0.0])
ext2 = newEndPointsLineCurve("ext2", [1.0, 0.5, 0.0], [1.0, 0.0, 0.0])
sym_bndy = newEndPointsLineCurve(":symmetry", [1.0, 0.0, 0.0], [9.0, 0.0, 0.0])
right = newEndPointsLineCurve("Right", [9.0 , 0.0, 0.0], [9.0, 5.0, 0.0])
Top = newEndPointsLineCurve("Top", [9.0 , 5.0, 0.0], [0.12940952, 5.0, 0.0])
Left = newEndPointsLineCurve("Left", [0.12940952, 5.0, 0.0], [0.12940952, 1.2126222, 0.0])

add!(p, ext1)
add!(p, ext2)
add!(p, sym_bndy)
add!(p, right)
add!(p, Top)
add!(p, Left)

# add refinement region (only line first)
refineLine = newRefinementLine("refine_line", "sharp",
                                [-1.5, 0.0, 0.0],
                                [3.0, 0.0, 0.0],
                                0.1,
                                2.0)
add!(p, refineLine)

# add refinement center to remove triangle like element
refineCent1 = newRefinementCenter("refine_center", "sharp",
                                 [0.15, 1.2, 0.0],
                                 0.05,
                                 0.25)
add!(p, refineCent1)


# plot project command
plotProject!(p, MODEL + GRID + REFINEMENTS)
