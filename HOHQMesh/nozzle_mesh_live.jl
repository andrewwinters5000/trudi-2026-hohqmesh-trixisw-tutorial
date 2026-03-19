
using HOHQMesh, GLMakie

# Make a new project
p = newProject("nozzle", "out")


# Reset params (optional)


# Set the background with bounds [TOP, LEFT, BOTTOM, RIGHT]
# and desired number of elements in each direction
addBackgroundGrid!(p, [6.0, -3.0, -1.0, 9.0], [24, 14, 0])


# Plot project command
plotProject!(p, MODEL + GRID + REFINEMENTS)


# Create the outer boundary for the nozzle (must be counter-clockwise)
# Helpful to first the spline data, then build the curves


# Add refinement regions (if desired)


# Generate the mesh
