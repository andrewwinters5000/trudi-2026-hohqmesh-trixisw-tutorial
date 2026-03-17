# TRUDI 2026 -- HOHQMesh tutorial

## Agenda

Broadly, the HOHQMesh portion of the tutorial is structured as:

1. Overview of HOHQMesh strategy and control file structure
2. Building a first control file
3. Read-in with HOHQMesh.jl and generate with HOHQMesh.jl
4. Interactive nozzle mesh creation.

## HOHQMesh.jl Quick Cheat Sheet

For easier reference during the HOHQMesh.jl interactive session
here is a cheat sheet of the basic functionality.
```
   p = newProject(<projectName>, <folder>)
   plotProject!(p, options)
   updatePlot!(p, options)
   generate_mesh(p)
   remove_mesh!(p)

   c = new(name, startLocation [x,y,z], endLocation [x,y,z])   *Straight Line*
   c = new(name, center [x,y,z], radius, startAngle, endAngle) *Circular Arc*
   c = new(name, xEqn, yEqn, zEqn)                             *Parametric equation*
   c = new(name, dataFile)                                     *Spline with data from a file*
   c = new(name, nKnots, knotsMatrix)                          *Spline with given knot values*
   r = newRefinementCenter(name, type, center, gridSize, radius)
   r = newRefinementLine(name, type, startPoint, endPoint, gridSize, width)

   add!(p, c)                      *Add outer boundary curve*
   add!(p, c, <InnerBoundaryName>) *Add curve to an inner boundary*
   add!(p, r)                      *Add refinement region*

   addBackgroundGrid!(p, [top, left, bottom, right], [nX, nY, nZ]) *No outer boundary*
   addBackgroundGrid!(p, [dx, dy, dz])                             *If an outer boundary is present*

   crv          = getCurve(p, curveName)               *Get a curve in the outer boundary*
   crv          = getCurve(p, curveName, boundaryName) *Get a curve in an inner boundary*
   index, chain = getChain(p, boundaryName)            *Get a complete inner boundary curve*
   r            = getRefinementRegion(p, name)

   removeOuterBoundary!(p)                    *Entire outer boundary curve*
   removeInnerBoundary!(p, innerBoundaryName) *Entire inner boundary curve*
   remove!(p, name)                           *Curve in outer boundary*
   remove!(p, name, innerBoundaryName)        *Curve in inner boundary*
   removeRefinementRegion!(p, name)
```
