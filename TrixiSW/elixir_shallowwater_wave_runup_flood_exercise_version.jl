
# Include packages
using TrixiBottomTopography
using JuMP
using HiGHS
using OrdinaryDiffEqSSPRK
using Trixi
using TrixiShallowWater


# Point to the two dimensional bottom data
bathymetry_data = joinpath(@__DIR__, "dgm_merged.txt")

# Shape preserving L1 spline interpolation of the data
spline_struct = LaverySpline2D(bathymetry_data)
spline_func(x::Float64, y::Float64) = spline_interpolation(spline_struct, x, y)

# H0 has an elevation 41 m
equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 41.0,
                                    threshold_desingularization = 1e-8)

function initial_condition_well_balanced(x, t, equations::ShallowWaterEquations2D)
    # # Calculate primitive variables
    v1 = 0
    v2 = 0
    b = spline_func(x[1], x[2])

    # It is mandatory to shift the water level at dry areas to make sure the water height h
    # stays positive. The system would not be stable for h set to a hard zero due to division by h in
    # the computation of velocity, e.g., (h v) / h. Therefore, a small dry state threshold
    # with a default value of 500*eps() ≈ 1e-13 in double precision, is set in the constructor above
    # for the ShallowWaterEquations and added to the initial condition if h = 0.
    H = max(equations.H0, b + equations.threshold_limiter)

    return prim2cons(SVector(H, v1, v2, b), equations)
end

# This function defines the time dependent elevation for a soliton wave
function flood_wave(t, equations)
    # Parameters for a soliton type wave
    A = 100
    scal = 10 * 0.25 / (1 * sqrt(log(2)))
    kx = 1 / sqrt(30)
    t0 = 26.0
    arg = ( kx * (t - t0) ) / scal
    eta = A * sech(arg)^2
    return eta
end


function boundary_function_flood_wave(x, t, equations::ShallowWaterEquations2D)
   b = spline_func(x[1], x[2])
   H = max(equations.H0, b + equations.threshold_limiter)
   v1 = 0
   v2 = 0
   
   # TODO: Add the `flood_wave` elevation to the water level `H` in the area of the river mouth.
   return prim2cons(SVector(H, v1, v2, b), equations)
end

# Setting initial condition
initial_condition = initial_condition_well_balanced

# Setting the boundary conditions
boundary_condition = (; Bottom = BoundaryConditionDirichlet(initial_condition),
                        Left = BoundaryConditionDirichlet(initial_condition),
                        Right = BoundaryConditionDirichlet(initial_condition),
                        # TODO: Introduce a Dirichlet-type boundary condition that prescribes the
                        # `boundary_condition_flood_wave`.
                        Top = BoundaryConditionDirichlet(initial_condition))

###############################################################################
# Get the DG approximation space

volume_flux = (flux_wintermeyer_etal, flux_nonconservative_wintermeyer_etal)
surface_flux = (FluxHydrostaticReconstruction(flux_hll_chen_noelle,
                                              hydrostatic_reconstruction_chen_noelle),
                flux_nonconservative_chen_noelle)

basis = LobattoLegendreBasis(7)

indicator_sc = IndicatorHennemannGassnerShallowWater(equations, basis,
                                                     alpha_max = 0.5,
                                                     alpha_min = 0.001,
                                                     alpha_smooth = true,
                                                     variable = waterheight)
volume_integral = VolumeIntegralShockCapturingHG(indicator_sc;
                                                 volume_flux_dg = volume_flux,
                                                 volume_flux_fv = surface_flux)

solver = DGSEM(basis, surface_flux, volume_integral)

###############################################################################
# Get the mesh

mesh_file = Trixi.download("https://gist.githubusercontent.com/andrewwinters5000/18cad369281c6a91ab3e2a8303a4e500/raw/c3a55010a13147fea822220c1c890e3932918203/CartesianBox.inp",
                           joinpath(@__DIR__, "CartesianBox.inp"))
mesh = P4estMesh{2}(mesh_file)

# create the semi discretization object
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
                                    boundary_conditions = boundary_condition)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 60.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
                                     save_analysis = false)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

save_solution = SaveSolutionCallback(dt = 3.0,
                                     save_initial_solution = true,
                                     save_final_solution = true)

stepsize_callback = StepsizeCallback(cfl = 0.7)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        stepsize_callback)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

###############################################################################
# run the simulation

# use a Runge-Kutta method with error based time step size control
sol = solve(ode, SSPRK43(stage_limiter!);
            ode_default_options()...,
            callback = callbacks,
            adaptive = false, dt = 1.0);