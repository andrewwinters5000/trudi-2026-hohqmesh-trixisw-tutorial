
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

# H0 based on Cologne has elevation 37 m and the Rhine is about 4 m deep
equations = ShallowWaterEquations2D(gravity = 9.81, H0 = 41.0,
                                    threshold_desingularization = 1e-8)

function initial_condition_wave(x, t, equations::ShallowWaterEquations2D)
    # # Calculate primitive variables
    v1 = 0
    v2 = 0
    b = spline_func(x[1], x[2])
    H = max(equations.H0, b + equations.threshold_limiter)

    return prim2cons(SVector(H, v1, v2, b), equations)
end


function boundary_condition_wave_maker(u_inner, orientation, direction, 
                                       x, t,
                                       surface_flux_functions, equations::ShallowWaterEquations2D)
    surface_flux_function, _ = surface_flux_functions

    # Use an exponential window function to put the tsunami in the river mouth
    # Parameters for a soliton type wave
    A = 100
    scal = 10 * 0.25 / (1 * sqrt(log(2)))
    kx = 1 / sqrt(30)
    t0 = 26.0
    arg = ( kx * (t - t0) ) / scal
    L = 357250
    w = 225
    window = exp(-((x[1] - L) / w)^2)
    eta = window * A * sech(arg)^2
    H_ext = equations.H0 + eta

    # Create the external water height and internal reference height
    h_ext = max(equations.threshold_limiter, H_ext - u_inner[4])
    h0 = max(equations.threshold_limiter, equations.H0 - u_inner[4])

    # Compute the incoming velocity as in Eq. (10) of the paper
    # "A limiter-based well-balanced discontinuous Galerkin method for shallow-water flows
    #  with wetting and drying: Triangular grids" by Vater, Beisiegel, and Behrens
    v2_ext = 2 * (sqrt(equations.gravity * h_ext) - sqrt(equations.gravity * h0))

    # Create the external solution state in the conservative variables
    u_outer = SVector(h_ext, zero(eltype(x)), -h_ext * v2_ext, u_inner[4])

    # Calculate the boundary flux
    flux = surface_flux_function(u_inner, u_outer, orientation, equations)

    # Return the conservative and nonconservative fluxes.
    # The nonconservative part is zero as we assume a constant bottom topography at the boundary.
    return (flux, zero(u_inner))
end

# Setting initial condition
initial_condition = initial_condition_wave

# Setting the boundary conditions
boundary_condition = (; y_neg = BoundaryConditionDirichlet(initial_condition),
                        x_neg = BoundaryConditionDirichlet(initial_condition),
                        x_pos = BoundaryConditionDirichlet(initial_condition),
                        y_pos = boundary_condition_wave_maker)

###############################################################################
# Get the DG approximation space

# Use entropy conservative fluxes in the volume and entropy stable fluxes with hydrostatic reconstruction
# on the surface to ensure positivity and stability in the presence of wetting and drying.
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

coordinates_min = (spline_struct.x[1], spline_struct.y[1])
coordinates_max = (spline_struct.x[end], spline_struct.y[end])

mesh = TreeMesh(coordinates_min, coordinates_max,
                initial_refinement_level = 6,
                n_cells_max = 10_000, periodicity = false)

# Create the semi discretization object
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

save_solution = SaveSolutionCallback(dt = 1.0,
                                     save_initial_solution = true,
                                     save_final_solution = true)

time_series = TimeSeriesCallback(semi, [(357000.0, 5.645e6)])

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
                        analysis_callback,
                        alive_callback,
                        save_solution,
                        time_series,
                        stepsize_callback)

stage_limiter! = PositivityPreservingLimiterShallowWater(variables = (waterheight,))

###############################################################################
# run the simulation

# use a Runge-Kutta method with error based time step size control
sol = solve(ode, SSPRK43(stage_limiter!);
            ode_default_options()...,
            callback = callbacks,
            adaptive = false, dt = 1.0);
