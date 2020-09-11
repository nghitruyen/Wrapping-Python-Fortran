import numpy as np
import os
import matplotlib.pyplot as plt

import dassflow1d
import dassflow1d.m_control as m_control
import dassflow1d.m_mesh as m_mesh
import dassflow1d.m_sw_mono as m_sw_mono
import dassflow1d.m_obs as m_obs
from dassflow1d.post.results import plot_BZQ


def run_direct():

    #----------------------------------------------------------------------------------------------
    # MESH SETUP
    #----------------------------------------------------------------------------------------------
    # Read mesh
    mesh = dassflow1d.read_mesh("mesh.geo")
    # Set Strickler type (power_law_h:K=alpha*h^beta)
    mesh.set_strickler_type("powerlaw_h")
    # Set uniform parameters for Strickler
    mesh.set_uniform_strickler_parameters([25.0, 0.0])
    
    #----------------------------------------------------------------------------------------------
    # MODEL SETUP
    #----------------------------------------------------------------------------------------------
    # Initialise model from mesh
    model = m_sw_mono.Model(mesh)
    # Set upstream boundary condition (id:inflow discharge)
    model.bc[0].id = "discharge"
    # Set inflow timeseries
    model.bc[0].set_timeseries(t=[0, 1800, 3600, 5400, 7200], y=[1.0, 2.0, 5.0, 3.0, 1.0])
    # Set downstream boundary condition (id:normal_depth)
    model.bc[1].id = "normal_depth"
    # Set scheme
    model.set_scheme("preissmann")
    
    #----------------------------------------------------------------------------------------------
    # SIMULATION SETTINGS
    #----------------------------------------------------------------------------------------------
    # Start time
    model.ts = 0
    # End time
    model.te = 7200
    # Computation timestep
    model.dt = 60
    # Results timestep
    model.dtout = 600
    # Minimum water depth
    model.heps = 0.001
    # Implicit coefficient for Preissmann scheme
    model.theta_preissmann = 0.667

    #----------------------------------------------------------------------------------------------
    # DIRECT RUN
    #----------------------------------------------------------------------------------------------
    # Run time loop
    model.time_loop()
    
    #----------------------------------------------------------------------------------------------
    # POST-PROCESSING
    #----------------------------------------------------------------------------------------------
    # Retrieve field 'x' (distance from outlet) for segment 0
    x0 = mesh.get_segment_field(0, "x")
    # Retrieve field 'bathy' (bathymetry) for segment 0
    b0 = mesh.get_segment_field(0, "bathy")
    # Retrieve field 'z' (water surface height) for segment 0
    z0 = model.dof.get_segment_field(mesh, 0, "z")
    # Retrieve field 'q' (discharge) for segment 0
    q0 = model.dof.get_segment_field(mesh, 0, "q")
    # Create and show plot
    plot_BZQ(x0, z0, q0, time_label="final", bathy=b0)


# To run tests without pytest (debug)
if __name__ == "__main__":
    run_direct()
