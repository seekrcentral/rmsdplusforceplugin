from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

import rmsdplusforceplugin

prmtop = AmberPrmtopFile("alanine_dipeptide.prmtop")
inpcrd = AmberInpcrdFile("alanine_dipeptide.inpcrd")

system = prmtop.createSystem(
    implicitSolvent=OBC2, implicitSolventKappa=1.0/nanometer, 
    constraints=HBonds)

energy_string = "k*(RMSD - value)^2"
ref_positions = inpcrd.positions
#particle_selection = list(range(22))
#rmsd_force = RMSDForce(ref_positions, particle_selection)
alignment_particle_selection = list(range(11))
rmsd_particle_selection = list(range(11, 22))

rmsd_force = rmsdplusforceplugin.RMSDPlusForce(
    ref_positions[0], alignment_particle_selection, rmsd_particle_selection)

#rmsd_force = RMSDForce(
#    ref_positions, alignment_particle_selection)

exit()

cv_force = CustomCVForce(energy_string)
cv_force.addCollectiveVariable("RMSD", rmsd_force)
cv_force.addGlobalParameter("k", 1.0)
cv_force.addGlobalParameter("value", 0.0*angstrom)

system.addForce(cv_force)

# Create a Langevin integrator for constant temperature.
integrator = LangevinIntegrator(300.0*kelvin, 1/picosecond, 0.002*picoseconds)

simulation = Simulation(prmtop.topology, system, integrator)

# Uncomment the following line to create a simulation that uses one of the GPUs.
#platform = Platform.getPlatformByName('CUDA')
#properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
#simulation = Simulation(prmtop.topology, system, integrator, platform, 
#                        properties)

# Assign atomic positions
simulation.context.setPositions(inpcrd.positions)

# Set the system's atomic velocities to a random distribution defined by a 
#  Maxwell-Boltzmann distribution.
simulation.context.setVelocitiesToTemperature(300.0*kelvin)
    
# Minimize system energy
simulation.minimizeEnergy()

# Create a reporter to write the trajectory to a file.
#simulation.reporters.append(PDBReporter('alanine_dipeptide_traj.pdb', 1000))

# Create another reporter to display system info as the simulation runs.
simulation.reporters.append(StateDataReporter(
    stdout, 100, step=True, potentialEnergy=True, temperature=True, 
    volume=True))

# Advance time by 10000 timesteps (40 picoseconds for Langevin integrator).
simulation.step(1000)
