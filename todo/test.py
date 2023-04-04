from mpi4py import MPI
import sibelius.galform as galform
import sibelius.velociraptor as velociraptor

#comm = MPI.COMM_WORLD
#
#d="/cosma7/data/dp004/jch/SibeliusOutput/Sibelius_200Mpc_1/Galform/models/Lacey16/output/"
#g = galform.read_galform(d, 10, 1, comm=comm)
#
#g.load_galaxies(['xgal', 'ygal', 'zgal', 'mag_SDSS-g_o_tot'])
#g.link_sibelius(compute_distance=True)
#g.gather_galaxies()
#print(g.data)

#import sibelius.hbt as hbt
#from mpi4py import MPI
#
#comm = MPI.COMM_WORLD
#h = hbt.read_hbt_subhaloes('/cosma6/data/dp004/rttw52/SibeliusOutput/Sibelius_200Mpc_1/hbt_refactored/', comm=comm)
#
#h.load_haloes(199, what_to_load=['Mbound'])
#h.gather_haloes()
#print(h.data)

prop_file = '/cosma8/data/dp004/rttw52/Sibelius/Sibelius_200Mpc_512_phase2/vr/stf'
x = velociraptor.read_velociraptor(prop_file)
print(x.HEADER)
x.load_galaxies(['SO_Mass_200_rhocrit', 'Xc', 'Yc', 'Zc'])
x.link_sibelius(compute_distance=True)
print(x.data)

