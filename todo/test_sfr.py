import sibelius.galform as galform

#comm = MPI.COMM_WORLD
d="/cosma7/data/dp004/jch/SibeliusOutput/Sibelius_200Mpc_1/Galform/models/Lacey16_icool4_vcut_table/output/"
d="/cosma7/data/dp004/jch/SibeliusOutput/Sibelius_200Mpc_1/Galform/models/Lacey16/output/"
g = galform.read_galform(d, 10, 1, comm=comm)

g.load_galaxies(['xgal', 'ygal', 'zgal', 'mag_SDSS-g_o_tot'])
g.link_sibelius(compute_distance=True)
g.gather_galaxies()
print(g.data)

