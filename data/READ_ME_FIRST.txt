WARNING: subhalo_scf/ directory has 11523 files in it. They are the SCF matrices of all 11523 subhalos, as found in the z=0 snapshot of VL-2 (see Ngan et al, 2015).

The subhalos' information is in subhalo_catalog.dat. The columns are
1: ID
2: Virial mass (solar masses)
3: x (kpc)
4: y (kpc)
5: z (kpc)
6: vx (km/s)
7: vy (km/s)
8: vz (km/s)
9: Virial radius (kpc)
10: Scale radius (kpc)

Each subhalo is identified by its ID in the catalogue. You may notice that the IDs are not sequential and have gaps between numbers. This is because these subhalos are only the first-level subhalos (ie. subhalos that are immediate children of the main VL-2 halo, and not sub-subhalos). The ID can be used to connect to the SCF file inside subhalo_scf/.

read_scf.py is a simple python tool to display the contents of an SCF coefficient file.

