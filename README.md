# Package to read Sibelius-DARK simulation data

Reads the catalog information from HBT, Velociraptor and GALFORM outputs from the Sibelius-DARK simulation.

Installation via:

```bash
git clone https://github.com/stuartmcalpine/pysibeliusdark.git
cd pysibeliusdark
python3 -m pip install .
```

## Reading GALFORM galaxies

Output data in the GALFORM directory is expected to be split over multiple `ivol_XXX` folders and `ivol_XXX_mags` directories.

In the MPI reading case, each rank reads its own subset of the `ivol_XXX` directories. These arrays can then be reduced to rank=0 calling `gather_galaxies()` (see example below).

### Input params to ``read_galform()``

| Input | Description | Is optional? | Default option |
| ----- | ----------- | --------- | ------- | 
| data_dir | path to GALFORM data (base directory containing the `ivol_XXX` directories) | No | - |
| num_files | number of files GALFORM data is split over (number of `ivol_XXX` directories) | No | - |
| output_no | output/snapshot number (int, 1 is redshift=0) | No | - | 
| comm= | MPI4PY communicator | Yes | None |
| verbose= | True for more stdout output | Yes | False |
| ivol_list= | Option to read only a subset of ivols (list of ints) | Yes | None |
| to_convert_out_h= | True to convert out h-values | Yes | True |

### Example usage (No MPI case)

```python
import pysibeliusdark.galform as galform

# Set up read_galform object.
data_dir = "/path/to/galform/folder/"
num_files = 1024
output_no = 1 # z=0
g = galform.read_galform(data_dir, num_files, output_no)

# Load galaxy data.
what_to_load = ["SubhaloID", "xgal", "ygal", "zgal"]
g.load_galaxies(what_to_load)

# Link Sibelius-DARK specific properties (ra, dec, vr, etc).
g.link_sibelius()

# Data gets stored in the data dict. Access the data:
print(g.data)
```

### Example usage (MPI case)

```python
from mpi4py import MPI
import pysibeliusdark.galform as galform

# MPI communicator.
comm = MPI.COMM_WORLD

# Set up read_galform object.
data_dir = "/path/to/galform/folder/"
num_files = 1024
output_no = 1 # z=0
g = galform.read_galform(data_dir, num_files, output_no, comm=comm)

# Load galaxy data.
what_to_load = ["SubhaloID", "xgal", "ygal", "zgal"]
g.load_galaxies(what_to_load)

# Link Sibelius-DARK specific properties (ra, dec, vr, etc).
g.link_sibelius()

# Reduce all galaxies to rank 0.
g.gather_galaxies()

# Data gets stored in the data dict. Access the data:
# Only rank 0 will have any data after gather.
print(g.data["xgal"])
```

## Compute Sibelius-DARK properties

When loading galaxies/subhaloes from HBT, Velociraptor or GALFORM the data gets stored in a ``data`` dict in the loader object (see examples above). To append Sibelius-DARK specific quantities (such as ra, dec, z) use the ``link_sibelius()`` function.

This appends the ``data`` dict with:

| Attribute in dict | Description |
| ----- | ----------- |
| distance : ndarray | Distance from each galaxy to observer (MW) |
| ra_deg, dec_deg : ndarray | Right ascention and declination of galaxies in degrees |
| ra_rad, dec_rad : ndarray | Right ascention and declination of galaxies in radians |
| vt, vr : ndarray | Tangential and radial velocities of galaxies relative to observer |
| mag_*_app : ndarray | Apparent magnitude of galaxies (if you loaded the absolute magnitudes |
| l_deg, b_deg : ndarray | Galactic longitude and latitude of galaxies in degrees |
| l_rad, b_rad : ndarray | Galactic longitude and latitude of galaxies in radians |
