import numpy as xp
import constant_dictionaries as phys
import scipy.ndimage
import os
# import pytom.basic.functions

# image display
# import matplotlib.pyplot as plt
# from pylab import *
# matplotlib.use('Qt5Agg')

V_WATER = 4.5301 # potential value of low density amorphous ice (from Vulovic et al., 2013)

def extend_volume(vol, increment, pad_value=0, symmetrically=False, true_center=False):
    """

    @param vol:
    @param increment: list with increment value for each dimension
    @param pad_value: Float, value to use as padding
    @param symmetrically: Boolean, if False (default) the volume is just padded with zeros.
    @param true_center: Boolean, if True interpolate to true center
    @param filter: Boolean, prefilter for scipy interpolation
    @return:
    """
    if symmetrically:
        if true_center:
            vol = xp.pad(vol, tuple([(0, x) for x in increment]), 'constant', constant_values=pad_value)
            # TODO Use voltools for this as well!
            # Filter removes potential artifacts from the interpolation
            return scipy.ndimage.interpolation.shift(vol, tuple([x / 2 for x in increment]), order=3)
        else:
            from pytom.tompy.tools import paste_in_center
            new = xp.zeros([a + b for a, b in zip(vol.shape, increment)])
            if pad_value:
                new += pad_value
            return paste_in_center(vol, new)
    else:
        return xp.pad(vol, tuple([(0, x) for x in increment]), 'constant', constant_values=pad_value)


def create_gaussian_low_pass(size, radius, center=None):
    """Create a 3D mask with gaussian edges.

    @param size: size of the resulting volume.
    @param cutoff: radius of the sphere inside the volume.
    @param center: sigma of the Gaussian.
    @param gpu: center of the sphere.

    @return: sphere inside a volume.
    """
    assert len(size) == 3

    if center is None:
        center = [size[0]/2, size[1]/2, size[2]/2]

    [x,y,z] = xp.mgrid[0:size[0], 0:size[1], 0:size[2]]
    r = xp.sqrt((x-center[0])**2+(y-center[1])**2+(z-center[2])**2)

    filter = xp.exp(-(r ** 2) / (2 * radius ** 2))

    return filter


def reduce_resolution(map, voxel_size, resolution):
    """
    Accepts only square input maps to easy low-pass filtering

    @param map: Volume
    @type map: 3d numpy array
    @param voxel_size: Voxel size of the volume in A
    @type voxel_size: float
    @param to_resolution: Desired resolution after applying low-pass filter
    @type to_resolution: float
    @param sharpness: ratio
    @type sharpness: float

    @return: Filtered volume
    @rtype: 3d numpy array

    @author: Marten Chaillet
    """
    from pytom.tompy.transform import fourier_filter

    assert resolution >= voxel_size, "the requested resolution is non-valid as it is smaller than the voxel size"
    assert len(set(map.shape)) == 1, "dimensions of input are not equal"

    # resolution reduction factor
    nr_pixels_fourier = (map.shape[0] * voxel_size ) / resolution
    # create full gaussian mask
    mask = create_gaussian_low_pass(map.shape, nr_pixels_fourier/2)
    # apply mask in fourier space
    result = fourier_filter(map, mask, human=True)
    return result


def low_pass_filter(map, voxel_size, resolution):
    """
    Accepts only square input maps to easy low-pass filtering

    @param map: Volume
    @type map: 3d numpy array
    @param voxel_size: Voxel size of the volume in A
    @type voxel_size: float
    @param to_resolution: Desired resolution after applying low-pass filter
    @type to_resolution: float
    @param sharpness: ratio
    @type sharpness: float

    @return: Filtered volume
    @rtype: 3d numpy array

    @author: Marten Chaillet
    """
    from pytom.tompy.transform import fourier_filter
    from pytom.tompy.tools import create_sphere

    assert resolution > voxel_size, "the requested resolution is non-valid as it is smaller than the voxel size"
    assert len(set(map.shape)) == 1, "dimensions of input are not equal"

    # calculate radius of mask from pixels needed in fourier space for the desired resolution
    nr_pixels_fourier = (map.shape[0] * voxel_size) / resolution
    mask = create_sphere(map.shape, radius=nr_pixels_fourier, sigma=15, num_sigma=10)
    # apply mask in fourier space
    result = fourier_filter(map, mask, human=True)
    return result


def scale_old(potential, voxel_size_in, voxel_size_out):
    """
    Scale volume with scipy ndimage zoom.
    @param potential: 3d array (numpy)
    @param voxel_size_in: Original voxel size in A
    @param voxel_size_out: Desired voxel size in A
    @return:
    """
    # Downsample volume with scipy.ndimage
    factor = voxel_size_in / voxel_size_out
    # VOLTOOLS?
    potential = scipy.ndimage.zoom(potential, factor, order=3)
    return potential


def scale(volume, factor):
    """
    Scale volumes with skimage.
    @param potential: ndimage, numpy array
    @param factor: float, tuple of floats. Single float for same scaling in each dimension. Tuple for different factors
    @return: scaled multidimensional array (numpy)
    """
    # Skimage scale could be better for downsampling than scipy zoom
    import skimage
    # order 3 for splines, preserve_range otherwise image is returned as float between -1 and 1
    return skimage.transform.rescale(volume, factor, mode='constant', order=3, preserve_range=True, multichannel=False,
                                     anti_aliasing=False)

def bin(potential, factor):
    """

    @param potential:
    @param factor: integer value
    @return:
    """
    assert type(factor) is int and factor > 1, print('non-valid binning factor, should be integer above 1')

    s = (potential.shape[0] % factor) // 2
    d = (potential.shape[0] % factor) % 2

    potential = potential[s:potential.shape[0] - s - d, s:potential.shape[0] - s - d, s:potential.shape[0] - s - d]

    ds = int(potential.shape[0] // factor)
    image_size = potential.shape[0]

    binned = potential.reshape(ds, image_size // ds,
                               ds, image_size // ds, ds, image_size // ds).mean(-1).mean(1).mean(-2)

    return binned


def call_chimera(pdb_folder, pdb_id):
    """
    Run chimera for pdb file in order to add hydrogens and add symmetry units. The resulting structure is stored as a
    new pdb file with the extension {id}_symmetry.pdb

    @param pdb_folder: Path to a folder where pdb files are stored
    @type pdb_folder: string
    @param pdb_id: ID of pdb file
    @type pdb_id: string

    @return: Returns the name of file where the new pdb stucture is stored in. This can differ for pdb depending on
    symmetry thus we need to return it.
    @rtype: string

    @author: Marten Chaillet
    """
    print(f' - Calling chimera for {pdb_id}')

    # Command 'sym' in chimera crashes when there is no BIOMT symmetry in the pdb file. We need to make sure sym is only
    # executed when the BIOMT information specifies symmetrical units.
    symmetry = []
    try:
        with open(f'{pdb_folder}/{pdb_id}.pdb','r') as pdb:
            line = pdb.readline().split()
            while line[0] != 'REMARK':
                line = pdb.readline().split()
            while line[0] == 'REMARK':
                if line[1] == '350' and len(line) > 3:
                    if 'BIOMT' in line[2]:
                        symmetry.append(int(line[3]))
                line = pdb.readline().split()
    except Exception as e:
        print(e)
        raise Exception('Could not read pdb file.')

    print(f'{pdb_id} has {len(set(symmetry))} symmetrical {"unit" if len(set(symmetry)) == 1 else "units"}.')

    extension = 'py'
    if len(set(symmetry)) > 1:
        scriptname = f'_rem-solvent_sym_{pdb_id}'
        try:
            with open(f'{pdb_folder}/{scriptname}.{extension}', 'w') as chimera_script:
                chimera_script.write(f'# Open chimera for {pdb_id} then execute following command:\n'
                                     f'# (i) remove solvent and ions (ii) add hydrogens (iii) add symmetry units\n'
                                     f'from chimera import runCommand as rc\n'
                                     f'rc("open {pdb_folder}/{pdb_id}.pdb")\n'
                                     f'rc("delete solvent")\n'
                                     f'rc("delete ions")\n'
                                     # f'rc("addh")\n'
                                     f'rc("sym group biomt")\n'             # group biomt is also the default
                                     f'rc("combine all modelId 10")\n'
                                     f'rc("write format pdb #10 {pdb_folder}/{pdb_id}_rem-solvent_sym.pdb")\n'
                                     f'rc("stop")\n')
        except Exception as e:
            print(e)
            raise Exception('Could not create chimera script.')
    else:
        scriptname = f'_rem-solvent_{pdb_id}'
        try:
            with open(f'{pdb_folder}/{scriptname}.{extension}', 'w') as chimera_script:
                chimera_script.write(f'# Open chimera for {pdb_id} then execute following command:\n'
                                     f'# (i) remove solvent and ions (ii) add hydrogens (iii) add symmetry units\n'
                                     f'from chimera import runCommand as rc\n'
                                     f'rc("open {pdb_folder}/{pdb_id}.pdb")\n'
                                     f'rc("delete solvent")\n'
                                     f'rc("delete ions")\n'
                                     # f'rc("addh")\n'
                                     f'rc("write format pdb #0 {pdb_folder}/{pdb_id}_rem-solvent.pdb")\n'
                                     f'rc("stop")\n')
        except Exception as e:
            print(e)
            raise Exception('Could not create chimera script.')
    # module chimera should be loaded here...
    try:
        os.system(f'chimera --nogui --script {pdb_folder}/{scriptname}.{extension}')
    except Exception as e:
        print(e)
        raise Exception('Chimera is likely not on your current path.')

    if len(set(symmetry)) > 1:
        return f'{pdb_id}_rem-solvent_sym' # returns new pdb name
    else:
        return f'{pdb_id}_rem-solvent'


def modify_structure_file(filename, pattern, replacement, line_start=''):
    """
    Function required to make pqr files with large negative coordinates readible for APBS. APBS can only parse
    columns in the file when they are properly separated by white spaces. Input file will be overwritten.

    @param filename: File path of pqr type file (with extension)
    @type filename: string
    @param pattern: Pattern to be replaced
    @type pattern: string
    @param replacement: Replacement string for pattern
    @type replacement: string
    @param line_start: Keyword argument, only modify line starting with line_start
    #type line_start: string

    @return: empty
    @rtype:

    @author: Marten Chaillet
    """
    from tempfile import mkstemp
    from shutil import move, copymode

    try:
        # Create temp file
        fh, abs_path = mkstemp()
        with os.fdopen(fh, 'w') as new_file:
            with open(filename) as old_file:
                for line in old_file:
                    if line_start=='':
                        new_file.write(line.replace(pattern,replacement))
                    elif line.split()[0] == line_start:
                        new_file.write(line.replace(pattern,replacement))
                    else:
                        new_file.write(line)
        # Copy the file permissions from the old file to the new file
        copymode(filename, abs_path)
        # Remove original file
        os.remove(filename)
        # Move new file
        move(abs_path, filename)
    except Exception as e:
        print(e)
        raise Exception('Unsuccessful in modifying pqr file with white space delimiters.')
    return


def call_apbs(pdb_folder, structure, apbs_folder, force_field='amber', ph=7.):
    """
    Calls external programs pdb2pqr and apbs to execute on pdb structure. References:

    @param pdb_folder: Folder where pdb structures are stored
    @type pdb_folder: string
    @param structure: Name of pdb file with coordinates of the protein structure
    @type structure: string
    @param apbs_folder: Folder to store output of pdb2pqr and apbs
    @type apbs_folder: string
    @param force_field: Force field for parameterizing atoms (option: amber, ...)
    @type force_field: string
    @param ph: pH value of solvent surrounding the protein
    @type ph: float

    @return: empty, output of programs called is stored in apbs_folder
    @rtype:

    @author: Marten Chaillet
    """
    # pdb2pqr and apbs should be on path for this function to run
    print(f' - Running pdb2pqr and APBS on {pdb_folder}/{structure}.pdb')
    cwd = os.getcwd()
    input_file = f'{pdb_folder}/{structure}.pdb'
    output_folder = f'{apbs_folder}/{structure.split("_")[0]}'
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    output_file = f'{output_folder}/{structure}.pqr'
    apbs_in = f'{output_folder}/{structure}.in'
    try:
        # Also PDB2PKA ph calculation method. Requires PARSE force field, can take very long for large proteins.
        os.system(f'pdb2pqr.py --ff={force_field} --ph-calc-method=propka --with-ph={ph} --apbs-input {input_file} {output_file}')
        print(' - Add white space delimiters to pqr file.')
        modify_structure_file(output_file, '-', ' -', line_start='ATOM')
        # APBS needs to execute from the folder where the structure is present
        os.chdir(output_folder)
        os.system(f'apbs {apbs_in}')
        # os.system(f'module load apbs_mc/1.5; apbs {apbs_in}')
        # ? subprocess.run(['apbs', f'{structure}.in'], cwd=output_folder) # executes in cwd
    except Exception as e:
        print(e)
        raise Exception('pdb2pqr or APBS potentially not on path.')

    # Change back to original directory
    os.chdir(cwd)
    return


def read_structure(filename):

    x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = [], [], [], [], [], []

    filetype = filename.split('.')[-1]

    if filetype == 'pdb':
        try:
            with open(filename, 'r') as pdb:
                lines = pdb.readlines()
                # TODO Whay about HETATM lines?
                atoms = [line for line in lines if line[:4] == 'ATOM']
                for line in atoms:
                    '''
        PDB example
        ATOM   4366  OXT SER I 456      10.602  32.380  -1.590  1.00 53.05           O
                    '''
                    x_coordinates.append(float(line[30:38]))
                    y_coordinates.append(float(line[38:46]))
                    z_coordinates.append(float(line[46:54]))
                    elements.append(line[76:78].strip())
                    b_factors.append(float(line[60:66]))
                    occupancies.append(float(line[54:60]))
        except Exception as e:
            print(e)
            raise Exception('Could not read pdb file.')
    elif filetype == 'pqr':
        try:
            with open(filename, 'r') as pqr:
                lines = pqr.readlines()
                for line in lines:
                    split_line = line.split()
                    # TODO Whay about HETATM lines?
                    if split_line[0] == 'ATOM':
                        '''
            PQR example
            ATOM   5860  HA  ILE   379      26.536  13.128  -3.443  0.0869 1.3870
                        '''
                        x_coordinates.append(float(split_line[5]))
                        y_coordinates.append(float(split_line[6]))
                        z_coordinates.append(float(split_line[7]))
                        elements.append(split_line[2][0])  # first letter of long atom id is the element
                        b_factors.append(0.0)
                        occupancies.append(1.0)
        except Exception as e:
            print(e)
            raise Exception('Could not read pdb file.')
    else:
        print('invalid filetype in iasa_potential, return 0')
        return 0

    return x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies


def iasa_potential(filename, voxel_size=1., oversampling=0, low_pass_filter=True): # add params voxel_size, oversampling?
    """
    interaction_potential: Calculates interaction potential map to 1 A volume as described initially by
    Rullgard et al. (2011) in TEM simulator, but adapted from matlab InSilicoTEM from Vulovic et al. (2013).

    @param filename: ID of pdb file as present in pdb folder
    @type filename: string
    @param voxel_size: Size (A) of voxel in output map, default 1 A
    @type voxel_size: float
    @param oversampling: Increased sampling of potential (multiple of 1), default 1 i.e. no oversampling
    @type oversampling: int

    @return: A volume with interaction potentials
    @rtype: 3d numpy/cupy array[x,y,z], float

    @author: Marten Chaillet
    """
    extra_pixels = 10  # extend volume by 10 A

    print(f' - Calculating IASA potential from {filename}')

    x_coordinates, y_coordinates, z_coordinates, elements, b_factors, occupancies = read_structure(filename)

    x_coordinates = x_coordinates - xp.min(x_coordinates) + extra_pixels
    y_coordinates = y_coordinates - xp.min(y_coordinates) + extra_pixels
    z_coordinates = z_coordinates - xp.min(z_coordinates) + extra_pixels
    # Define the volume of the protein
    szx = xp.abs(xp.max(x_coordinates) - xp.min(x_coordinates)) + 2 * extra_pixels
    szy = xp.abs(xp.max(y_coordinates) - xp.min(y_coordinates)) + 2 * extra_pixels
    szz = xp.abs(xp.max(z_coordinates) - xp.min(z_coordinates)) + 2 * extra_pixels

    if (oversampling > 1) and ( type(oversampling) is int ):
        spacing = voxel_size / oversampling
    else:
        spacing = voxel_size
    print(f'IASA: volume will be sampled at {spacing}A')

    sz = [ x / spacing for x in [szx,szy,szz]] # could increase size for oversampling here
    # sz = [x / voxel_size for x in [szx, szy, szz]]
    # C = 2132.8 A^2 * V; 1E20 is a conversion factor for Angstrom^2
    C = 4 * xp.sqrt(xp.pi) * phys.constants['h']**2 / (phys.constants['el'] * phys.constants['me']) * 1E20

    potential = xp.zeros(tuple(xp.round(sz).astype(int)))
    print(f'#atoms to go over is {len(x_coordinates)}.')

    for i in range(len(elements)):
        if xp.mod(i, 5000) == 0:
            print(f'Potential for atom number {i}.')

        atom = elements[i]
        b_factor = b_factors[i]
        occupancy = occupancies[i]

        sf = xp.array(phys.scattering_factors[atom]['g'])
        a = sf[0:5]
        b = sf[5:10]

        # b is used for calculating the radius of the potential. See Rullgard et al. (2011) for addition of 16 R^2.
        # Incorporates low pass filtering by extending the gaussian radius. NOTE: b_factor is zero anyway
        b += (b_factor) # + 16 * spacing**2)

        r2 = 0
        b1 = xp.zeros(5)
        for j in range(5):
            # Find the max radius over all gaussians (assuming symmetrical potential to 4.5 sigma truncation
            # (corresponds to 10).
            b1[j] = 4 * xp.pi ** 2 / b[j] * spacing ** 2
            r2 = xp.maximum(r2, 10/b1[j])
        # Radius of gaussian sphere
        r = xp.sqrt(r2 / 3)
        xc1, yc1, zc1 = x_coordinates[i] / spacing, y_coordinates[i] / spacing, z_coordinates[i] / spacing
        # Center of the gaussian sphere.
        rc = [xc1, yc1, zc1]
        # Calculate the absolute indexes for the potential matrix.
        kmin = [xp.maximum(0,x).astype(int) for x in xp.ceil(rc-r)]
        kmax = [xp.minimum(xp.floor(x)-1,xp.floor(y+r)).astype(int) for x,y in zip(sz,rc)]
        kmm = max([x-y for x,y in zip(kmax,kmin)])
        # Determine the coordinates for sampling from the sum of Gaussians.
        x = xc1 - xp.arange(kmin[0], kmin[0]+kmm+1, 1)
        y = yc1 - xp.arange(kmin[1], kmin[1]+kmm+1, 1)
        z = zc1 - xp.arange(kmin[2], kmin[2]+kmm+1, 1)

        test = 0
        if test:
            # test in 1D
            print(atom)
            potential_1d = 0
            for j in range(5):
                print(f'{a[j] / (b[j] ** (3 / 2)) * C:.3f}, a_i = {a[j]:.3f}, b_i = {b[j]:.3f}, C = {C:.3f}')
                gaussian = xp.exp(-b1[j] * x ** 2)
                print([f'{c:.2f}' for c in gaussian])
                tmp = a[j] / ( b[j]**(3/2) ) * C * gaussian
                print([f'{c:.2f}' for c in tmp])
                potential_1d += tmp
            print('\n')
            print([f'{c:.2f}' for c in x])
            print([f'{p:.2f}' for p in potential_1d])
            return None

        atom_potential = 0
        for j in range(5):
            x_matrix = xp.tile(xp.exp(-b1[j] * x**2)[:,xp.newaxis,xp.newaxis], [1, kmm+1, kmm+1])
            y_matrix = xp.tile(xp.exp(-b1[j] * y**2)[xp.newaxis,:,xp.newaxis], [kmm+1, 1, kmm+1])
            z_matrix = xp.tile(xp.exp(-b1[j] * z**2)[xp.newaxis,xp.newaxis,:], [kmm+1, kmm+1, 1])
            tmp = a[j] / b[j]**(3/2) * C * x_matrix * y_matrix * z_matrix
            atom_potential += tmp
        atom_potential *= occupancy
        # add the potential of element i to the full potential map
        potential[kmin[0]:kmin[0]+kmm+1, kmin[1]:kmin[1]+kmm+1, kmin[2]:kmin[2]+kmm+1] = \
            potential[kmin[0]:kmin[0]+kmm+1, kmin[1]:kmin[1]+kmm+1, kmin[2]:kmin[2]+kmm+1] + atom_potential

    # Only need to downscale if the volume was oversampled.
    if spacing != voxel_size:
        # Volume needs to be cubic before applying a low pass filter and binning (??)
        size = max(potential.shape)
        increment = [size - a for a in potential.shape]
        # using spline interpolation here does not decrease our accuracy as the volume is already oversampled!
        potential = extend_volume(potential, increment, pad_value=0, symmetrically=True, true_center=True)
        if low_pass_filter:
            potential = reduce_resolution(potential, spacing, voxel_size)
        # potential = bin(potential, oversampling)
        # Other option is to use scale instead of binning, which is proper downsampling.
        potential = scale(potential, 1/oversampling)

    return potential


def parse_apbs_output(filename):
    """
    parse_apbs_output: Parses output file from APBS.
    @param filename: Path and name of .pqr.dx file produced by APBS software
    @return: data, dxnew, dynew, dznew are in A, thickness in m, data is the reshaped apbs output
    @author: Marten Chaillet
    """
    # Read file header
    try:
        with open(filename, 'r') as apbs_file:
            grid, origin, spacing = [],[],[]
            delta_count = 0
            for i in range(11):
                line = apbs_file.readline()
                if 'object 1 class gridpositions counts' in line:
                    grid.extend([int(x) for x in line.split()[5:8]])
                elif 'origin' in line:
                    origin.extend([float(x) for x in line.split()[1:4]])
                elif 'delta' in line:
                    spacing.append(float(line.split()[1+delta_count]))
                    delta_count+=1
            print(f'grid\t= {grid}')
            print(f'origin\t= {origin}')
            print(f'spacing\t= {spacing}')

            line = apbs_file.readline().split()
            data = []
            while line[0] != 'attribute':
                data.extend(line)
                line = apbs_file.readline().split()
    except Exception as e:
        print(e)
        raise Exception('Could not open APBS data file.')

    # Reshape to sequence
    data = xp.array(data,dtype=float)
    # Form to 3D grid
    data = xp.reshape(data,(grid[2],grid[1],grid[0]),order='F').transpose(2,1,0) # x, y, z needs to be done as such to correspond to iasa
    dx = spacing[0]
    dy = spacing[1]
    dz = spacing[2]
    return data, dx, dy, dz


def resample_apbs(filename, voxel_size=1.0, low_pass_filter=False):
    """
    resample_APBS: First calls parse_abps_output to read an apbs output file, then scales voxels to 1 A voxels and
    refactors the values to volts.

    @param filename: Full filenamepath
    @type filename: string

    @return: Resampled and scaled V_bond potential
    @rtype: 3d numpy/cupy array[x,y,z], float

    @author: Marten Chaillet
    """
    # from pytom_volume import vol
    # from pytom_numpy import vol2npy
    from voltools import transform

    print(f' - Parsing and resampling APBS file {filename}')

    # Parse APBS data file
    potential, dxnew, dynew, dznew = parse_apbs_output(filename)
    # Check if the voxel size is allowed, and adjust if neccessary
    smallest_possible = max([dxnew, dynew, dznew])
    if not(voxel_size >= smallest_possible):
        print(f'Requested voxel size is smaller than voxel size of the apbs map. Adjust to smallest possible voxel size '
              f'of {smallest_possible}.')
        voxel_size = smallest_possible
    # Make voxels same size by scaling to dx
    factor = (dxnew/dxnew, dynew/dxnew, dznew/dxnew)
    potential = scale(potential, factor)
    if low_pass_filter:
        size = max(potential.shape)
        increment = [size-a for a in potential.shape]
        potential = extend_volume(potential, increment, pad_value=0, symmetrically=True, true_center=True)
        potential = reduce_resolution(potential, dxnew, voxel_size)
    # Scale the volume to voxels with voxel_size
    potential = scale(potential, dxnew/voxel_size)

    # TODO use voltools

    print(f'Data after reshaping to {voxel_size} A voxels: {potential.shape}')

    # Convert to from kT/e to volts
    temperature = 291 # [K] = 18 C (room temperature)
    convert_to_volts = phys.constants['kb'] * temperature / phys.constants['el']
    potential = convert_to_volts * potential
    return potential


def combine_potential(iasa_potential, bond_potential, voxel_size):
    """
    Combine isolated atom and bond potential.

    @param iasa_potential:
    @type iasa_potential: 3d numpy/cupy array[x,y,z], float
    @param bond_potential:
    @type bond_potential: 3d numpy/cupy array[x,y,z], float

    @return: Combined IASA and bond potential (APBS)
    @rtype: 3d numpy/cupy array[x,y,z], float

    @author: Marten Chaillet
    """
    print(' - Combining iasa and bond potential')
    # Reduce interpolation! Only shift one of the volumes to true center. Preferably V_bond.
    # This code could be more elegant.
    assert max(iasa_potential.shape) >= max(bond_potential.shape), print('Break! Attempting to extend IASA volume.')
    try:
        # Because volumes are shifted to true center (by interpolation) the precision is reduced slightly.
        difference = [max(iasa_potential.shape) - x for x in bond_potential.shape]
        bond_potential = extend_volume(bond_potential, difference, symmetrically=True, true_center=True)
        # Apply filter to remove the artifacts from interpolating
        bond_potential = reduce_resolution(bond_potential, voxel_size, voxel_size)
        full_potential = iasa_potential + bond_potential
    except Exception as e:
        print(e)
        raise Exception('Could not fit atom and bond potential together.')
    return iasa_potential, bond_potential, full_potential


def wrapper(pdb_id, pdb_folder, apbs_folder, iasa_folder, bond_folder, map_folder, voxel_size=1.0, ph=7.0):
    import pytom.tompy.io
    # TODO function can be executed in parallel for multiple structures
    # pdb_folder = '/data2/mchaillet/structures/pdb'
    # apbs_folder = '/data2/mchaillet/structures/apbs'
    # iasa_folder = '/data2/mchaillet/structures/potential_iasa'
    # bond_folder = '/data2/mchaillet/structures/potential_bond'
    # map_folder = '/data2/mchaillet/structures/potential_map'

    # Call external programs for structure preparation and PB-solver
    structure = call_chimera(pdb_folder, pdb_id) # output structure name is dependent on modification by chimera
    call_apbs(pdb_folder, structure, apbs_folder, ph=ph)

    outfile = f'{structure}_ph{ph:.1f}_{voxel_size:.2f}A'
    # Calculate atom and bond potential, and store them
    # 4 times oversampling of IASA yields accurate potentials
    v_atom = iasa_potential(f'{apbs_folder}/{structure.split("_")[0]}/{structure}.pqr', voxel_size=voxel_size,
                            oversampling=4, low_pass_filter=True)
    # extension with even numbers does not require interpolation
    v_atom = extend_volume(v_atom, [10, 10, 10], symmetrically=True, true_center=False)
    pytom.tompy.io.write(f'{iasa_folder}/{outfile}.mrc', v_atom)
    # Volume for creating mask and calculating correlation scores
    pytom.tompy.io.write(f'{iasa_folder}/{outfile}.em', v_atom)

    v_bond = resample_apbs(f'{apbs_folder}/{structure.split("_")[0]}/{structure}.pqr.dx', voxel_size=voxel_size,
                           low_pass_filter=False)
    _, v_bond, map = combine_potential(v_atom, v_bond, voxel_size)
    pytom.tompy.io.write(f'{bond_folder}/{outfile}.mrc', v_bond)
    pytom.tompy.io.write(f'{map_folder}/{outfile}.mrc', map)

    # Bin volume 10 times
    resolution = voxel_size * 10
    map = reduce_resolution(map, voxel_size, resolution)
    map = bin(map, 10)

    outfile = f'{structure.split("_")[0]}_{resolution:.2f}A.mrc'
    pytom.tompy.io.write(f'{map_folder}/{outfile}', map)
    return


if __name__ == '__main__':
    # DEFAULT FOLDERS FOR WRITING INPUT/OUTPUT
    pdb_folder = '/data2/mchaillet/structures/pdb'
    apbs_folder = '/data2/mchaillet/structures/apbs'
    iasa_folder = '/data2/mchaillet/structures/potential_iasa'
    bond_folder = '/data2/mchaillet/structures/potential_bond'
    map_folder = '/data2/mchaillet/structures/potential_map'

    # IN ORDER TO FUNCTION, SCRIPT REQUIRES INSTALLATION OF PYTOM (and dependencies), CHIMERA, PDB2PQR (modified), APBS

    # LIST OF PDBS TO EXECUTE ON
    # pdb_ids = ['3cf3', '1s3x', '1u6g', '4cr2', '1qvr', '3h84', '2cg9', '3qm1', '3gl1', '3d2f', '4d8q', '1bxn']
    pdb_ids = ['6RGQ']

    for id in [id.upper() for id in pdb_ids]:
        if not os.path.exists(f'{pdb_folder}/{id}.pdb'):
            print(f'Skipping {id} because the pdb file does not exist in {pdb_folder}.')
            continue
        # elif os.path.exists(f'{map_folder}/{id}_1.0A.mrc') or os.path.exists(f'{map_folder}/{id}_sym_addh_1.0A.mrc'):
        #     print(f'{id} already has a map in folder {map_folder}.')
        #     continue
        else:
            try:
                wrapper(id, pdb_folder, apbs_folder, iasa_folder, bond_folder, map_folder, voxel_size=0.81)
            except Exception as e:
                print(e)
                print(f'Something when wrong while creating map for {id}. Continuing with next pdb file in list.')
                continue
