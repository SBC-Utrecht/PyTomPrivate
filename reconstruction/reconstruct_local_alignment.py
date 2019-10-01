#!/usr/bin/env pytom
# -*- coding: utf-8 -*-

"""
Created on Sep 02, 2019

@author: dschulte
"""


def write_subtomograms(particlelist, projection_directory, offset, binning, vol_size, tilt_angles,
                       reconstruction_method='WBP', infr_iterations=1):
    """
    To create subtomograms of all particles in the list

    @param particlelist: list with all particles
    @type particlelist: ParticleList
    @param projection_directory: the directory of the projections
    @type projection_directory: str
    @param offset: the offset
    @type offset: list(int, int, int)
    @param binning: the binningfactor used
    @type binning: int
    @param vol_size: the size of the volume to be reconstructed (in pixels)
    @type vol_size: int
    @param tilt_angles: the tilt angles as belonging to the projections
    @type tilt_angles: list(int)
    @param reconstruction_method: the method used to reconstruct the volumes (choose "WBP" or "INFR")
    @type reconstruction_method: str
    @param infr_iterations: the amount of iterations used when choosing INFR (int)
    @type infr_iterations: int
    @return: nothing, writes the volumes to disk
    @returntype: void
    """
    import os.path
    assert isinstance(projection_directory, str)
    assert os.path.isdir(projection_directory)
    assert isinstance(offset, list) and len(offset) == 3
    assert isinstance(offset[0], int) and isinstance(offset[1], int) and isinstance(offset[2], int)
    assert isinstance(binning, int)
    assert isinstance(vol_size, int)
    assert isinstance(tilt_angles, int)
    assert reconstruction_method == "WBP" or reconstruction_method == "INFR"
    assert isinstance(infr_iterations, int) and infr_iterations > 0

    from math import cos, sin, pi
    from pytom.tompy.transform import cut_from_projection
    from nufft.reconstruction import fourier_2d1d_iter_reconstruct, fourier_2d1d_gridding_reconstruct
    from pytom.tools.ProgressBar import FixedProgBar
    from pytom.tompy.io import read, write
    from pytom.reconstruction.reconstructionStructures import ProjectionList
    import numpy as np

    projection_list = ProjectionList()
    projection_list.loadDirectory(projection_directory)
    projection_list.sort()

    if reconstruction_method == 'WBP':
        # transform the cropping offset
        tmp = projection_list[0]
        sx = tmp.getXSize()  # here should be the size of original projection!
        sy = tmp.getYSize()

        offset[0] = -sx / 2 + offset[0] * binning
        offset[1] = -sy / 2 + offset[1] * binning
        offset[2] = -sx / 2 + offset[2] * binning

        from pytom.basic.structures import PickPosition
        for particle in particlelist:
            pick_position = particle.getPickPosition()
            x = (pick_position.getX() * binning + offset[0])
            y = (pick_position.getY() * binning + offset[1])
            z = (pick_position.getZ() * binning + offset[2])
            particle.setPickPosition(PickPosition(x=x, y=y, z=z))

        # Reconstruct all particles
        projection_list.reconstructVolumes(particles=particlelist, cubeSize=vol_size, binning=1, applyWeighting=True,
                                           showProgressBar=False, verbose=False, preScale=1, postScale=1)

    elif reconstruction_method == 'INFR':
        # NOT WORKING YET some bugs are still in

        print(len(projection_list))
        print("INFR NOT TESTED YET")
        dim_x = projection_list[0].getXSize()  # here should be the size of original projection!
        # dim_y = projection_list[0].getYSize()
        dim_z = dim_x  # make z dim the same as x!

        # reconstruct each particles

        prog = FixedProgBar(0, len(particlelist) - 1, '')

        i = 0
        for m, p in enumerate(particlelist[:1]):
            prog.update(i)
            i += 1

            # transfer the coordinate system
            x, y, z = p.getPickPosition().toVector()
            x = (x + offset[0]) * binning
            y = (y + offset[1]) * binning
            z = (z + offset[2]) * binning

            subregions = []
            n = 0
            for img, ang in zip(projection_list[:], tilt_angles[:]):
                n += 1
                # project the coordinate to 2D image
                yy = y  # assume the rotation axis is around y
                xx = (cos(ang * pi / 180) * (x - dim_x / 2) - sin(ang * pi / 180) * (z - dim_z / 2)) + dim_x / 2

                # cut the small patch out
                patch = cut_from_projection(img, [xx, yy], [vol_size, vol_size])
                patch = patch - np.mean(patch)

                # fill in the subregion
                subregions.append(patch)

            # reconstruct INFR
            v = fourier_2d1d_iter_reconstruct(subregions, tilt_angles, infr_iterations)
            write("Subtomograms/{:s}/particle_{:d}.em".format(particlelist.getFilename(), m), v)


def local_alignment(projections, vol_size, binning, offset, tilt_angles, particle_list_filename, projection_directory,
                    mpi, projection_method='WBP', infr_iterations=1, create_graphics=False, create_subtomograms=False,
                    averaged_subtomogram=False, number_of_particles=-1, skip_alignment=False, particle_list_name="",
                    start_glocal=True, fsc_path=""):
    """
    To polish a particle list based on (an) initial subtomogram(s).

    @param projections: a list with filenames of projections
    @type projections: list(str)
    @param vol_size: the size of the volume to build the new subtomogram in (in pixels)
    @type vol_size: int
    @param binning: the binning factor used
    @type binning: int
    @param offset: the offset used (x, y, z)
    @type offset: list(int, int, int)
    @param tilt_angles: the list of tiltangles used
    @type tilt_angles: list(int)
    @param particle_list_filename: the filename of the particlelist
    @type particle_list_filename: str
    @param projection_directory: the directory of the projections
    @type projection_directory: str
    @param projection_method: the projection method used when create_subtomograms is True
               (possible values: WBP and INFR)
    @type projection_method: str
    @param infr_iterations: the amount of iterations to use when INFR projection is used
    @type infr_iterations: int
    @param create_graphics: to create plots of major parts of the algorithm, mainly used for debugging
               and initial creation
    @type create_graphics: bool
    @param create_subtomograms: to flag if initial subtomograms of the particles should be made (takes a long time!)
    @type create_subtomograms: bool
    @param averaged_subtomogram: to give a path to an averaged subtomogram to be used instead of subtomograms of all
               particles separately (False for off, otherwise a string)
    @type averaged_subtomogram: bool or str
    @param number_of_particles: to use a subset of the particles for the particle polishing
    @type number_of_particles: int
    @param skip_alignment: skips the alignment phase, does not do particle polishing
    @type skip_alignment: bool
    @return: nothing, it writes everything to disk
    @returntype: void
    """
    assert number_of_particles == -1 or number_of_particles > 0
    assert infr_iterations > 0
    assert binning > 0
    assert vol_size > 0
    assert vol_size % 2 == 0
    assert projection_method == "WBP" or projection_method == "INFR"
    assert isinstance(projections, list)
    assert isinstance(vol_size, int)
    assert isinstance(binning, int)
    assert isinstance(offset, list) and len(offset) == 3
    assert isinstance(offset[0], int) and isinstance(offset[1], int) and isinstance(offset[2], int)
    assert isinstance(tilt_angles, list)
    assert isinstance(particle_list_filename, str)
    assert isinstance(projection_directory, str)
    assert isinstance(infr_iterations, int)
    assert isinstance(create_graphics, bool)
    assert isinstance(create_subtomograms, bool)
    assert isinstance(averaged_subtomogram, str) or averaged_subtomogram is False
    assert isinstance(number_of_particles, int)
    assert isinstance(skip_alignment, bool)

    import numpy as np
    import glob

    # load particle list
    from pytom.basic.structures import ParticleList

    particlelist1 = ParticleList()
    particlelist1.fromXMLFile(particle_list_filename)

    particlelist = ParticleList()

    for p in particlelist1:
        if particle_list_name in p.getFilename():
            particlelist.append(p)

    if number_of_particles != -1:
        particlelist = particlelist[:number_of_particles]

    print(len(particlelist))

    # If needed create all subtomograms
    if create_subtomograms:
        write_subtomograms(particlelist, projection_directory, offset, binning, vol_size, tilt_angles,
                           projection_method, infr_iterations)

    # Read all subtomograms
    subtomograms = []
    if not averaged_subtomogram:
        subtomograms = [f for f in glob.glob("Subtomograms/{:s}/*".format(particle_list_name))]

    print("{:s}> Creating the input array".format(gettime()))

    # from pytom.tools.ProgressBar import FixedProgBar
    # progressBar = FixedProgBar(0, len(input_to_processes), 'Particle volumes generated ')
    # progressBar.update(0)

    input_to_processes = []
    particle_number = -1
    for particle in particlelist:
        particle_number += 1
        if averaged_subtomogram:
            subtomogram = averaged_subtomogram
        else:
            subtomogram = subtomograms[particle_number]

        rot = (particle.getRotation().getZ1(), particle.getRotation().getX(), particle.getRotation().getZ2())
        # getX() Z1() Z2() (phi: z1, the: x, psi: z2) test forwards or backwards

        # loop over tiltrange, take patch and cross correlate with reprojected subtomogram
        for img, ang in zip(projections, tilt_angles):
            input_to_processes.append([ang, subtomogram, offset, vol_size, particle.getPickPosition().toVector(), rot,
                                       particle.getFilename(), particle_number, binning, img, create_graphics, fsc_path])

    print(len(input_to_processes))
    print("{:s}> Created the input array".format(gettime()))

    results_file = "local_alignment_results_{:s}.txt".format(particle_list_name)

    if not skip_alignment:
        print("{:s}> Started on running the process".format(gettime()))

        output = mpi.parfor(run_single_tilt_angle_unpack, input_to_processes)

        from pytom.basic.datatypes import fmtLAR, headerLocalAlignmentResults, LOCAL_ALIGNMENT_RESULTS
        np.savetxt(results_file, np.array(output, dtype=LOCAL_ALIGNMENT_RESULTS), fmt=fmtLAR,
                   header=headerLocalAlignmentResults)

        print("{:s}> Ran the processes".format(gettime()))

    if start_glocal:
        run_polished_subtomograms(particle_list_filename, projection_directory, results_file, binning, offset, vol_size, start_glocal)


def run_single_tilt_angle_unpack(inp):
    """
    unpack a list with arguments to "run_single_tilt_angle"
    @param inp: the arguments to "run_single_tilt_angle" in the same order in a single list (iterable)
    @type inp: list/something with indexing
    @return: the value from "run_single_tilt_angle"
    @returntype: list
    """
    return run_single_tilt_angle(
        inp[0], inp[1], inp[2], inp[3], inp[4], inp[5], inp[6], inp[7], inp[8], inp[9], inp[10], inp[11])


def run_single_tilt_angle(ang, subtomogram, offset, vol_size, particle_position, particle_rotation,  particle_filename,
                          particle_number, binning, img, create_graphics=False, fsc_path=""):
    """
    To run a single tilt angle to allow for parallel computing

    @param ang: the tilt angle
    @type ang: int
    @param subtomogram: the filename of the subtomogram
    @type subtomogram: str
    @param offset: the offset used (x,y,z)
    @type offset: list(int, int, int)
    @param vol_size: the size of the volume to be reconstructed (in pixels)
    @type vol_size: int
    @param particle_position: the position of the particle in vector format,
               as given by particle.pickPosition().toVector()
    @type particle_position: tuple
    @param particle_rotation: the rotation of the particle (Z1/phi, X/the, Z2/psi)
    @type particle_rotation: tuple
    @param particle_filename: the filename of the particle, as given by particle.getfilename()
    @type particle_filename: str
    @param particle_number: the number of the particle, to allow for unique mapping
    @type particle_number: int
    @param binning: the binning factor used
    @type binning: int
    @param img: the filename of the projection to be used
    @type img: str
    @param create_graphics: to flag if images should be created for human inspection of the work done
    @type create_graphics: bool
    @return: the newly found positions of the particle, as a list  in the LOCAL_ALIGNMENT_RESULTS format
    @returntype: list
    """
    # assert isinstance(ang, float)
    # assert isinstance(subtomogram, str)
    # assert isinstance(offset, list) and len(offset) == 3
    # assert isinstance(offset[0], int) and isinstance(offset[1], int) and isinstance(offset[2], int)
    # assert isinstance(vol_size, int)
    # assert isinstance(particle_position, tuple)
    # assert isinstance(particle_rotation, tuple)
    # assert isinstance(particle_filename, str)
    # assert isinstance(particle_number, int)
    # assert isinstance(binning, int)
    # assert binning > 0
    # assert isinstance(img, str)
    # assert isinstance(create_graphics, bool)

    print(ang)
    from pytom.tompy.transform import rotate3d
    import numpy as np
    from math import cos, sin, pi
    from pytom.tompy.transform import cut_from_projection
    from pytom.tompy.io import read

    subtomogram = read(subtomogram)
    img = read(img)

    # Get the size of the original projection
    dim_x = img.shape[0]
    # dim_y = img.shape[1]
    dim_z = dim_x  # make z dim the same as x!

    x, y, z = particle_position
    x = (x + offset[0]) * binning
    y = (y + offset[1]) * binning
    z = (z + offset[2]) * binning

    # Get template
    # first rotate towards orientation of the particle, then to the tilt angle
    rotated1 = rotate3d(subtomogram, phi=particle_rotation[0], the=particle_rotation[1], psi=particle_rotation[2])
    rotated2 = rotate3d(rotated1, the=ang)  # 'the' is the rotational axis
    template = rotated2.sum(axis=2)

    # Get coordinates of the paricle adjusted for the tilt angle
    yy = y  # assume the rotation axis is around y
    xx = (cos(ang * pi / 180) * (x - dim_x / 2) - sin(ang * pi / 180) * (z - dim_z / 2)) + dim_x / 2

    # cut the small patch out
    patch = cut_from_projection(img, [xx, yy], [vol_size, vol_size])
    patch = patch - np.mean(patch)
    patch = patch.squeeze()

    # filter using FSC
    fsc_mask = None
    import os
    if os.path.isfile(fsc_path):
        f = open(fsc_path, "r")
        fsc = map(lambda a: float(a), f.readlines())
        f.close()
        fsc_mask = create_fsc_mask(fsc, vol_size)
    elif fsc_path != "":
        print("Not an existing FSC file: " + fsc_path)

    # Cross correlate the template and patch, this should give the pixel shift it is after
    ccf = normalised_cross_correlation_numpy(template, patch, fsc_mask)
    points2d = find_sub_pixel_max_value_2d(ccf)

    x_diff = points2d[0] - vol_size / 2
    y_diff = points2d[1] - vol_size / 2

    if create_graphics:
        m_style = dict(color='tab:blue', linestyle=':', marker='o', markersize=5, markerfacecoloralt='tab:red')
        m_style_alt = dict(color='tab:red', linestyle=':', marker='o', markersize=5, markerfacecoloralt='tab:blue')

        points = find_sub_pixel_max_value(ccf)

        nx, ny, nz = particle_position
        nx += x_diff
        ny += y_diff

        npatch = cut_from_projection(img, [xx + x_diff, yy + y_diff], [vol_size, vol_size])
        npatch = npatch - np.mean(npatch)

        nccf = normalised_cross_correlation_numpy(template, npatch.squeeze())
        npoints = find_sub_pixel_max_value(nccf)
        npoints2d = find_sub_pixel_max_value_2d(nccf)

        import pylab as pp

        grid = pp.GridSpec(3, 3, wspace=0, hspace=0.35, left=0.05, right=0.95, top=0.90, bottom=0.05)

        ax_0_0 = pp.subplot(grid[0, 0])
        ax_0_1 = pp.subplot(grid[0, 1])
        ax_0_2 = pp.subplot(grid[0, 2])
        ax_1_0 = pp.subplot(grid[1, 0])
        ax_1_1 = pp.subplot(grid[1, 1])
        ax_1_2 = pp.subplot(grid[1, 2])
        ax_2_0 = pp.subplot(grid[2, 0])
        ax_2_1 = pp.subplot(grid[2, 1])
        ax_2_2 = pp.subplot(grid[2, 2])

        ax_0_0.axis('off')
        ax_0_1.axis('off')
        ax_0_2.axis('off')
        ax_1_0.axis('off')
        ax_1_1.axis('off')
        ax_1_2.axis('off')
        ax_2_0.axis('off')
        ax_2_1.axis('off')
        ax_2_2.axis('off')

        axis_title(ax_0_0, "Cutout")
        ax_0_0.imshow(patch)
        # np.savetxt("cutout_{:04d}_tiltimage_{:05.2f}.txt".format(particle_number, ang), patch)
        # bigpatch = cut_from_projection(img, [xx, yy], [vol_size*2, vol_size*2])
        # bigpatch = bigpatch - np.mean(bigpatch)
        # bigpatch = bigpatch.squeeze()
        # np.savetxt("bigcutout_{:04d}_tiltimage_{:05.2f}.txt".format(particle_number, ang), bigpatch)
        axis_title(ax_0_1, "Template")
        ax_0_1.imshow(template)
        # np.savetxt("template_{:04d}_tiltimage_{:05.2f}.txt".format(particle_number, ang), template)
        axis_title(ax_0_2, "Shifted Cutout\n(based on cross correlation)")
        ax_0_2.imshow(npatch.squeeze())

        axis_title(ax_1_0, u"Cross correlation\ncutout × template")
        ax_1_0.imshow(ccf)
        ax_1_0.plot([p[1] for p in points], [p[0] for p in points], fillstyle='none', **m_style)
        ax_1_0.plot([points2d[0]], [points2d[1]], fillstyle='none', **m_style_alt)
        ax_1_0.plot([vol_size / 2], [vol_size / 2], ",k")

        ax_1_1.text(0.5, 0.8,
                    "Red: 2D spline interpolation\nx: {:f}\ny: {:f}\nBlue: 1D spline interpolation\nx: {:f}\ny: {:f}"
                    "\nBlack: center".format(
                        x_diff, y_diff, points[0][0] - vol_size / 2, points[0][1] - vol_size / 2), fontsize=8,
                    horizontalalignment='center', verticalalignment='center', transform=ax_1_1.transAxes)

        axis_title(ax_1_2, u"Cross correlation\nshifted cutout × template")
        ax_1_2.imshow(nccf)
        ax_1_2.plot([p[1] for p in npoints], [p[0] for p in npoints], fillstyle='none', **m_style)
        ax_1_2.plot([npoints2d[0]], [npoints2d[1]], fillstyle='none', **m_style_alt)
        ax_1_2.plot([vol_size / 2], [vol_size / 2], ",k")

        axis_title(ax_2_0, u"Zoom into red peak\nin CC cutout × template")
        d = 10
        peak = ccf[int(points2d[0]) - d:int(points2d[0]) + d, int(points2d[1]) - d:int(points2d[1] + d)]
        ax_2_0.imshow(peak)

        axis_title(ax_2_1, u"Zoom into red peak\nin CC cutout × template\ninterpolated")
        ax_2_1.imshow(points2d[2])

        axis_title(ax_2_2, u"Cutout\nGaussian filter σ3")
        import scipy
        ax_2_2.imshow(scipy.ndimage.gaussian_filter(patch, 3))

        pp.savefig("polish_particle_{:04d}_tiltimage_{:05.2f}.png".format(particle_number, ang))

    # progressbar.increment_amount()
    return particle_number, x_diff, y_diff, ang, 0, 0, particle_filename


def axis_title(axis, title):
    axis.text(0.5, 1.05, title, fontsize=8, horizontalalignment='center', transform=axis.transAxes)


def run_polished_subtomograms(particle_list_filename, projection_directory, particle_polish_file, binning, offset,
                              vol_size, start_glocal):
    """
    Reconstructs subtomograms based on a polished particlelist, writes these to the places as specified in particlelist

    @param particle_list_filename: The name of the file of the particlelist
    @type particle_list_filename: str
    @param projection_directory: The name of the directory containing the projections
    @type projection_directory: str
    @param particle_polish_file: The name of the file containing the polished alignment results
    @type particle_polish_file: str
    @param binning: The binning factor
    @type binning: int
    @param offset: The reconstruction offset
    @type offset: list(int, int, int)
    @param vol_size: The size of the particle
    @type vol_size: int
    @return: void
    @returntype: void
    """
    assert isinstance(particle_list_filename, str)
    assert isinstance(projection_directory, str)
    assert isinstance(particle_polish_file, str)
    assert isinstance(binning, int)
    assert isinstance(offset, list) and len(offset) == 3
    assert isinstance(offset[0], int) and isinstance(offset[1], int) and isinstance(offset[2], int)
    assert isinstance(vol_size, int)

    import os.path
    assert os.path.isfile(particle_list_filename)
    assert os.path.isfile(particle_polish_file)
    assert os.path.isdir(projection_directory)

    cwd = os.getcwd()

    first_batchfile = """#!/usr/bin/bash
#SBATCH --time        12:00:00
#SBATCH -N 1
#SBATCH --partition defq
#SBATCH --ntasks-per-node 20
#SBATCH --job-name    polishedReconstruction                                                                     
#SBATCH --error="../LogFiles/%j-polished_subtomograms.err"
#SBATCH --output="../LogFiles/%j-polished_subtomograms.out"
#SBATCH --oversubscribe

module load openmpi/2.1.1 python/2.7 lib64/append pytom/dev/dschulte

cd {:s}

reconstructWB.py --particleList {:s} \
--projectionDirectory {:s} \
--coordinateBinning {:d} \
--size {:d} \
--applyWeighting \
--projBinning 1 \
--recOffset {:d},{:d},{:d} \
--particlePolishFile {:s} \
-n 20"""\
        .format(cwd, particle_list_filename, projection_directory, binning, vol_size, offset[0], offset[1], offset[2],
                particle_polish_file)

    f = open("polished_subtomograms.sh", "w+")
    f.write(first_batchfile)
    f.close()

    import subprocess
    out = subprocess.check_output(['sbatch', 'polished_subtomograms.sh'])
    print("Started reconstruction")

    if start_glocal and out.startswith("Submitted batch job "):
        pid = out.split(" ")[3]

        mask_filename = "/data2/dschulte/BachelorThesis/Data/VPP2/05_Subtomogram_Analysis/Alignment/FRM/test-11-09-19/FRM_mask_200_70_4.mrc"
        jobname = "fsc-filter-27-09-2019"
        iterations = 6
        pixelsize = 2.62
        particleDiameter = 300

        glocal_batchfile = """#!/usr/bin/bash
#SBATCH --time        36:00:00
#SBATCH -N 4
#SBATCH --partition defq
#SBATCH --ntasks-per-node 20
#SBATCH --job-name    pGLocalAlign                                                                       
#SBATCH --error="../LogFiles/%j-GLocal_polished_subtomograms.err"
#SBATCH --output="../LogFiles/%j-GLocal_polished_subtomograms.out"
#SBATCH --dependency=afterok:{:s}

module load openmpi/2.1.1 python/2.7 lib64/append pytom/dev/dschulte

cd {:s}

mpiexec -n 80 pytom GLocalJob.py \
    -p {:s}
    --mask {:s}
    --sphericalMask
    --destination Alignment/GLocal/{:s}/
    --numberIterations {:d}
    --pixelSize {:f}
    --particleDiameter {:d}
    --jobName /Alignment/Glocal/{:s}/job.xml""".format(pid, cwd, particle_list_filename, mask_filename, jobname, iterations, pixelsize, particleDiameter, jobname)
        f = open("glocal_align.sh", "w+")
        f.write(glocal_batchfile)
        f.close()

        if not os.path.isdir(cwd + "/Alignment/GLocal/" + jobname): os.mkdir(cwd + "/Alignment/GLocal/" + jobname)

        out = subprocess.check_output(['sbatch', 'glocal_align.sh'])

        if out.startswith("Submitted batch job "):
            print("Reconstruction and Glocal alignment scheduled")
        else:
            print("Could not start the Glocal alignment script:\n" + out)
            raise Exception("Could not start the Glocal alignment script: " + out)
    elif not start_glocal:
        print("Could not start the reconstruction script:\n" + out)
        raise Exception("Could not start the reconstruction script: " + out)


def normalised_cross_correlation_numpy(first, second, filter_mask=None):
    """
    Do a cross correlation based on numpy

    @param first: The first dataset (numpy 2D)
    @type first: numpy array 2D
    @param second: The second dataset (numpy 2D)
    @type second: numpy array 2D
    @param filter_mask: a filter which is used to filter image 'first'
    @type filter_mask: numpy array 2D
    @return: The cross correlation result
    @returntype: numpy array 2D

    @requires: the shape of first to be equal to the shape of second, and equal t the shape of the filter (if used of course)
    """
    assert first.shape == second.shape
    assert len(first.shape) == 2
    if not(filter_mask is None): assert first.shape == filter_mask.shape

    from numpy.fft import fftshift, ifftn, fftn
    import numpy as np

    prod = 1
    for d in first.shape:
        prod *= d

    if filter_mask is None:
        ffirst = fftn(first)
    else:
        ffirst = fftshift(fftshift(fftn(first)) * filter_mask)

    return np.real(fftshift(ifftn(np.multiply(fftn(second), np.conj(ffirst))))) / prod


def normalised_cross_correlation_cupy(first, second):
    import cupy.fft as nf
    import cupy as np

    prod = 1
    for d in first.shape:
        prod *= d

    return np.real(nf.fftshift(nf.ifftn(np.multiply(nf.fftn(second), np.conj(nf.fftn(first)))))) / prod


def normalised_cross_correlation_mask_numpy(first, second, mask):
    """
    Do cross correlation with a running mask based on numpy

    @param first: The first dataset (numpy 2D)
    @type first: numpy array 2D
    @param second: The second dataset (numpy 2D)
    @type second: numpy array 2D
    @param mask: The mask
    @type mask: numpy array 2D
    @return: The cross correlation result
    @returntype: numpy array 2D

    @requires: the shape of first to be equal to the shape of second and the shape of the mask
    """
    # assert first.shape == second.shape
    # assert first.shape == mask.shape
    # assert len(first.shape) == 2

    import numpy.fft as nf
    import numpy as np

    a = norm_inside_mask(first, mask)
    b = norm_inside_mask(second, mask)

    return np.real(nf.fftshift(nf.ifftn(np.multiply(nf.fftn(b), np.conj(nf.fftn(a)))))) / np.sum(mask)


def norm_inside_mask(inp, mask):
    """
    To normalise a 2D array within a mask

    @param inp: A 2D array to be normalized.
    @type inp: numpy array 2D
    @param mask: A 2D array of the same size as the input to mask of parts of the inp.
    @type mask: numpy array 2D
    @return: A normalized 2D array.
    @returntype: numpy array 2D

    @requires: the shape of inp to be equal to the shape of the mask
    """
    # assert inp.shape == mask.shape
    # assert len(inp.shape) == 2

    import numpy as np

    mea = np.divide(np.sum(np.multiply(inp, mask)), np.sum(mask))
    st = np.sqrt(np.sum((np.multiply(mask, mea) + np.multiply(inp, mask) - mea) ** 2) / np.sum(mask))
    return np.multiply((inp - mea) / st, mask)


def find_sub_pixel_max_value(inp, k=4):
    """
    To find the highest point in a 2D array, with subpixel accuracy based on 1D spline interpolation .

    @param inp: A 2D numpy array containing the data points.
    @type inp: numpy array 2D
    @param k: The smoothing factor used in the spline interpolation, must be 1 <= k <= 5.
    @type k: int
    @return: A list of all points of maximal value in the structure of tuples with the x position, the y position and
        the value.
    @returntype: list
    """
    # assert len(inp.shape) == 2
    # assert isinstance(k, int) and 1 <= k <= 5

    import numpy as np
    from scipy.interpolate import InterpolatedUnivariateSpline

    v = np.amax(inp)  # the max value
    result = np.where(inp == v)  # arrays of x and y positions of max values
    output = []

    for xp, yp in zip(result[0], result[1]):
        # Find the highest point for x (first check if on sides otherwise interpolate)
        if xp == 1 or xp == inp.shape[0]:
            x = xp
            xv = v
        else:
            f = InterpolatedUnivariateSpline(range(0, inp.shape[0]), inp[:, yp], k=k)  # spline interpolation
            cr_pts = f.derivative().roots()
            cr_vals = f(cr_pts)
            val = np.argmax(cr_vals)
            x = cr_pts[val]
            xv = cr_vals[val]

        # Find the highest point for y (first check if on sides otherwise interpolate)
        if yp == 1 or yp == inp.shape[1]:
            y = yp
            yv = v
        else:
            f = InterpolatedUnivariateSpline(range(0, inp.shape[1]), inp[xp, :], k=k)  # spline interpolation
            cr_pts = f.derivative().roots()
            cr_vals = f(cr_pts)
            val = np.argmax(cr_vals)
            y = cr_pts[val]
            yv = cr_vals[val]

        # Calculate the average of the max value to return a single value which is maybe more close to the true value
        output.append((x, y, (xv + yv) / 2))

    return output


def find_sub_pixel_max_value_2d(inp, interpolate_factor=10, smoothing=2, dim=10, border_size=2, ignore_border=75):
    """
    To find the highest point in a given numpy array based on 2d spline interpolation, returns the maximum with subpixel
    precision.

    @param inp: The input data array (numpy 2D)
    @type inp: numpy array 2D
    @param interpolate_factor: The amount of interpolation to be done
    @type interpolate_factor: int
    @param smoothing: The amount of smoothing in the spline interpolation
    @type smoothing: int
    @param dim: The dimensions of the peak cutout, which is interpolated to find the subpixel maximum (initial pixels)
    @type dim: int
    @param border_size: The amount of pixels (initial pixels) to disregard in the peak cutout
    @type border_size: int
    @param ignore_border: The amount of pixels (initial pixels) to disregard in the initial finding of the initial maximum, to force the
       found maximum to be more in the center
    @type ignore_border: int
    @return: The subpixel maximum (x, y, the interpolated peak (excluding the border area))
    @returntype: tuple

     <-------- Vol_Size ------->
    | ignore_border             |
    | |   <------ a ------>     |
    | -> |                 | <- |
    |    | Here the max is |    |
    |    | found           |    |
    |    |    d> <c> <d    |    |
    |    |   |.. max ..|   |    |
    |    |   |... * ...|   |    |
    |    |   |.........|   |    |
    |    |    <-- b -->    |    |
    |    -------------------    |
    |___________________________|

    a: vol_size - 2 * ignore_border     (original pixels)
    b: dim * 2                          (original pixels)
    c: b * interpolate_factor - 2 * d                        (interpolated pixels)
    d: border_size                      (interpolated pixels)
    ...: interpolated values
    *: peak found

    """
    # assert len(inp.shape) == 2
    # assert isinstance(interpolate_factor, int) and interpolate_factor > 0
    # assert isinstance(smoothing, float) and smoothing >= 0
    # assert isinstance(dim, int) and dim > 0
    # assert isinstance(border_size, int) and border_size > 0
    # assert isinstance(ignore_border, int) and ignore_border > 0

    import numpy as np
    from scipy import interpolate
    import warnings

    border_size = border_size * interpolate_factor

    # Get the position of the initial maximum
    inp_without_border = inp[ignore_border:-ignore_border, ignore_border:-ignore_border]
    initial_max = np.unravel_index(inp_without_border.argmax(), inp_without_border.shape)
    # Reset the coordinates to be relative to the original inp(ut)
    initial_max = (initial_max[0] + ignore_border, initial_max[1] + ignore_border)

    # Get the starting points of the peak cutout
    x_dim = inp.shape[0]
    y_dim = inp.shape[1]
    x_start = max([0, initial_max[0] - dim])
    x_end = min([x_dim, initial_max[0] + dim])
    y_start = max([0, initial_max[1] - dim])
    y_end = min([y_dim, initial_max[1] + dim])

    # Create a grid to save the original points and one to save the interpolated points
    x, y, = np.mgrid[x_start:x_end:complex(x_end - x_start), y_start:y_end:complex(y_end - y_start)]
    xnew, ynew = np.mgrid[x_start:x_end:complex((x_end - x_start) * interpolate_factor),
                          y_start:y_end:complex((y_end - y_start) * interpolate_factor)]

    # Interpolate the points
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tck = interpolate.bisplrep(x, y, inp[x_start:x_end, y_start:y_end], s=smoothing)
        interpolated_grid = interpolate.bisplev(xnew[:, 0], ynew[0, :], tck)
        cropped_inter_grid = interpolated_grid[border_size:-border_size, border_size:-border_size]
        result = np.unravel_index(cropped_inter_grid.argmax(), cropped_inter_grid.shape)

        # Reset the coordinates to point to a place in the original data array
        result = (float(result[1])/interpolate_factor + y_start, float(result[0])/interpolate_factor + x_start)

        return result[0], result[1], cropped_inter_grid


def gettime():
    from time import gmtime, strftime
    return strftime("%H:%M:%S", gmtime())


def create_fsc_mask(fsc, size):
    from numpy import meshgrid, arange, sqrt, zeros_like, float32

    X, Y, Z = meshgrid(arange(size), arange(size), arange(size))

    X -= size // 2
    Y -= size // 2
    Z -= size // 2

    R = sqrt(X ** 2 + Y ** 2 + Z ** 2).astype(int)

    out = zeros_like(R).astype(float32)

    for n, val in enumerate(fsc):
        out[R == n] = val

    return out[size // 2, :, :]
