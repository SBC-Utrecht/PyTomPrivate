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

    :param particlelist: list with all particles (ParticleList())
    :param projection_directory: the directory of the projections (str)
    :param offset: the offset (x,y,z)
    :param binning: the binningfactor used (int)
    :param vol_size: the size of the volume to be reconstructed (in pixels, int)
    :param tilt_angles: the tilt angles as belonging to the projections (list(int))
    :param reconstruction_method: the method used to reconstruct the volumes (choose "WBP" or "INFR")
    :param infr_iterations: the amount of iterations used when choosing INFR (int)
    :return: nothing, writes the volumes to disk
    """
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
    else:
        raise Exception(reconstruction_method + " is not a known reconstruction method, use WBP or INFR")


def local_alignment(projections, vol_size, binning, offset, tilt_angles, particle_list_filename, projection_directory,
                    projection_method='WBP', infr_iterations=1, create_graphics=True, create_subtomograms=False,
                    averaged_subtomogram=False):
    """
    To polish a particle list based on (an) initial subtomogram(s).

    :param projections: a list with filenames of projections (list(Str))
    :param vol_size: the size of the volume to build the new subtomogram in (in pixels)
    :param binning: the binning factor used
    :param offset: the offset used (x, y, z)
    :param tilt_angles: the list of tiltangles used (list(Int))
    :param particle_list_filename: the filename of the particlelist
    :param projection_directory: the directory of the projections
    :param projection_method: the projection method used when create_subtomograms is True
                (possible values: WBP and INFR)
    :param infr_iterations: the amount of iterations to use when INFR projection is used
    :param create_graphics: to create plots of major parts of the algorithm, mainly used for debugging
                and initial creation
    :param create_subtomograms: to flag if initial subtomograms of the particles should be made (takes a long time!)
    :param averaged_subtomogram: to give a path to an averaged subtomogram to be used instead of subtomograms of all
                particles separately (False for off, otherwise a string)
    :return: nothing, it writes everything to disk
    """
    import numpy as np
    import glob
    from pytom.tompy.mpi import MPI

    mpi = MPI()
    mpi.begin()

    # load particle list
    from pytom.basic.structures import ParticleList

    particlelist = ParticleList()
    particlelist.fromXMLFile(particle_list_filename)

    particle_list_name = "particleList_TM_tomogram_010_WBP"

    # If needed create all subtomograms
    if create_subtomograms:
        write_subtomograms(particlelist, projection_directory, offset, binning, vol_size, tilt_angles,
                           projection_method, infr_iterations)

    # Read all subtomograms
    subtomograms = []
    if not averaged_subtomogram:
        subtomograms = [f for f in glob.glob("Subtomograms/{:s}/*".format(particle_list_name))]

    print("Creating the input array")

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

        # loop over tiltrange, take patch and cross correlate with reprojected subtomogram
        for img, ang in zip(projections, tilt_angles):
            input_to_processes.append([ang, subtomogram, offset, vol_size, particle.getPickPosition().toVector(),
                                       particle.getFilename(), particle_number, binning, img, create_graphics])

    print(len(input_to_processes))
    print("Created the input array\nStarting on running the processes")

    output = mpi.parfor(run_single_tilt_angle_unpack, input_to_processes)

    print("Ran the processes")

    results_file = "local_alignment_results.txt"

    from pytom.basic.datatypes import fmtLAR, headerLocalAlignmentResults, LOCAL_ALIGNMENT_RESULTS
    np.savetxt(results_file, np.array(output, dtype=LOCAL_ALIGNMENT_RESULTS), fmt=fmtLAR,
               header=headerLocalAlignmentResults)

    run_polished_subtomograms(particle_list_filename, projection_directory, results_file, binning, offset, vol_size)
    mpi.end()


def run_single_tilt_angle_unpack(inp):
    """
    unpack a list with arguments to "run_single_tilt_angle"
    :param inp: the arguments to "run_single_tilt_angle" in the same order in a single list (iterable)
    :return: the value from "run_single_tilt_angle"
    """
    return run_single_tilt_angle(inp[0], inp[1], inp[2], inp[3], inp[4], inp[5], inp[6], inp[7], inp[8])


def run_single_tilt_angle(ang, subtomogram, offset, vol_size, particle_position, particle_filename, particle_number,
                          binning, img, create_graphics=False):
    """
    To run a single tilt angle to allow for parallel computing

    :param ang: the tilt angle
    :param subtomogram: the filename of the subtomogram
    :param offset: the offset used (x,y,z)
    :param vol_size: the size of the volume to be reconstructed (in pixels)
    :param particle_position: the position of the particle in vector format,
                as given by particle.pickPosition().toVector()
    :param particle_filename: the filename of the particle, as given by particle.getfilename()
    :param particle_number: the number of the particle, to allow for unique mapping
    :param binning: the binning factor used
    :param img: the filename of the projection to be used
    :param create_graphics: to flag if images should be created for human inspection of the work done
    :return: the newly found positions of the particle, as a list  in the LOCAL_ALIGNMENT_RESULTS format
    """
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
    rotated = rotate3d(subtomogram, the=ang)  # 'the' is the rotational axis
    template = rotated.sum(axis=2)

    # Get coordinates of the paricle adjusted for the tilt angle
    yy = y  # assume the rotation axis is around y
    xx = (cos(ang * pi / 180) * (x - dim_x / 2) - sin(ang * pi / 180) * (z - dim_z / 2)) + dim_x / 2

    # cut the small patch out
    patch = cut_from_projection(img, [xx, yy], [vol_size, vol_size])
    patch = patch - np.mean(patch)
    patch = patch.squeeze()

    # Cross correlate the template and patch, this should give the pixel shift it is after
    ccf = normalized_cross_correlation_numpy(template, patch)
    points2d = find_sub_pixel_max_value_2d(ccf)

    if create_graphics:
        m_style = dict(color='tab:blue', linestyle=':', marker='o', markersize=5, markerfacecoloralt='tab:red')
        m_style_alt = dict(color='tab:red', linestyle=':', marker='o', markersize=5, markerfacecoloralt='tab:blue')
        v_m_style = dict(color='tab:blue', linestyle=':', marker='o', markersize=2, markerfacecoloralt='tab:red')
        v_m_style_alt = dict(color='tab:red', linestyle=':', marker='o', markersize=2, markerfacecoloralt='tab:blue')

        points = find_sub_pixel_max_value(ccf)

        nx, ny, nz = particle_position
        nx += (points2d[0] - vol_size / 2) / binning
        ny += (points2d[1] - vol_size / 2) / binning

        npatch = cut_from_projection(img, [xx + points2d[0] - vol_size / 2, yy + points2d[1] - vol_size / 2],
                                     [vol_size, vol_size])
        npatch = npatch - np.mean(npatch)

        nccf = normalized_cross_correlation_numpy(template, npatch.squeeze())
        npoints = find_sub_pixel_max_value(nccf)
        npoints2d = find_sub_pixel_max_value_2d(nccf)

        import pylab as pp

        vr = 7
        grid = pp.GridSpec(vr * 3, 3, wspace=0.4, hspace=0.3)

        ax_0_0 = pp.subplot(grid[0:vr - 1, 0])
        ax_0_1 = pp.subplot(grid[0:vr - 1, 1])
        ax_0_2 = pp.subplot(grid[0:vr - 1, 2])
        ax_1_0 = pp.subplot(grid[vr:2 * vr - 1, 0])
        ax_1_1 = pp.subplot(grid[vr:2 * vr - 1, 1])
        ax_1_2 = pp.subplot(grid[vr:2 * vr - 1, 2])
        ax_2_0 = pp.subplot(grid[2 * vr:3 * vr - 1, 0])
        ax_2_1 = pp.subplot(grid[2 * vr:3 * vr - 1, 1])
        ax_2_2 = []
        for n in range(2 * vr, 3 * vr):
            ax_2_2.append(pp.subplot(grid[n, 2]))

        ax_0_0.text(0, -3, "Cutout")
        ax_0_0.imshow(patch)
        ax_0_1.text(0, -3, "Template (projected subtomogram)")
        ax_0_1.imshow(template)
        ax_0_2.text(0, -3, "Shifted Cutout (based on cross correlation)")
        ax_0_2.imshow(npatch.squeeze())

        ax_1_0.text(0, -3, u"Cross correlation cutout × template")
        ax_1_0.imshow(ccf)
        ax_1_0.plot([p[1] for p in points], [p[0] for p in points], fillstyle='none', **m_style)
        ax_1_0.plot([points2d[0]], [points2d[1]], fillstyle='none', **m_style_alt)
        ax_1_0.plot([vol_size / 2], [vol_size / 2], ",k")

        ax_1_1.text(0, 0,
                    "Red: 2D spline interpolation, x: {:f} y: {:f}\nBlue: 1D spline interpolation, x: {:f} y: {:f}"
                    "\nBlack: center".format(
                        points2d[0] - vol_size / 2, points2d[1] - vol_size / 2, points[0][0] - vol_size / 2,
                        points[0][1] - vol_size / 2))
        ax_1_1.axis('off')

        ax_1_2.text(0, -3, u"Cross correlation shifted cutout × template")
        ax_1_2.imshow(nccf)
        ax_1_2.plot([p[1] for p in npoints], [p[0] for p in npoints], fillstyle='none', **m_style)
        ax_1_2.plot([npoints2d[0]], [npoints2d[1]], fillstyle='none', **m_style_alt)
        ax_1_2.plot([vol_size / 2], [vol_size / 2], ",k")

        ax_2_0.text(0, -3, u"Zoom into red peak in CC cutout × template")
        ax_2_0.get_xaxis().set_visible(False)
        ax_2_0.get_yaxis().set_visible(False)
        d = 10
        peak = ccf[int(points2d[0]) - d:int(points2d[0]) + d, int(points2d[1]) - d:int(points2d[1] + d)]
        ax_2_0.imshow(peak)

        ax_2_1.text(0, -3, u"Zoom into red peak in CC cutout × template, interpolated")
        ax_2_1.get_xaxis().set_visible(False)
        ax_2_1.get_yaxis().set_visible(False)
        ax_2_1.imshow(points2d[2])

        ax_2_2[0].text(0, -3, "Sections of the peak and interpolated figures")
        for n, ax in enumerate(ax_2_2):
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            stripe = (n + 0.5) / (vr + 1)
            peak_stripe = peak[:, int(stripe * peak.shape[1])]
            ax.plot(np.arange(0.0, 1.0, 1.0 / peak.shape[1]), peak_stripe / np.amax(peak_stripe),
                    **v_m_style)
            points2d_stripe = points2d[2][:, int(stripe * points2d[2].shape[1])]
            ax.plot(np.arange(0.0, 1.0, 1.0 / points2d[2].shape[1]), points2d_stripe / np.amax(points2d_stripe),
                    **v_m_style_alt)

        pp.savefig("polish_particle_{:04d}_tiltimage_{:02d}.png".format(particle_number, ang))

    # progressbar.increment_amount()
    return particle_number, points2d[0] - vol_size / 2, points2d[1] - vol_size / 2, ang, 0, 0, particle_filename


def run_polished_subtomograms(particle_list_filename, projection_directory, particle_polish_file, binning, offset, vol_size):
    """
    Reconstructs subtomograms based on a polished particlelist, writes these to the places as specified in particlelist

    :param particle_list_filename: The name of the file of the particlelist
    :param projection_directory: The name of the directory containing the projections
    :param particle_polish_file: The name of the file containing the polished alignment results
    :param binning: The binning factor
    :param offset: The reconstruction offset
    :param vol_size: The size of the particle
    :return: void
    """
    import os
    cwd = os.getcwd()

    batchfile = """#!/usr/bin/bash
#SBATCH --time        1:00:00
#SBATCH -N 1
#SBATCH --partition defq
#SBATCH --ntasks-per-node 1
#SBATCH --job-name    polishedReconstr                                                                       
#SBATCH --error="../LogFiles/%j-polished_subtomograms.err"
#SBATCH --output="../LogFiles/%j-polished_subtomograms.out"
#SBATCH --oversubscribe   

module load openmpi/2.1.1 python/2.7 lib64/append pytom/dev/dschulte

cd {:s}

reconstructWB.py --particleList {:s} \
--projectionDirectory {:s} \
--coordinateBinning {:d} \
--size {:d} \
--applyWeighting 0 \
--projBinning 1 \
--recOffset {:d},{:d},{:d} \
--particlePolishFile {:s}""".format(cwd, particle_list_filename, projection_directory, binning, vol_size, offset[0],
                                    offset[1], offset[2], particle_polish_file)

    f = open("polished_subtomograms.sh", "w+")
    f.write(batchfile)
    f.close()
    # os.system('sbatch polished_subtomograms.sh')


def normalized_cross_correlation_numpy(first, second):
    """
    Do a cross correlation based on numpy

    :param first: The first dataset (numpy 2D)
    :param second: The second dataset (numpy 2D)
    :return: The cross correlation result (shape == first.shape == second.shape)
    """
    import numpy.fft as nf
    import numpy as np

    prod = 1
    for d in first.shape:
        prod *= d

    return np.real(nf.fftshift(nf.ifftn(np.multiply(nf.fftn(second), np.conj(nf.fftn(first)))))) / prod


def normalized_cross_correlation_mask_numpy(first, second, mask):
    """
    Do cross correlation with a running mask based on numpy

    :param first: The first dataset (numpy 2D)
    :param second: The second dataset (numpy 2D)
    :param mask: The mask (numpy 2D) same size as first and second
    :return: The cross correlation result (shape == first.shape == second.shape)
    """
    import numpy.fft as nf
    import numpy as np

    a = norm_inside_mask(first, mask)
    b = norm_inside_mask(second, mask)

    return np.real(nf.fftshift(nf.ifftn(np.multiply(nf.fftn(b), np.conj(nf.fftn(a)))))) / np.sum(mask)


def norm_inside_mask(inp, mask):
    """
    To normalise a 2D array within a mask

    :param inp: A 2D array to be normalized.
    :param mask: A 2D array of the same size as the input to mask of parts of the inp.
    :return: A normalized 2D array.
    """
    import numpy as np

    mea = np.divide(np.sum(np.multiply(inp, mask)), np.sum(mask))
    st = np.sqrt(np.sum((np.multiply(mask, mea) + np.multiply(inp, mask) - mea) ** 2) / np.sum(mask))
    return np.multiply((inp - mea) / st, mask)


def find_sub_pixel_max_value(inp, k=4):
    """
    To find the highest point in a 2D array, with subpixel accuracy based on 1D spline interpolation .

    :param inp: A 2D numpy array containing the data points.
    :param k: The smoothing factor used in the spline interpolation, must be 1 <= k <= 5.
    :return: A list of all points of maximal value in the structure of tuples with the x position, the y position and
        the value.
    """
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


def find_sub_pixel_max_value_2d(inp, interpolate_factor=10, smoothing=2, dim=10, border_size=10, ignore_border=20):
    """
    To find the highest point in a given numpy array based on 2d spline interpolation, returns the maximum with subpixel
    precision.

    :param inp: The input data array (numpy 2D)
    :param interpolate_factor: The amount of interpolation to be done
    :param smoothing: The amount of smoothing in the spline interpolation
    :param dim: The dimensions of the peak cutout, which is interpolated to find the subpixel maximum
    :param border_size: The amount of pixels (after interpolation) to disregard in the peak cutout
    :param ignore_border: The amount of pixels to disregard in the initial finding of the initial maximum, to force the
        found maximum to be more in the center
    :return: The subpixel maximum (x, y, the interpolated peak (excluding the border area))
    """
    import numpy as np
    from scipy import interpolate
    import warnings

    # Get the position of the initial maximum
    inp_without_border = inp[ignore_border:-ignore_border, ignore_border:-ignore_border]
    initial_max = np.unravel_index(inp_without_border.argmax(), inp_without_border.shape)
    # Reset the coordinates to include the border
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
