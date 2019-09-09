#!/usr/bin/env pytom
# -*- coding: utf-8 -*-

"""
Created on Sep 02, 2019

@author: dschulte
"""


def write_subtomograms(particlelist, projectionslist, offset, binning, vol_size, iterations, tilt_angles, projection_method = 'WBP'):
    from math import cos, sin, pi, ceil
    from pytom.tompy.transform import cut_from_projection
    from nufft.reconstruction import fourier_2d1d_iter_reconstruct, fourier_2d1d_gridding_reconstruct
    from pytom.tools.ProgressBar import FixedProgBar
    from pytom.tompy.io import write
    import numpy as np

    if projection_method == 'WBP':
        # transform the cropping offset
        tmp = projectionslist[0]
        sx = tmp.getXSize()  # here should be the size of original projection!
        sy = tmp.getYSize()

        offset[0] = -sx / 2 + offset[0] * binning
        offset[1] = -sy / 2 + offset[1] * binning
        offset[2] = -sx / 2 + offset[2] * binning

        from pytom.basic.structures import PickPosition
        for particle in particlelist:
            pickPosition = particle.getPickPosition()
            x = (pickPosition.getX() * binning + offset[0])
            y = (pickPosition.getY() * binning + offset[1])
            z = (pickPosition.getZ() * binning + offset[2])
            particle.setPickPosition(PickPosition(x=x, y=y, z=z))

        # Reconstruct all particles
        projectionslist.reconstructVolumes(particles=particlelist, cubeSize=vol_size, \
                                       binning=1, applyWeighting=True, \
                                       showProgressBar=False, verbose=False, \
                                       preScale=1, postScale=1)

    elif projection_method == 'INFR':
        # NOT WORKING YET some bugs are still in

        print(len(projectionslist))
        print("INFR NOT TESTED YET")
        dim_x = projectionslist[0].getXSize()  # here should be the size of original projection!
        dim_y = projectionslist[0].getYSize()
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
            for img, ang in zip(projectionslist[:], tilt_angles[:]):
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
            v = fourier_2d1d_iter_reconstruct(subregions, tilt_angles, iterations)
            write("Subtomograms/{:s}/particle_{:d}.em".format(particlelist.getFilename(), m) , v)


def local_alignment(particlelist, projections, projectionslist, vol_size, binning, offset, tilt_angles, iterations, particlelistfilename, projection_directory, projection_method = 'WBP', create_graphics = True):
    from pytom.tompy.io import read
    import numpy as np
    import glob
    import multiprocessing as mp

    particlelistname = "particleList_TM_tomogram_010_WBP"

    # If needed create all subtomograms
    #write_subtomograms(particlelist, projectionslist, offset, binning, vol_size, iterations, tilt_angles, projection_method)

    # Read all subtomograms
    subtomograms = [read(file) for file in glob.glob("Subtomograms/{:s}/*".format(particlelistname))]

    # Get the size of the original projection
    dim_x = projectionslist[0].getXSize()
    dim_y = projectionslist[0].getYSize()
    dim_z = dim_x  # make z dim the same as x!

    output = []
    particlenumber = -1
    for particle, subtomogram in zip(particlelist, subtomograms):
        particlenumber += 1

        print("At particlenumber {:d}".format(particlenumber))

        if create_graphics:
            from pylab import savefig, subplots

            fig, ax = subplots(2, 3, figsize=(15, 10))

            ax[1, 0].imshow(subtomogram[:,:,0])
            ax[1, 1].imshow(subtomogram[:, 0, :])
            ax[1, 2].imshow(subtomogram[0, :, :])
            ax[0, 0].imshow(subtomogram.sum(axis=0))
            ax[0, 1].imshow(subtomogram.sum(axis=1))
            ax[0, 2].imshow(subtomogram.sum(axis=2))
            ax[0, 0].text(0, -3, "Projections")
            ax[1, 0].text(0, -3, "Slices")

            savefig("projections_and_slices.png")

            marker_style = dict(color='tab:blue', linestyle=':', marker='o', markersize=5, markerfacecoloralt='tab:red')
            marker_style_alt = dict(color='tab:red', linestyle=':', marker='o', markersize=5, markerfacecoloralt='tab:blue')
            vertical_marker_style = dict(color='tab:blue', linestyle=':', marker='o', markersize=2, markerfacecoloralt='tab:red')
            vertical_marker_style_alt = dict(color='tab:red', linestyle=':', marker='o', markersize=2, markerfacecoloralt='tab:blue')

        #from pylab import imread
        #mask = np.array(imread("/data2/dschulte/BachelorThesis/Scripts/FourrierTests/mask_softedges.png"))
        #mask_int = mask[:, :, 0] / 3 + mask[:, :, 1] / 3 + mask[:, :, 2] / 3

        x, y, z = particle.getPickPosition().toVector()
        x = (x + offset[0]) * binning
        y = (y + offset[1]) * binning
        z = (z + offset[2]) * binning

        # loop over tiltrange, take patch and crosscorrelate with reprojected subtomogram
        n = 0
        for img, ang in zip(projections, tilt_angles):
            print(ang)
            from pytom.tompy.transform import rotate3d

            n += 1

            # Get template
            rotated = rotate3d(subtomogram, the=ang) # 'the' is the rotational axis
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

            output.append((particlenumber, points2d[0] - vol_size / 2, points2d[1] - vol_size / 2, ang, 0, 0, particle.getFilename()))

            if create_graphics:
                points = find_sub_pixel_max_value(ccf)

                nx, ny, nz = particle.getPickPosition().toVector()
                nx += (points2d[0] - vol_size / 2) / binning
                ny += (points2d[1] - vol_size / 2) / binning

                npatch = cut_from_projection(img, [xx + points2d[0] - vol_size / 2, yy + points2d[1] - vol_size / 2], [vol_size, vol_size])
                npatch = npatch - np.mean(npatch)

                nccf = normalized_cross_correlation_numpy(template, npatch.squeeze())
                npoints = find_sub_pixel_max_value(nccf)
                npoints2d = find_sub_pixel_max_value_2d(nccf)

                import pylab as pp

                vr = 7;
                grid = pp.GridSpec(vr*3, 3, wspace=0.4, hspace=0.3)

                ax_0_0 = pp.subplot(grid[0:vr-1, 0])
                ax_0_1 = pp.subplot(grid[0:vr-1, 1])
                ax_0_2 = pp.subplot(grid[0:vr-1, 2])
                ax_1_0 = pp.subplot(grid[vr:2*vr-1, 0])
                ax_1_1 = pp.subplot(grid[vr:2*vr-1, 1])
                ax_1_2 = pp.subplot(grid[vr:2*vr-1, 2])
                ax_2_0 = pp.subplot(grid[2*vr:3*vr-1, 0])
                ax_2_1 = pp.subplot(grid[2*vr:3*vr-1, 1])
                ax_2_2 = []
                for n in range(2*vr, 3*vr):
                    ax_2_2.append(pp.subplot(grid[n, 2]))

                ax_0_0.text(0, -3, "Cutout")
                ax_0_0.imshow(patch)
                ax_0_1.text(0, -3, "Template (projected subtomogram)")
                ax_0_1.imshow(template)
                ax_0_2.text(0, -3, "Shifted Cutout (based on cross correlation)")
                ax_0_2.imshow(npatch.squeeze())

                ax_1_0.text(0, -3, u"Cross correlation cutout × template")
                ax_1_0.imshow(ccf)
                ax_1_0.plot([p[1] for p in points], [p[0] for p in points], fillstyle='none', **marker_style)
                ax_1_0.plot([points2d[0]], [points2d[1]], fillstyle='none', **marker_style_alt)
                ax_1_0.plot([vol_size / 2], [vol_size / 2], ",k")

                ax_1_1.text(0, 0, "Red: 2D spline interpolation, x: {:f} y: {:f}\nBlue: 1D spline interpolation, x: {:f} y: {:f}\nBlack: center".format(points2d[0] - vol_size / 2, points2d[1] - vol_size / 2, points[0][0] - vol_size / 2, points[0][1] - vol_size / 2))
                ax_1_1.axis('off')

                ax_1_2.text(0, -3, u"Cross correlation shifted cutout × template")
                ax_1_2.imshow(nccf)
                ax_1_2.plot([p[1] for p in npoints], [p[0] for p in npoints], fillstyle='none', **marker_style)
                ax_1_2.plot([npoints2d[0]], [npoints2d[1]], fillstyle='none', **marker_style_alt)
                ax_1_2.plot([vol_size / 2], [vol_size / 2], ",k")

                ax_2_0.text(0, -3, u"Zoom into red peak in CC cutout × template")
                ax_2_0.get_xaxis().set_visible(False)
                ax_2_0.get_yaxis().set_visible(False)
                d = 10
                peak = ccf[int(points2d[0])-d:int(points2d[0])+d, int(points2d[1])-d:int(points2d[1]+d)]
                ax_2_0.imshow(peak)

                ax_2_1.text(0, -3, u"Zoom into red peak in CC cutout × template, interpolated")
                ax_2_1.get_xaxis().set_visible(False)
                ax_2_1.get_yaxis().set_visible(False)
                ax_2_1.imshow(points2d[2])

                ax_2_2[0].text(0, -3, "Sections of the peak and interpolated figures")
                for n, ax in enumerate(ax_2_2):
                    ax.get_xaxis().set_visible(False)
                    ax.get_yaxis().set_visible(False)
                    stripe = (n+0.5)/(vr+1)
                    peak_stripe = peak[:,int(stripe*peak.shape[1])]
                    ax.plot(np.arange(0.0, 1.0, 1.0/peak.shape[1]), peak_stripe/np.amax(peak_stripe), **vertical_marker_style)
                    points2d_stripe = points2d[2][:,int(stripe*points2d[2].shape[1])]
                    ax.plot(np.arange(0.0, 1.0, 1.0/points2d[2].shape[1]), points2d_stripe/np.amax(points2d_stripe), **vertical_marker_style_alt)

                np.savetxt("data_{:04d}_tiltimage_{:02d}.txt".format(particlenumber, n), ccf)
                np.savetxt("data_interpolated_{:04d}_tiltimage_{:02d}.txt".format(particlenumber, n), points2d[2])
                pp.savefig("ccf_particle_{:04d}_tiltimage_{:02d}.png".format(particlenumber, n))

    from pytom.basic.datatypes import *
    np.savetxt("local_alignment_results.txt", np.array(output, dtype=LOCAL_ALIGNMENT_RESULTS), fmt=fmtLAR, header=headerLocalAlignmentResults)

    run_polished_subtomograms(particlelistfilename,projection_directory, "local_alignment_results.txt", binning, offset, vol_size)


def run_single_tilt_angle(ang, subtomogram, x, y, z, vol_size, dim_x, dim_z, particle, particlenumber, binning, img, create_graphics = False):
    print(ang)
    from pytom.tompy.transform import rotate3d
    import numpy as np
    from math import cos, sin, pi
    from pytom.tompy.transform import cut_from_projection

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
        marker_style = dict(color='tab:blue', linestyle=':', marker='o', markersize=5, markerfacecoloralt='tab:red')
        marker_style_alt = dict(color='tab:red', linestyle=':', marker='o', markersize=5, markerfacecoloralt='tab:blue')
        vertical_marker_style = dict(color='tab:blue', linestyle=':', marker='o', markersize=2, markerfacecoloralt='tab:red')
        vertical_marker_style_alt = dict(color='tab:red', linestyle=':', marker='o', markersize=2, markerfacecoloralt='tab:blue')

        points = find_sub_pixel_max_value(ccf)

        nx, ny, nz = particle.getPickPosition().toVector()
        nx += (points2d[0] - vol_size / 2) / binning
        ny += (points2d[1] - vol_size / 2) / binning

        npatch = cut_from_projection(img, [xx + points2d[0] - vol_size / 2, yy + points2d[1] - vol_size / 2],
                                     [vol_size, vol_size])
        npatch = npatch - np.mean(npatch)

        nccf = normalized_cross_correlation_numpy(template, npatch.squeeze())
        npoints = find_sub_pixel_max_value(nccf)
        npoints2d = find_sub_pixel_max_value_2d(nccf)

        import pylab as pp

        vr = 7;
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
        ax_1_0.plot([p[1] for p in points], [p[0] for p in points], fillstyle='none', **marker_style)
        ax_1_0.plot([points2d[0]], [points2d[1]], fillstyle='none', **marker_style_alt)
        ax_1_0.plot([vol_size / 2], [vol_size / 2], ",k")

        ax_1_1.text(0, 0,
                    "Red: 2D spline interpolation, x: {:f} y: {:f}\nBlue: 1D spline interpolation, x: {:f} y: {:f}\nBlack: center".format(
                        points2d[0] - vol_size / 2, points2d[1] - vol_size / 2, points[0][0] - vol_size / 2,
                        points[0][1] - vol_size / 2))
        ax_1_1.axis('off')

        ax_1_2.text(0, -3, u"Cross correlation shifted cutout × template")
        ax_1_2.imshow(nccf)
        ax_1_2.plot([p[1] for p in npoints], [p[0] for p in npoints], fillstyle='none', **marker_style)
        ax_1_2.plot([npoints2d[0]], [npoints2d[1]], fillstyle='none', **marker_style_alt)
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
                    **vertical_marker_style)
            points2d_stripe = points2d[2][:, int(stripe * points2d[2].shape[1])]
            ax.plot(np.arange(0.0, 1.0, 1.0 / points2d[2].shape[1]), points2d_stripe / np.amax(points2d_stripe),
                    **vertical_marker_style_alt)

        #np.savetxt("data_{:04d}_tiltimage_{:02d}.txt".format(particlenumber, n), ccf)
        #np.savetxt("data_interpolated_{:04d}_tiltimage_{:02d}.txt".format(particlenumber, n), points2d[2])
        #pp.savefig("ccf_particle_{:04d}_tiltimage_{:02d}.png".format(particlenumber, n))

    return particlenumber, points2d[0] - vol_size / 2, points2d[1] - vol_size / 2, ang, 0, 0, particle.getFilename()


def run_polished_subtomograms(particlelistfilename, projectiondirectory, particlepolishfile, binning, offset, vol_size):
    """
    Reconstructs subtomograms based on a polished particlelist, writes these to the places as specified in the particlelist

    :param particlelistfilename: The name of the file of the particlelist
    :param projectiondirectory: The name of the directory containing the projections
    :param particlepolishfile: The name of the file containing the polished alignment results
    :param binning: The binning factor
    :param offset: The reconstruction offset
    :param vol_size: The size of the particle
    :return: void
    """
    import os
    cwd = os.getcwd()

    file = """#!/usr/bin/bash
#SBATCH --time        1:00:00
#SBATCH -N 1
#SBATCH --partition fastq
#SBATCH --ntasks-per-node 1
#SBATCH --job-name    polishedReconstr                                                                       
#SBATCH --error="{:s}/polished_subtomograms-%j.err"
#SBATCH --output="{:s}/polished_subtomograms-%j.out"
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
--particlePolishFile {:s}""".format(cwd, cwd, cwd, cwd, particlelistfilename, projectiondirectory, binning, vol_size, offset[0], offset[1], offset[2], particlepolishfile)
    f = open("polished_subtomograms.sh", "w+")
    f.write(file)
    f.close()
    #os.system('sbatch polished_subtomograms.sh')


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
        x = y = xv = yv = 0

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
    x, y, = np.mgrid[x_start:x_end:complex(x_end-x_start), y_start:y_end:complex(y_end-y_start)]
    xnew, ynew = np.mgrid[x_start:x_end:complex((x_end-x_start)*interpolate_factor), y_start:y_end:complex((y_end-y_start)*interpolate_factor)]

    # Interpolate the points
    tck = interpolate.bisplrep(x, y, inp[x_start:x_end, y_start:y_end], s=smoothing)
    interpolated_grid = interpolate.bisplev(xnew[:,0], ynew[0,:], tck)
    cropped_inter_grid = interpolated_grid[border_size:-border_size, border_size:-border_size]
    result = np.unravel_index(cropped_inter_grid.argmax(), cropped_inter_grid.shape)

    # Reset the coordinates to point to a place in the original data array
    result = (float(result[1])/interpolate_factor + y_start, float(result[0])/interpolate_factor + x_start)

    return result[0], result[1], cropped_inter_grid