#!/usr/bin/env pytom

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
        print(offset)
        offset[0] = -sx / 2 + offset[0] * binning
        offset[1] = -sy / 2 + offset[1] * binning
        offset[2] = -sx / 2 + offset[2] * binning
        print(offset)
        from pytom.basic.structures import PickPosition
        for particle in particlelist:
            pickPosition = particle.getPickPosition()
            x = (pickPosition.getX() * binning + offset[0])
            y = (pickPosition.getY() * binning + offset[1])
            z = (pickPosition.getZ() * binning + offset[2])
            particle.setPickPosition(PickPosition(x=x, y=y, z=z))

        projectionslist.reconstructVolumes(particles=particlelist[:1], cubeSize=vol_size, \
                                       binning=1, applyWeighting=True, \
                                       showProgressBar=False, verbose=False, \
                                       preScale=1, postScale=1)

    elif projection_method == 'INFR':

        print(len(projectionslist))
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

            from pylab import *
            from matplotlib.colors import LogNorm
            fig, ax = subplots(1, 1, figsize=(5, 5))

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


def local_alignment(particleListName, particlelist, projections, projectionslist, vol_size, binning, offset, tilt_angles, iterations, projection_method = 'WBP'):
    from math import cos, sin, pi, ceil
    from pytom.tompy.transform import cut_from_projection
    from pytom.tools.ProgressBar import FixedProgBar
    from pytom.reconstruction.reconstructionStructures import ProjectionList
    from pytom.tompy.io import read
    import numpy as np
    import glob

    #write_subtomograms(particlelist, projectionslist, offset, binning, vol_size, iterations, tilt_angles, projection_method)
    subtomograms = [read(file) for file in glob.glob("/Subtomograms/{:s}/*.em".format(particlelist.getDirectory()))]

    print(len(projections))
    dim_x = projections[0].shape[0]
    dim_y = projections[0].shape[1]
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

        from pylab import *
        from matplotlib.colors import LogNorm
        fig, ax = subplots(1, 1, figsize=(5, 5))

        subregions = ProjectionList()
        n = 0
        for img, ang in zip(projections[:], tilt_angles[:]):

            n += 1
            # project the coordinate to 2D image
            yy = y  # assume the rotation axis is around y
            xx = (cos(ang * pi / 180) * (x - dim_x / 2) - sin(ang * pi / 180) * (z - dim_z / 2)) + dim_x / 2


            # cut the small patch out
            patch = cut_from_projection(img, [xx, yy], [vol_size, vol_size])
            patch = patch - np.mean(patch)

            # fill in the subregion
            subregions.append(patch)

            # visualise
            patch = patch[:, :, 0]
            # f = abs(fftshift(fftn(region)))**2
            # print region.shape
            ax.imshow(patch, cmap='gray')
            ax.text(0, -3, "x: {:f} y: {:f} ang: {:f}".format(xx, yy, ang))
            # ax[1].imshow(f, cmap='binary', norm=LogNorm() )
            savefig("image_particle_{:04d}_tiltimage_{:02d}.png".format(m, n))
            ax.clear()

            print("Finished {:02d}".format(n))

        # img = np.array(imread("/data2/dschulte/BachelorThesis/Data/VPP/05_Subtomogram_Analysis/Smallcctest_0.png"))
        # img_int = img[:, :, 0] / 3 + img[:, :, 1] / 3 + img[:, :, 2] / 3
        # img2 = np.array(imread("/data2/dschulte/BachelorThesis/Data/VPP/05_Subtomogram_Analysis/Smallcctest_1.png"))
        # img2_int = img2[:, :, 0] / 3 + img2[:, :, 1] / 3 + img2[:, :, 2] / 3
        # mask = np.array(imread("/data2/dschulte/BachelorThesis/Data/VPP/05_Subtomogram_Analysis/mask.png"))
        # mask_int = mask[:, :, 0] / 3 + mask[:, :, 1] / 3 + mask[:, :, 2] / 3
        #
        # autoccf = normalized_cross_correlation_numpy_mask(img_int, img2_int, mask_int)
        # print(autoccf)
        # fig, ax = subplots(1, 2, figsize=(10, 5))
        # ax[0].imshow(autoccf, cmap='gray')
        # d = 5
        # ax[1].imshow(autoccf[img.shape[0]/2-d:img.shape[0]/2+d,img.shape[1]/2-d:img.shape[1]/2+d], cmap='gray')
        # savefig("autocorrelation_0_0.png")




        fig, ax = subplots(2, 3, figsize=(15, 10))

        ax[1, 0].imshow(v[:,:,0])
        ax[1, 1].imshow(v[:, 0, :])
        ax[1, 2].imshow(v[0, :, :])
        ax[0, 0].imshow(v.sum(axis=0))
        ax[0, 1].imshow(v.sum(axis=1))
        ax[0, 2].imshow(v.sum(axis=2))
        ax[0, 0].text(0, -3, "Projections")
        ax[1, 0].text(0, -3, "Slices")

        savefig("projections_and_slices.png")

        print("Finished subtomogram assembly")

        template = v.sum(axis=2)

        fig, ax = subplots(1, 2, figsize=(10, 5))

        marker_style = dict(color='tab:blue', linestyle=':', marker='o', markersize=5, markerfacecoloralt='tab:red')
        mask = np.array(imread("/data2/dschulte/BachelorThesis/Scripts/FourrierTests/mask_softedges.png"))
        mask_int = mask[:, :, 0] / 3 + mask[:, :, 1] / 3 + mask[:, :, 2] / 3

        for n, img in enumerate(subregions):
            img = img.squeeze()
            print(img.shape, template.shape)
            ccf = normalized_cross_correlation_mask_numpy(template, img, mask_int)
            ccf = ccf / np.amax(ccf)
            #print(ccf, np.amax(ccf), np.amin(ccf))
            ax[0].imshow(ccf, cmap='gray')
            points = find_sub_pixel_max_value(ccf)
            ax[0].plot([p[1] for p in points], [p[0] for p in points], fillstyle='none', **marker_style)
            d = 10;
            ax[1].imshow(ccf[int(points[0][0])-d:int(points[0][0])+d, int(points[0][1])-d:int(points[0][1]+d)], cmap='gray')
            savefig("ccf_particle_{:04d}_tiltimage_{:02d}.png".format(m, n))
            ax[0].clear()
            ax[1].clear()


        # ySliceThroughCenter = object3d[:,dimy//2,:]
        # projectionZ = object3d.sum(axis=0)

        #print(len(subregions))
        # use the subregion array to try to locally optimize the alignment
        #for n, region in enumerate(subregions):
        #    region= region[:,:,0]
        #    #f = abs(fftshift(fftn(region)))**2
        #    print region.shape
        #    ax.imshow(region, cmap='gray')
        #    #ax[1].imshow(f, cmap='binary', norm=LogNorm() )
        #    savefig("image_particle_{:04d}_tiltimage_{:02d}.png".format(m, n) )
        #    ax.clear()
        # reconstruct
        #v = fourier_2d1d_iter_reconstruct(subregions, tilt_angles, iter)

        # write to the disk
        #write(p.getFilename(), v)


def normalized_cross_correlation_numpy(first, second):
    # tested, seems to work see /data2/dschulte/pytom-develop/Tests/
    import numpy.fft as nf
    import numpy as np

    prod = 1;
    for d in first.shape:
        prod *= d

    return np.real(nf.fftshift(nf.ifftn(np.multiply(nf.fftn(second), np.conj(nf.fftn(first)))))) / prod

def normalized_cross_correlation_mask_numpy(first, second, mask):
    # tested, seems to work see /data2/dschulte/pytom-develop/Tests/
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
    To find the highest point in a 2D array, with subpixel accuracy based on spline interpolation .

    :param inp: A 2D numpy array containing the data points.
    :param k: The smoothing factor used in the spline interpolation, must be 1 <= k <= 5.
    :return: A list of all points of maximal value in the structure of tuples with the x position, the y position and the value.
    """
    import numpy as np
    from scipy.interpolate import InterpolatedUnivariateSpline

    v = np.amax(inp)  # the max value
    result = np.where(inp == v)  # arrays of x and y positions of max values
    output = []

    for xp, yp in zip(result[0], result[1]):
        x = y = xv = yv = 0;

        # Find the highest point for x (first check if on sides otherwise interpolate)
        if xp == 1 or xp == inp.shape[0]:
            x = xp
            xv = v
        else:
            f = InterpolatedUnivariateSpline(range(0, inp.shape[0]), inp[:, yp], k=k);  # spline interpolation
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
            f = InterpolatedUnivariateSpline(range(0, inp.shape[1]), inp[xp, :], k=k);  # spline interpolation
            cr_pts = f.derivative().roots()
            cr_vals = f(cr_pts)
            val = np.argmax(cr_vals)
            y = cr_pts[val]
            yv = cr_vals[val]

        # Calculate the average of the max value to return a single value which is maybe more close to the true value
        output.append((x, y, (xv + yv) / 2))

    return output
