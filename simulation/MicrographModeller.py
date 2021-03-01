"""
Marten Chaillet's cryoET simulator - updated version of the simulator used in SHREC2019 and SHREC2020
Contributions from Ilja Gubins and Gijs van der Schot.
Original simulator (used in SHREC2019) written by Gijs van der Schot, which was loosely based on the simulator in the
TOM toolbox for matlab.
"""

# IO related modules
# from pytom.gui.mrcOperations import *
import configparser
import tracemalloc
import logging
import os
import datetime
import sys
from tqdm import tqdm

# math
# from pytom.basic.files import *
import numpy as xp
import random
import pytom.simulation.physics as physics

# Plotting, use Qt5Agg to prevent conflict with tkinter in pylab on cluster
# import matplotlib
# matplotlib.use('Qt5Agg')
# import matplotlib.pylab as plt


class ConfigLogger(object):
    """
    Facilitates writing the conf file to a .log file in the output_folder for reference of settings.
    """
    def __init__(self, log):
        self.__log = log

    def __call__(self, config):
        self.__log.info("Config:")
        config.write(self)

    def write(self, data):
        # stripping the data makes the output nicer and avoids empty lines
        line = data.strip()
        self.__log.info(line)


def downscale_class_mask(volume, binning, order=0):
    """
    Downscale class mask for binned reconstructions.
    """
    from scipy.ndimage import zoom

    if binning==1:
        return volume

    if volume.dtype != 'int':
        print('Dtype of class mask/occupancy mask was not int, forcing to int.')
        volume.astype(int)

    return zoom(volume, 1/binning, order=order)


def draw_range(range, datatype, name):
    """
    Input parsing from config file.
    """
    if type(range) == list and len(range) == 2:
        xp.random.seed(seed)
        random.seed(seed)
        if datatype == int:
            return xp.random.randint(range[0], range[1])
        elif datatype == float:
            return xp.random.uniform(range[0], range[1])
    elif type(range) == list and len(range) == 1:
        if datatype == int:
            return int(range[0])
        elif datatype == float:
            return float(range[0])
    elif type(range) == float or type(range) == int:
        if datatype == int:
            return int(range)
        elif datatype == float:
            return float(range)
    else:
        print(f'invalid data range or input type for parameter {name}')
        sys.exit(0)


def generate_model(particle_folder, save_path, listpdbs, listmembranes, pixel_size = 1E-10,
                   size=1024, thickness=200,
                   solvent_potential=physics.V_WATER, solvent_factor=1.0, number_of_particles=1000,
                   placement_size=512, retries=5000, number_of_markers=0,
                   absorption_contrast=False, voltage=300E3, number_of_membranes=0, sigma_motion=8):
    # IMPORTANT: We assume the particle models are in the desired voxel spacing for the pixel size of the simulation!
    from pytom.simulation.potential import create_gold_marker
    from pytom.voltools import transform
    from pytom.tompy.io import read_mrc, write
    from pytom.simulation.support import reduce_resolution

    # outputs
    X, Y, Z = size, size, thickness
    cell_real = xp.zeros((X, Y, Z))
    if absorption_contrast: cell_imag = xp.zeros_like(cell_real)

    # occupancy_bbox_mask = xp.zeros_like(cell_real)
    occupancy_accurate_mask = xp.zeros_like(cell_real)
    # class_bbox_mask = xp.zeros_like(cell_real)
    class_accurate_mask = xp.zeros_like(cell_real)
    ground_truth_txt_file = ''

    # load pdb volumes and pad them
    volumes_real = []
    if absorption_contrast: volumes_imag = []
    for pdb in listpdbs:
        try:
            # file paths for particles
            vol_base_name   = f'{pdb}_{pixel_size*1E10:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}V'
            filename_real       = os.path.join(particle_folder, f'{vol_base_name}_real.mrc')
            filename_imag       = os.path.join(particle_folder, f'{vol_base_name}_imag_{voltage*1E-3:.0f}V.mrc')
            # load the particle
            vol_real = read_mrc(filename_real)
            if absorption_contrast:
                vol_imag = read_mrc(filename_imag)
                # make sure real and imaginary part are the same size
                if vol_real.shape != vol_imag.shape:
                    print(f'real and imaginary interaction potential not the same shape for {pdb}, skipping model')
                    continue
                volumes_imag.append(vol_imag)
            volumes_real.append(vol_real)
        except Exception as ee:
            print(ee)
            raise Exception(f'Could not open pdb {pdb}, skipping the model')

    # attributes
    number_of_classes = len(listpdbs)
    names_of_classes = listpdbs
    particles_by_class = [0, ] * number_of_classes
    particle_nr = 1
    default_tries_left = retries
    skipped_particles = 0

    difference = size - placement_size
    if not difference:
        loc_x_start = loc_y_start = 0
        loc_x_end = loc_y_end = size
    else:
        loc_x_start = loc_y_start = difference // 2
        loc_x_end = loc_y_end = int(size - xp.ceil(difference / 2))

    # Add large cell structures, such as membranes first!
    if number_of_membranes:
        number_of_classes += 1
        names_of_classes.append('vesicles')
        particles_by_class += [0]

        # class id of cell structures is always the same
        cls_id = number_of_classes - 1

        for _ in tqdm(range(number_of_membranes), desc='Placing membranes and micelles'):

            # Give membranes a numbered names to randomly index them
            membrane_model = listmembranes[xp.random.randint(0, len(listmembranes))]

            # file names membranes
            vol_base_name   = f'{membrane_model}_{pixel_size*1E10:.2f}A_solvent-{solvent_potential*solvent_factor:.3f}V'
            filename_real       = os.path.join(particle_folder, f'{vol_base_name}_real.mrc')
            filename_imag       = os.path.join(particle_folder, f'{vol_base_name}_imag_{voltage*1E-3:.0f}V.mrc')
            # load the vesicle
            membrane_vol_real = read_mrc(filename_real)
            if absorption_contrast:
                membrane_vol_imag = read_mrc(filename_imag)
                if membrane_vol_real.shape != membrane_vol_imag.shape:
                    print(f'skipped mebrane model {membrane_model} due to real and imaginary shape not matching')
                    continue

            u = xp.random.uniform(0.0, 1.0, (2,))
            theta = xp.arccos(2 * u[0] - 1)
            phi = 2 * xp.pi * u[1]
            psi = xp.random.uniform(0.0, 2 * xp.pi)
            p_angles = xp.rad2deg([theta, phi, psi])

            # randomly mirror and rotate particle
            try:
                membrane_vol_real = transform(membrane_vol_real, rotation=p_angles,
                                        rotation_order='szxz', interpolation='linear', device='cpu')
                if absorption_contrast: membrane_vol_imag = transform(membrane_vol_imag, rotation=p_angles,
                                        rotation_order='szxz', interpolation='linear', device='cpu')
                # do something to remove edge artifacts of rotation? linear instead of filt_bspline
            except Exception as e:
                print(e)
                print('Something went wrong while rotating?')
                continue

            threshold = 0.001
            membrane_vol_real[membrane_vol_real < threshold] = 0
            if absorption_contrast: membrane_vol_imag[membrane_vol_imag < threshold] = 0

            # thresholded particle
            accurate_particle_occupancy = membrane_vol_real > 0

            # find random location for the particle
            # allow membranes to be placed half outside of the grandcell
            xx, yy, zz = membrane_vol_real.shape
            tries_left = default_tries_left
            # x_cut_left, x_cut_right = 0,0
            x_cut_left, x_cut_right, y_cut_left, y_cut_right, z_cut_low, z_cut_high = 0,0,0,0,0,0
            while tries_left > 0:
                # loc_x = xp.random.randint(loc_x_start + xx // 2 + 1, loc_x_end - xx // 2 - 1)
                loc_x = xp.random.randint(loc_x_start, loc_x_end)
                loc_y = xp.random.randint(loc_y_start, loc_y_end)
                loc_z = xp.random.randint(0, Z)

                tries_left -= 1

                x_cut_left = abs(loc_x - xx // 2) if (loc_x - xx // 2) < 0 else 0
                x_cut_right = abs(X - (loc_x + xx // 2 + xx % 2)) if (loc_x + xx // 2 + xx % 2) > Y else 0
                # adjust for bbox indexing for membrane sticking out of box
                # for x limit between x_start and x_end
                y_cut_left = abs(loc_y - yy//2) if (loc_y - yy//2) < 0 else 0
                y_cut_right = abs(Y - (loc_y + yy//2 + yy%2)) if (loc_y + yy//2 + yy%2) > Y else 0

                z_cut_low = abs(loc_z - zz//2) if (loc_z - zz//2) < 0 else 0
                z_cut_high = abs(Z - (loc_z + zz//2 + zz%2)) if (loc_z + zz//2 + zz%2) > Z else 0

                # crop the occupancy by removing overhang
                accurate_particle_occupancy_crop = accurate_particle_occupancy[x_cut_left:xx-x_cut_right,
                                                            y_cut_left:yy-y_cut_right, z_cut_low:zz-z_cut_high]

                # calculate coordinates of bbox for the newly rotated particle
                # bbox_x = [loc_x - xx // 2, loc_x + xx // 2 + xx % 2]
                bbox_x = [loc_x - xx // 2 + x_cut_left, loc_x + xx // 2 + xx % 2 - x_cut_right]
                bbox_y = [loc_y - yy // 2 + y_cut_left, loc_y + yy // 2 + yy % 2 - y_cut_right]
                bbox_z = [loc_z - zz // 2 + z_cut_low, loc_z + zz // 2 + zz % 2 - z_cut_high]

                # print(bbox_x, bbox_y, bbox_z)
                # print(xx,yy, zz)
                # print(y_cut_left, y_cut_right, z_cut_low, z_cut_high)

                # create masked occupancy mask
                masked_occupancy_mask = occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1],
                                        bbox_z[0]:bbox_z[1]]
                masked_occupancy_mask = masked_occupancy_mask * accurate_particle_occupancy_crop

                # if the location fits (masked occupancy pixel-wise mask is empty), break the loop and use this location
                if masked_occupancy_mask.sum() == 0:
                    break

            # however if still can't fit, ignore this particle (also adds variance in how many particles are
            # actually put)
            if tries_left < 1:
                skipped_particles += 1
                continue

            # populate occupancy volumes
            # occupancy_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = particle_nr
            occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                accurate_particle_occupancy_crop * particle_nr

            # populate class masks
            # class_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = (cls_id + 1)
            class_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                accurate_particle_occupancy_crop * (cls_id + 1)

            # populate density volume
            cell_real[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                membrane_vol_real[x_cut_left:xx-x_cut_right, y_cut_left:yy-y_cut_right, z_cut_low:zz-z_cut_high]
            if absorption_contrast: cell_imag[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                membrane_vol_imag[x_cut_left:xx-x_cut_right, y_cut_left:yy-y_cut_right, z_cut_low:zz-z_cut_high]

            # update stats
            particle_nr += 1
            particles_by_class[cls_id] += 1

            # update text
            ground_truth_txt_file += f'vesicle {int(loc_x - loc_x_start)} {int(loc_y - loc_y_start)} {int(loc_z)} ' \
                                     f'NaN NaN NaN\n'

    # GOLD MARKERS WILL ALSO BE COUNTED TOWARDS TOTAL PARTICLE NUMBER
    # There are also added as an additional class
    if number_of_markers:
        number_of_classes += 1
        names_of_classes.append('fiducials')
        particles_by_class += [0]

        # class id of gold markers is always the same
        cls_id = number_of_classes - 1

        for _ in tqdm(range(number_of_markers), desc='Placing gold markers'):

            # create the gold marker in other function
            if absorption_contrast:
                gold_real, gold_imag = create_gold_marker(pixel_size, solvent_potential, oversampling=2,
                                                            solvent_factor=solvent_factor,
                                                            imaginary=True, voltage=voltage)
            else:
                gold_real = create_gold_marker(pixel_size, solvent_potential, oversampling=2,
                                                 solvent_factor=solvent_factor)

            u = xp.random.uniform(0.0, 1.0, (2,))
            theta = xp.arccos(2 * u[0] - 1)
            phi = 2 * xp.pi * u[1]
            psi = xp.random.uniform(0.0, 2 * xp.pi)
            p_angles = xp.rad2deg([theta, phi, psi])

            # randomly mirror and rotate particle
            try:
                gold_real = transform(gold_real, rotation=p_angles,
                                        rotation_order='szxz', interpolation='linear', device='cpu')
                if absorption_contrast: gold_imag = transform(gold_imag, rotation=p_angles,
                                        rotation_order='szxz', interpolation='linear', device='cpu')
                # do something to remove edge artifacts of rotation? linear instead of filt_bspline
            except Exception as e:
                print(e)
                print('Something went wrong while rotating?')
                continue

            threshold = 0.001
            gold_real[gold_real < threshold] = 0
            if absorption_contrast: gold_imag[gold_imag < threshold] = 0

            # thresholded particle
            accurate_particle_occupancy = gold_real > 0

            # find random location for the particle
            xx, yy, zz = gold_real.shape
            tries_left = default_tries_left
            while tries_left > 0:
                loc_x = xp.random.randint(loc_x_start + xx // 2 + 1, loc_x_end - xx // 2 - 1)
                loc_y = xp.random.randint(loc_y_start + yy // 2 + 1, loc_y_end - yy // 2 - 1)
                loc_z = xp.random.randint(zz // 2 + 1, Z - zz // 2 - 1)

                tries_left -= 1

                # calculate coordinates of bbox for the newly rotated particle
                bbox_x = [loc_x - xx // 2, loc_x + xx // 2 + xx % 2]
                bbox_y = [loc_y - yy // 2, loc_y + yy // 2 + yy % 2]
                bbox_z = [loc_z - zz // 2, loc_z + zz // 2 + zz % 2]

                # create masked occupancy mask
                masked_occupancy_mask = occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1],
                                        bbox_z[0]:bbox_z[1]]
                masked_occupancy_mask = masked_occupancy_mask * accurate_particle_occupancy

                # if the location fits (masked occupancy pixel-wise mask is empty), break the loop and use this location
                if masked_occupancy_mask.sum() == 0:
                    break

            # however if still can't fit, ignore this particle (also adds variance in how many particles are
            # actually put)
            if tries_left < 1:
                skipped_particles += 1
                continue

            # populate occupancy volumes
            # occupancy_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = particle_nr
            occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                accurate_particle_occupancy * particle_nr

            # populate class masks
            # class_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = (cls_id + 1)
            class_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                accurate_particle_occupancy * (cls_id + 1)

            # populate density volume
            cell_real[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += gold_real
            if absorption_contrast: cell_imag[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += gold_imag

            # update stats
            particle_nr += 1
            particles_by_class[cls_id] += 1

            # update text
            ground_truth_txt_file += f'fiducial {int(loc_x - loc_x_start)} {int(loc_y - loc_y_start)} {int(loc_z)} ' \
                                     f'NaN NaN NaN\n'

    for _ in tqdm(range(number_of_particles), desc='Placing particles'):

        # select random class but correct for artifact class if adding gold particles

        cls_id = xp.random.randint(0, number_of_classes - bool(number_of_markers) - bool(number_of_membranes))

        # generate random rotation
        # to be properly random, use uniform sphere sampling
        # https://math.stackexchange.com/a/442423/72032
        # https://en.wikipedia.org/wiki/Rotation_matrix#Uniform_random_rotation_matrices
        # http://corysimon.github.io/articles/uniformdistn-on-sphere/
        u = xp.random.uniform(0.0, 1.0, (2,))
        theta = xp.arccos(2 * u[0] - 1)
        phi = 2 * xp.pi * u[1]
        psi = xp.random.uniform(0.0, 2 * xp.pi)
        p_angles = xp.rad2deg([theta, phi, psi])

        # randomly mirror and rotate particle
        try:
            particle_real = volumes_real[cls_id]
            if absorption_contrast: particle_imag = volumes_imag[cls_id]
            if xp.random.randint(2): # Generate true/false randomly
                # Mirror the particle to cover both left and right handedness of the proteins
                ax = xp.random.randint(3)
                particle_real = xp.flip(particle_real, axis=ax)
                if absorption_contrast: particle_imag = xp.flip(particle_imag, axis=ax)
            rotated_particle_real = transform(particle_real, rotation=p_angles,
                                         rotation_order='szxz', interpolation='filt_bspline', device='cpu')
            if absorption_contrast: rotated_particle_imag = transform(particle_imag, rotation=p_angles,
                                         rotation_order='szxz', interpolation='filt_bspline', device='cpu')
        except Exception as e:
            print(e)
            print('Something went wrong while rotating?')
            continue

        # remove particle rotation artifacts
        # threshold = min(volumes[cls_id][volumes[cls_id] > 0]) / 10
        # the approach above doesn't work well with PDB particles, there are often voxels with values of ^-11
        threshold = 0.01
        rotated_particle_real[rotated_particle_real < threshold] = 0
        if absorption_contrast: rotated_particle_imag[rotated_particle_imag < threshold]

        # thresholded rotated particle
        accurate_particle_occupancy = rotated_particle_real > 0

        # find random location for the particle
        xx, yy, zz = rotated_particle_real.shape
        tries_left = default_tries_left
        while tries_left > 0:
            loc_x = xp.random.randint(loc_x_start + xx // 2 + 1, loc_x_end - xx // 2 - 1)
            loc_y = xp.random.randint(loc_y_start + yy // 2 + 1, loc_y_end - yy // 2 - 1)
            loc_z = xp.random.randint(zz // 2 + 1, Z - zz // 2 - 1)

            tries_left -= 1

            # calculate coordinates of bbox for the newly rotated particle
            # bbox_x = [loc_x - dims[cls_id][0] // 2, loc_x + dims[cls_id][0] // 2 + dims[cls_id][0] % 2]
            # bbox_y = [loc_y - dims[cls_id][1] // 2, loc_y + dims[cls_id][1] // 2 + dims[cls_id][1] % 2]
            # bbox_z = [loc_z - dims[cls_id][2] // 2, loc_z + dims[cls_id][2] // 2 + dims[cls_id][2] % 2]
            bbox_x = [loc_x - xx // 2, loc_x + xx // 2 + xx % 2]
            bbox_y = [loc_y - yy // 2, loc_y + yy // 2 + yy % 2]
            bbox_z = [loc_z - zz // 2, loc_z + zz // 2 + zz % 2]

            # create masked occupancy mask
            masked_occupancy_mask = occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]]
            masked_occupancy_mask = masked_occupancy_mask * accurate_particle_occupancy

            # if the location fits (masked occupancy pixel-wise mask is empty), break the loop and use this location
            if masked_occupancy_mask.sum() == 0:
                break

        # however if still can't fit, ignore this particle (also adds variance in how many particles are actually put)
        if tries_left < 1:
            skipped_particles += 1
            continue

        # populate occupancy volumes
        # occupancy_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = particle_nr
        occupancy_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
            accurate_particle_occupancy * particle_nr

        # populate class masks
        # class_bbox_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] = (cls_id + 1)
        class_accurate_mask[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
            accurate_particle_occupancy * (cls_id + 1)

        # populate density volume
        cell_real[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += rotated_particle_real
        if absorption_contrast: cell_imag[bbox_x[0]:bbox_x[1], bbox_y[0]:bbox_y[1], bbox_z[0]:bbox_z[1]] += \
                                                    rotated_particle_imag

        # update stats
        particle_nr += 1
        particles_by_class[cls_id] += 1

        # update text
        ground_truth_txt_file += f'{listpdbs[cls_id]} {int(loc_x - loc_x_start)} {int(loc_y - loc_y_start)} {int(loc_z)} ' \
                                 f'{p_angles[0]:.4f} {p_angles[1]:.4f} {p_angles[2]:.4f}\n'

    # # Add structural noise
    # print('Adding structural noise to grand model cell')
    # noisy_cell = cell #+ xp.random.normal(0, sigma_structural, cell.shape)
    # pytom.tompy.io.write(f'{save_path}/grandmodel_noisefree_original.mrc', cell)

    # motion blur is simply a gaussian in fourier space to the required frequency component
    cell_real = reduce_resolution(cell_real, pixel_size * 1E10, sigma_motion)
    if absorption_contrast: cell_imag = reduce_resolution(cell_imag, pixel_size * 1E10, sigma_motion)

    # grandmodel names
    filename_gm_real    = os.path.join(save_path, 'grandmodel.mrc')
    filename_gm_imag    = os.path.join(save_path, 'grandmodel_imag.mrc')
    filename_cm         = os.path.join(save_path, 'class_mask.mrc')
    filename_om         = os.path.join(save_path, 'occupancy_mask.mrc')
    filename_loc        = os.path.join(save_path, 'particle_locations.txt')
    filename_con        = os.path.join(save_path, 'class_conversion_table.txt')

    # save grandmodels
    print('Saving grandmodels')
    write(filename_gm_real, cell_real)
    if absorption_contrast: write(filename_gm_imag, cell_imag)
    # save class masks
    print('Saving class volumes')
    # pytom.tompy.io.write(f'{save_path}/class_bbox.mrc', class_bbox_mask)
    write(filename_cm, class_accurate_mask)
    # save occupancy masks
    print('Saving occupancy volumes')
    # pytom.tompy.io.write(f'{save_path}/occupancy_bbox.mrc', occupancy_bbox_mask)
    write(filename_om, occupancy_accurate_mask)
    # save particle text file
    with open(filename_loc, 'w') as f:
        f.write(ground_truth_txt_file)
    # save conversion table from pdb to class
    with open(filename_con, 'w') as f:
        conversion_table = 'background 0\n'
        for i in range(number_of_classes):
            conversion_table += f'{names_of_classes[i]} {i+1}\n'
        f.write(conversion_table)

    # reporting
    print(f'Total number of particles in the tomogram: {particle_nr - 1}\n'
          f'Skipped {skipped_particles} particles ({default_tries_left} random location retries)\n'
          f'Particles by class: {particles_by_class}')

    # # debug: inspect all generated volumes
    # import napari
    # with napari.gui_qt():
    #     viewer = napari.Viewer()
    #     viewer.add_image(occupancy_bbox_mask, name='occupancy bbox mask', interpolation='bicubic')
    #     viewer.add_image(occupancy_accurate_mask, name='occupancy accurate mask', interpolation='bicubic')
    #     viewer.add_image(class_bbox_mask, name='class bbox mask', interpolation='bicubic')
    #     viewer.add_image(class_accurate_mask, name='class accurate mask', interpolation='bicubic')
    #     viewer.add_image(cell_real, name='cell_real', interpolation='bicubic')

    # trying to reduce intermediate memory usage
    del cell_real, class_accurate_mask, occupancy_accurate_mask  # , class_bbox_mask, occupancy_bbox_mask
    if absorption_contrast: del cell_imag
    return


def create_ice_layer(shape, angle, width, value=1.0, sigma=0.0):
    """
    Efficiently create an ice layer at specified angle within the volume shape.
    Assumes the angle rotates perpendicular to the y-axis of the volume.

    @param volume: volume to add the ice to
    @param angle: rotation angle of ice layer
    @param width: width of the ice layer in number of pixels
    @param value: value of ice layer
    @param sigma: std of gaussian filter for smoothing edges

    @author: Marten Chaillet
    """
    from scipy.ndimage import gaussian_filter
    assert xp.abs(angle) <= 90, print('rotation angle of ice layer cannot be larger than +- 90 degrees.')

    xsize = shape[0]
    ysize = shape[1]
    zsize = shape[2]

    x = xp.arange(-xsize / 2, xsize / 2, 1, dtype=xp.float32)
    z = xp.arange(-zsize / 2, zsize / 2, 1, dtype=xp.float32)
    zm = xp.tile(z[xp.newaxis, :], [xsize, 1])

    xline = x * xp.tan(-angle * xp.pi / 180)
    nwidth = width / xp.cos(angle * xp.pi / 180)
    xmin = xline - nwidth / 2
    xmax = xline + nwidth / 2

    square_min = xp.tile(xmin[:, xp.newaxis], [1, zsize])
    square_max = xp.tile(xmax[:, xp.newaxis], [1, zsize])

    # Create smooth edge for side 1 of the layer
    grad1 = zm - square_min
    range1 = xp.abs(grad1.min() - grad1.max())
    c1 = (range1 / grad1.shape[0]) / 0.5
    grad1[grad1 > c1] = c1
    grad1[grad1 < -c1] = -c1
    grad1 = (grad1 - grad1.min()) / (grad1.max() - grad1.min())
    # Create smooth edge for side 2 of the layer
    grad2 = square_max - zm
    range2 = xp.abs(grad2.min() - grad2.max())
    c2 = (range2 / grad2.shape[0]) / 0.5
    grad2[grad2 > c2] = c2
    grad2[grad2 < -c2] = -c2
    grad2 = (grad2 - grad2.min()) / (grad2.max() - grad2.min())

    # if a sigma is provided apply gaussian filter
    if not (sigma == 0.0):
        layer = gaussian_filter(grad1 * grad2 * value, sigma=sigma)
    else:
        layer = grad1 * grad2 * value

    # before returning tile the 2d layer to a volume
    return xp.tile(layer[:, xp.newaxis, :], [1, ysize, 1])


def microscope_single_projection(noisefree_projection, dqe, mtf, dose, pixel_size, binning=1):
    """
    Inspired by InSilicoTEM (Vulovic et al., 2013)
    @author: Marten Chaillet
    """
    # from scipy.ndimage import shift

    ntf = xp.sqrt(mtf ** 2 / dqe) # square root because dqe = mtf^2 / ntf^2
    ntf = xp.maximum(ntf, 1E-7) # ensure we do not divide by zero
    mtf_shift = xp.fft.ifftshift(mtf)
    ntf_shift = xp.fft.ifftshift(ntf)

    # NUMBER OF ELECTRONS PER PIXEL
    dose_per_pixel = dose * (pixel_size*1E10)**2 / binning**2 # from square A to square nm (10A pixels)
    print(f'Number of electrons per pixel (before binning and sample absorption): {dose_per_pixel}')

    # Fourier transform and multiply with sqrt(dqe) = mtf/ntf
    projection_fourier = xp.fft.fftn(xp.fft.ifftshift(noisefree_projection)) * mtf_shift / ntf_shift
    # projection_fourier = projection_fourier
    # Convert back to real space
    projection = xp.real(xp.fft.fftshift(xp.fft.ifftn(projection_fourier)))
    projection[projection<0] = 0
    # Draw from poisson distribution and scale by camera's conversion factor
    # conversion_factor = 100  # in ADU/e- , this is an arbitrary unit. Value taken from Vulovic et al., 2010

    # Apply poissonian noise
    # poisson_mean = projection * dose_per_pixel
    # sigma_motion = 1 / (pixel_size*1E9)
    projection_poisson = xp.zeros(projection.shape)
    for _ in range(binning**2):
        # generate shifts per projection to add motion blur
        # translation = (xp.random.normal(0, sigma_motion), xp.random.normal(0, sigma_motion))
        # poisson_intermediate = xp.random.poisson(lam=shift(projection, translation, mode='nearest') * dose_per_pixel)
        poisson_intermediate = xp.random.poisson(lam=projection * dose_per_pixel)
        projection_poisson += xp.real(xp.fft.fftshift(xp.fft.ifftn(xp.fft.fftn(xp.fft.ifftshift(poisson_intermediate))
                                                                   * ntf_shift))) / binning**2
        # imshow(projection_poisson)
        # show()

    # Image values are now in ADU
    # Apply the camera's noise transfer function to the noisy image
    # projection_fourier = xp.fft.fftn(xp.fft.ifftshift(projection_poisson)) * ntf_shift

    # Add additional noise from digitization process, less relevant for modern day cameras.
    # readout noise standard deviation can be 7 ADUs, from Vulovic et al., 2010
    # sigma_readout = 7
    # readsim = xp.random.normal(0, sigma_readout, projection.shape) # readout noise has a gaussian distribution
    # darksim = 0     # dark current noise has a poisson distribution, usually an order of magnitude smaller than
    #                 # readout noise and can hence be neglected

    # Add readout noise and dark noise in real space
    # return xp.real(xp.fft.fftshift(xp.fft.ifftn(projection_fourier))) # + readsim + darksim
    return projection_poisson


def parallel_project(grandcell, frame, image_size, pixel_size, msdz, n_slices, ctf, dose, dqe, mtf, voltage,
                     binning=1, translation=(.0,.0,.0), rotation=(.0,.0,.0), solvent_potential=physics.V_WATER,
                     solvent_absorption=.0, sigma_damage=.0, ice_voxels=None):
    """
    Only use multislice for this.
    @param model:
    @param frame:
    @param translation:
    @param rotation:
    @return:
    """
    from pytom.voltools import transform
    from pytom.simulation.microscope import transmission_function, fresnel_propagator

    print('Transforming sample for tilt/frame ', frame)

    sample = grandcell.copy()

    max_tilt_radians = abs(rotation[1]) * xp.pi / 180
    max_tilt_radians_opp = (90 - abs(rotation[1])) * xp.pi / 180
    rotation_height = int(xp.ceil(xp.sin(max_tilt_radians) * image_size +
                                           xp.sin(max_tilt_radians_opp) * ice_voxels))
    print(f'Reduced rotation height for relevant specimens: {rotation_height}')
    if rotation_height % 2: rotation_height += 1
    diff = sample.shape[2]-rotation_height
    i = diff // 2

    # model parameter is a complex volume or real volume depending on the addition of absorption contrast
    # first transform the volume
    if sample.dtype == 'complex64':
        transform(sample[...,i:-i].real, translation=translation, rotation=rotation, rotation_order='sxyz',
                                interpolation='filt_bspline', device='cpu', output=sample[...,i:-i].real)
        transform(sample[...,i:-i].imag, translation=translation, rotation=rotation, rotation_order='sxyz',
                                interpolation='filt_bspline', device='cpu', output=sample[...,i:-i].imag)
        # remove rotation artifacts
        threshold = 0.001
        sample.real[sample.real < threshold] = 0
        sample.imag[sample.imag < threshold] = 0
    elif sample.dtype == 'float32':
        transform(sample[...,i:-i], translation=translation, rotation=rotation, rotation_order='sxyz',
                                interpolation='filt_bspline', device='cpu', output=sample[...,i:-i])
        # remove rotation artifacts
        threshold = 0.001
        sample[sample < threshold] = 0
    else:
        print('Invalid dtype of sample.')
        sys.exit(0)

    print('Simulating projection with multislice method for frame/tilt ', frame)
    # Start with projection, first prepare parameters
    box_size = sample.shape[0]
    box_height = sample.shape[2]

    # add the ice layer to the sample
    if ice_voxels is None: # this does not seem like the best condition 'or not sample.dtype=='complex64' '
        ice_layer = 1
    else:
        ice_layer = create_ice_layer(sample.shape, rotation[1], ice_voxels, value=1.0, sigma=0.0)
        # ice layer datatype at this point should be np.float32
        # print(f'data type of ice: {ice_layer.dtype}')
    # apply structural deterioration due to beam damage via random noise with increment based on frame number
    # incremental damage based on frame number
    # if not (sigma_damage == 0.0):
    #     # this does not add noise to imaginary part, but maybe not needed
    #     sample += ( ice_layer * xp.random.normal(0, sigma_damage, sample.shape) )

    if sample.dtype == 'complex64':
        ice_layer *= solvent_potential
        sample.real += ice_layer # Maybe we can remove this? Are we only interested in absorption of ice layer?
        ice_layer *= (solvent_absorption / solvent_potential)
        sample.imag += ice_layer
    else:
        ice_layer *= solvent_potential
        sample += ice_layer

    del ice_layer # free memory

    if n_slices==box_height and image_size==box_size:
        projected_potent_ms = sample
        px_per_slice = 1
        num_px_last_slice = box_height % px_per_slice
    else:
        ileft = (box_size - image_size) // 2
        iright = -int(xp.ceil((box_size - image_size) / 2))

        # Allocate space for multislice projection
        projected_potent_ms = xp.zeros((image_size, image_size, n_slices), dtype=xp.complex64)

        px_per_slice = int(msdz/pixel_size)
        num_px_last_slice = box_height % px_per_slice  # consider last slice might have a different number of pixels

        # Project potential for each slice (phase grating)
        for ii in range(n_slices):
            if ileft==0 and iright==0:
                projected_potent_ms[:, :, ii] = sample[:,:, ii * px_per_slice: (ii + 1) * px_per_slice].mean(axis=2) #.get()
            else:
                projected_potent_ms[:, :, ii] = sample[ileft:iright, ileft:iright,
                                                ii * px_per_slice: (ii + 1) * px_per_slice].mean(axis=2) #.get()
    # at this point we no longer need the sample, because all information in now contained in the projected slices
    # free the memory to accomodate space
    del sample

    # calculate the transmission function for each slice
    psi_t = transmission_function(projected_potent_ms, voltage, msdz)

    # calculate the fresnel propagator (identical for same dz)
    propagator = fresnel_propagator(image_size, pixel_size, voltage, msdz)

    # Wave propagation with MULTISLICE method, psi_multislice is complex
    psi_multislice = xp.zeros((image_size,image_size), dtype=xp.complex64) + 1 # +1 for initial probability

    # Loop over all the slices, except the last one if it has a different slice size
    for ii in range(n_slices-min(1,num_px_last_slice)):
        wave_field = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, ii]) )
        psi_multislice = xp.fft.fftshift( xp.fft.ifftn((wave_field * xp.fft.ifftshift(propagator) )) )

    # Calculate propagation through last slice in case the last slice contains a different number of pixels
    if num_px_last_slice:
        msdz_end = num_px_last_slice * pixel_size
        psi_t[:, :, -1] = transmission_function(projected_potent_ms[:, :, -1], voltage, msdz_end)
        propagator_end = fresnel_propagator(image_size, pixel_size, voltage, msdz_end)
        wave_field = xp.fft.fftn( xp.fft.ifftshift(psi_multislice) * xp.fft.ifftshift(psi_t[:, :, -1]) )
        psi_multislice = xp.fft.fftshift( xp.fft.ifftn( wave_field * xp.fft.ifftshift(propagator_end) ) )

    # Multiple by CTF for microscope effects on electron wave
    wave_ctf = xp.fft.ifftshift(ctf) * xp.fft.fftn(xp.fft.ifftshift(psi_multislice) )
    # Intensity in image plane is obtained by taking the absolute square of the wave function
    noisefree_projection = xp.abs(xp.fft.fftshift(xp.fft.ifftn(wave_ctf))) ** 2

    # DEBUGGING
    # test = xp.log(xp.abs(xp.fft.fftshift(xp.fft.fftn(noisefree_projection))))
    # r1, m1 = radial_average(test)
    # fig, ax = plt.subplots()
    # ax.plot(r1, m1)
    # ax.legend()
    # plt.savefig(f'{folder}/radial.png')

    # Apply the microscope function
    projection = microscope_single_projection(noisefree_projection, dqe, mtf, dose, pixel_size, binning=binning)
    # Write the projection to the projection folder
    # pytom.tompy.io.write(f'{folder}/synthetic_{frame+1}.mrc', projection)
    # Return noisefree and projection as tuple for writing as mrc stack in higher function
    return (noisefree_projection, projection)


def generate_tilt_series_cpu(save_path, angles, nodes=1, image_size=None, rotation_box_height=None,
                             pixel_size=1E-9, binning=1, dose=80, voltage=300E3, spherical_aberration=2.7E-3,
                             chromatic_aberration=2.7E-3, energy_spread=0.7, illumination_aperture=0.030E-3,
                             objective_diameter=100E-6, focus_length=4.7E-3, astigmatism=0.0E-9, astigmatism_angle=0,
                             msdz=1E-9, defocus=2E-6, sigma_damage=0.0, camera_type='K2SUMMIT', camera_folder='',
                             solvent_potential=physics.V_WATER, absorption_contrast=False):
    from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS as dar
    from pytom.basic.datatypes import fmtAlignmentResults, HEADER_ALIGNMENT_RESULTS
    from pytom.gui.guiFunctions import savestar
    from pytom.simulation.microscope import create_detector_response, create_complex_ctf
    from pytom.tompy.io import read_mrc, write
    from joblib import Parallel, delayed
    # NOTE; Parameter defocus specifies the defocus at the bottom of the model!

    # grab model
    filename_gm_real = os.path.join(save_path, 'grandmodel.mrc')
    filename_gm_imag = os.path.join(save_path, 'grandmodel_imag.mrc')
    grandcell = read_mrc(filename_gm_real)
    if absorption_contrast:
        xp_type = xp.complex64
        grandcell = grandcell.astype(xp_type)
        grandcell.imag = read_mrc(filename_gm_imag)
        # calculate the absorption for amorphous ice at the specified voltage
        solvent_amplitude = physics.potential_amplitude(physics.AMORPHOUS_ICE_DENSITY, physics.WATER_MW, voltage)
        print(f'solvent absorption = {solvent_amplitude:.3f}')
    else:
        xp_type = xp.float32
        grandcell = grandcell.astype(xp_type)
        solvent_amplitude = 0.0

    box_size = grandcell.shape[0]
    box_height = grandcell.shape[2]

    # For sample set the arrays specifically to np.complex64 datatype to save memory space
    if rotation_box_height is not None:
        # Use the specified rotation box height
        rotation_volume = xp.zeros((box_size, box_size, rotation_box_height), dtype=xp_type)  # complex or real
        offset = (rotation_box_height - box_height) // 2

        if (rotation_box_height - box_height) % 2:
            rotation_volume[:, :, offset + 1:-offset] = grandcell[:, :, :]
        else:
            rotation_volume[:, :, offset:-offset] = grandcell[:, :, :]
    else:
        # Calculate maximal box height to contain all the information for the rotation
        max_tilt = max([abs(a) for a in angles])
        max_tilt_radians = max_tilt * xp.pi / 180
        rotation_box_height = int(xp.ceil(xp.tan(max_tilt_radians) * image_size +
                                          box_height / xp.cos(max_tilt_radians)))
        if rotation_box_height % 2: rotation_box_height += 1
        # Place grandcell in rotation volume
        rotation_volume = xp.zeros((box_size, box_size, rotation_box_height), dtype=xp_type)
        offset = (rotation_box_height - box_height) // 2

        if (rotation_box_height - box_height) % 2:
            rotation_volume[:, :, offset + 1:-offset] = grandcell[:, :, :]
        else:
            rotation_volume[:, :, offset:-offset] = grandcell[:, :, :]
    del grandcell

    if image_size is None:
        image_size = box_size

    # confirm image_size is valid
    assert image_size <= box_size, 'Specified projection image size is invalid as it is larger than the model dimension.'

    # adjust defocus because the value specifies defocus at the bottom of the box. the input expects defocus at the
    # center of the sample, therefore subtract half the box size.
    zheight = rotation_box_height * pixel_size  # thickness of rotation volume in nm
    defocus -= (zheight / 2)  # center defocus value at tilt angle

    # Check if msdz is viable, else correct it
    if msdz % pixel_size != 0:
        # Round the msdz to the closest integer mutltiple of the pixel size
        round_up = xp.round(msdz % pixel_size / pixel_size)
        if round_up:
            msdz += (pixel_size - msdz % pixel_size)
        else:
            msdz -= (msdz % pixel_size)
        print(f'We adjusted msdz to {msdz*1E9:.3f} nm to make it an integer multiple of pixel size.')

    # Determine the number of slices
    if msdz > zheight:
        n_slices = 1
        msdz = zheight
        print('The multislice step can not be larger than volume thickness. Use only one slice.\n')
    elif msdz < pixel_size:
        n_slices = box_height
        msdz = pixel_size
        print(
            'The multislice steps can not be smaller than the pixel size. Use the pixel size instead for the step size.\n')
    else:
        n_slices = int(xp.ceil(xp.around(zheight / msdz, 3)))
    print('Number of slices for multislice: ', n_slices)

    # determine dose per frame
    dose_per_tilt = dose / len(angles)

    # create motion blur by introducing random  x,y translation for each tilt
    sigma_motion = 0.0 # nm
    translations = []
    for i in angles:
        x = xp.random.normal(0, sigma_motion) / (pixel_size*1E9)
        y = xp.random.normal(0, sigma_motion) / (pixel_size*1E9)
        # print((x,y,0.0))
        translations.append((x,y,0.0))

    # defocus_series = [xp.random.normal(defocus, 0.2E-6) for a in angles]
    ctf_series = []
    for i in angles:
        dz = xp.random.normal(defocus, 0.2E-6)
        ctf = create_complex_ctf((image_size, image_size), pixel_size, dz, voltage=voltage,
                                     Cs=spherical_aberration, Cc=chromatic_aberration, energy_spread=energy_spread,
                                     illumination_aperture=illumination_aperture, objective_diameter=objective_diameter,
                                     focus_length=focus_length, astigmatism=astigmatism,
                                     astigmatism_angle=astigmatism_angle, display=False)
        ctf_series.append(ctf)

    dqe = create_detector_response(camera_type, 'DQE', xp.zeros((image_size, image_size)), voltage=voltage,
                                            folder=camera_folder)
    mtf = create_detector_response(camera_type, 'MTF', xp.zeros((image_size, image_size)), voltage=voltage,
                                            folder=camera_folder)

    # joblib automatically memory maps a numpy array to child processes
    print(f'Projecting the model with {nodes} processes')

    verbosity = 55  # set to 55 for debugging, 11 to see progress, 0 to turn off output
    results = Parallel(n_jobs=nodes, verbose=verbosity, prefer="threads") \
        (delayed(parallel_project)(rotation_volume, i, image_size, pixel_size, msdz, n_slices, ctf,
                                   dose_per_tilt, dqe, mtf, voltage, binning=binning, translation=translation,
                                   rotation=(.0, angle, .0), solvent_potential=solvent_potential,
                                   solvent_absorption=solvent_amplitude, sigma_damage=sigma_damage,
                                   ice_voxels=box_height)
         for i, (angle, translation, ctf) in enumerate(zip(angles, translations, ctf_series)))

    sys.stdout.flush()

    if results.count(None) == 0:
        print('All projection processes finished successfully')
    else:
        print(f'{results.count(None)} rotation processes did not finish successfully')

    # write (noisefree) projections as mrc stacks
    filename_nf = os.path.join(save_path, 'noisefree_projections.mrc')
    filename_pr = os.path.join(save_path, 'projections.mrc')
    write(filename_nf, xp.stack([n for (n, p) in results], axis=2))
    write(filename_pr, xp.stack([p for (n, p) in results], axis=2))

    # store alignment information
    # len(angles) is the number of files that we have
    alignment                       = xp.zeros(len(angles), dtype=dar)
    # IMPORTANT: correct angles by -1 to get the right reconstruction
    alignment['TiltAngle']          = -1.0 * xp.array(angles)
    alignment['Magnification']      = xp.repeat(1.0, len(angles))
    for i in range(len(angles)):
        alignment['FileName'][i]    = os.path.join(save_path, 'projections', f'synthetic_{i+1}.mrc')

    # write the alignment file
    filename_align                      = os.path.join(save_path, 'alignment_simulated.txt')
    savestar(filename_align, alignment, fmt=fmtAlignmentResults, header=HEADER_ALIGNMENT_RESULTS)

    return


def generate_frame_series_cpu(save_path, n_frames=20, nodes=1, image_size=None, pixel_size=1E-9,
                              binning=1, dose=80, voltage=300E3, spherical_aberration=2.7E-3,
                              chromatic_aberration=2.7E-3, energy_spread=0.7, illumination_aperture=0.030E-3,
                              objective_diameter=100E-6, focus_length=4.7E-3, astigmatism=0.0E-9, astigmatism_angle=0,
                              msdz=1E-9, defocus=2E-6, sigma_damage=0.0, camera_type='K2SUMMIT', camera_folder='',
                              solvent_potential=physics.V_WATER, absorption_contrast=False):
    """
    Creating a frame series for the initial grand model by applying stage drift (translation) and some rotations for
    each frame, and the calculating the sample projection in the microscope. Additionally apply increasing random noise
    to the images to imitate particle degradation through beam interaction. I cannot do this incremental as I want
    to parallelize the rotation+projection, but the overall effect on signal should be identical.

    @param save_path
    @param n_frames:
    @param image_size:
    @param pixel_size:
    @param binning:
    @param dose:
    @param voltage:
    @param spherical_aberration:
    @param multislice:
    @param msdz:
    @param defocus:
    @param sigma_decay_CTF:
    @param camera_type:
    @param camera_folder:
    @param random_gradient:
    @param solvent_potential:
    @param solvent_factor:
    @param absorption_contrast:
    @return:
    """
    from pytom.basic.datatypes import DATATYPE_ALIGNMENT_RESULTS as dar
    from pytom.basic.datatypes import fmtAlignmentResults, HEADER_ALIGNMENT_RESULTS
    from pytom.gui.guiFunctions import savestar
    from pytom.simulation.microscope import create_detector_response, create_complex_ctf
    from pytom.tompy.io import read_mrc, write
    from joblib import Parallel, delayed

    # NOTE; Parameter defocus specifies the defocus at the bottom of the model!

    # grab model
    filename_gm_real = os.path.join(save_path, 'grandmodel.mrc')
    filename_gm_imag = os.path.join(save_path, 'grandmodel_imag.mrc')
    grandcell = read_mrc(filename_gm_real)
    if absorption_contrast:
        grandcell = grandcell.astype(xp.complex64)
        grandcell.imag = read_mrc(filename_gm_imag)
        # calculate the absorption for amorphous ice at the specified voltage
        solvent_amplitude = physics.potential_amplitude(physics.AMORPHOUS_ICE_DENSITY, physics.WATER_MW, voltage)
        print(f'solvent absorption = {solvent_amplitude:.3f}')
    else:
        grandcell = grandcell.astype(xp.float32)
        solvent_amplitude = 0.0

    # extract size
    box_size = grandcell.shape[0]
    box_height = grandcell.shape[2]

    if image_size is None:
        image_size = box_size

    # confirm image_size is valid
    assert image_size <= box_size, 'Specified projection image size is invalid as it is larger than the model dimension.'

    # adjust defocus
    zheight = box_height * pixel_size  # thickness of rotation volume in nm
    defocus -= (zheight / 2)  # center defocus value at tilt angle

    # Check if msdz is viable, else correct it
    if msdz % pixel_size != 0:
        # Round the msdz to the closest integer mutltiple of the pixel size
        round_up = xp.round(msdz % pixel_size / pixel_size)
        if round_up:
            msdz += (pixel_size - msdz % pixel_size)
        else:
            msdz -= (msdz % pixel_size)
        print(f'We adjusted msdz to {msdz*1E9:.3f} nm to make it an integer multiple of pixel size.')

    # Determine the number of slices
    if msdz > zheight:
        n_slices = 1
        msdz = zheight
        print('The multislice step can not be larger than volume thickness. Use only one slice.\n')
    elif msdz < pixel_size:
        n_slices = box_height
        msdz = pixel_size
        print(
            'The multislice steps can not be smaller than the pixel size. Use the pixel size instead for the step size.\n')
    else:
        n_slices = int(xp.ceil(xp.around(zheight / msdz, 3)))
    print('Number of slices for multislice: ', n_slices)

    # determine dose per frame
    dose_per_frame = dose / n_frames

    # First generate stage drift and in-plane rotation, stage drift is a set of correlated translations across the
    # number of frames. In MotionCorr2 paper accumulated motion for 20S proteasome dataset is 11A across the whole
    # frame series.
    # First generate global motion and global direction of motion.
    global_motion = xp.random.normal(10, 3) # normal around mean 10 A and std 3A
    average_motion_per_frame = global_motion / n_frames
    global_angle = xp.random.uniform(0,360) # random angle from uniform
    translations, cumulative_translations, translations_voxel = [], [], []
    x, y = 0, 0
    for i in range(n_frames):
        # randomly vary the motion per frame and angle
        motion_i = xp.random.normal(average_motion_per_frame, average_motion_per_frame/2)
        angle_i = xp.random.normal(global_angle, 20)
        # decompose motion into x and y translation
        y_i = xp.sin(angle_i * xp.pi / 180) * motion_i
        x_i = xp.cos(angle_i * xp.pi / 180) * motion_i
        # only seem to need cumulative_translations
        translations.append((x_i, y_i, 0)) # append the translation for the frame as a tuple
        y += y_i
        x += x_i
        cumulative_translations.append((x,y, 0)) # translation for z coordinate as we are referring to volumes
        translations_voxel.append((x*1E-10 / pixel_size, y*1E-10 / pixel_size, 0))

    # write motion trajectory to a png file for debugging
    # fig, ax = plt.subplots(2)
    # ax[0].plot([x for (x,y,z) in cumulative_translations], [y for (x,y,z) in cumulative_translations], label='trajectory')
    # ax[0].set_xlabel('x (A)')
    # ax[0].set_ylabel('y (A)')
    # ax[0].legend()
    # ax[1].plot([x for (x,y,z) in translations_voxel], [y for (x,y,z) in translations_voxel], label='trajectory')
    # ax[1].set_xlabel('x (voxels)')
    # ax[1].set_ylabel('y (voxels)')
    # ax[1].legend()
    # plt.savefig(f'{save_path}/global_motion.png')
    # plt.close()

    # get the contrast transfer function
    ctf = create_complex_ctf((image_size, image_size), pixel_size, defocus, voltage=voltage,
                                   Cs=spherical_aberration, Cc=chromatic_aberration, energy_spread=energy_spread,
                                   illumination_aperture=illumination_aperture, objective_diameter=objective_diameter,
                                   focus_length=focus_length, astigmatism=astigmatism,
                                   astigmatism_angle=astigmatism_angle, display=False)

    dqe = create_detector_response(camera_type, 'DQE', xp.zeros((image_size, image_size)), voltage=voltage,
                                            folder=camera_folder)
    mtf = create_detector_response(camera_type, 'MTF', xp.zeros((image_size, image_size)), voltage=voltage,
                                            folder=camera_folder)

    # joblib automatically memory maps a numpy array to child processes
    print(f'Projecting the model with {nodes} processes')

    verbosity = 55  # set to 55 for debugging, 11 to see progress, 0 to turn off output
    results = Parallel(n_jobs=nodes, verbose=verbosity, prefer="threads") \
        (delayed(parallel_project)(grandcell, frame, image_size, pixel_size, msdz, n_slices, ctf,
                                   dose_per_frame, dqe, mtf, voltage, binning=binning, translation=shift, rotation=(.0,.0,.0),
                                   solvent_potential=solvent_potential, solvent_absorption=solvent_amplitude,
                                   sigma_damage=sigma_damage, ice_voxels=None)
         for frame, shift in enumerate(translations_voxel))

    sys.stdout.flush()
    if results.count(None) == 0:
        print('All projection processes finished successfully')
    else:
        print(f'{results.count(None)} rotation processes did not finish successfully')

    # write (noisefree) projections as mrc stacks
    filename_nf = os.path.join(save_path, 'noisefree_projections.mrc')
    filename_pr = os.path.join(save_path, 'projections.mrc')
    write(filename_nf, xp.stack([n for (n,p) in results], axis=2))
    write(filename_pr, xp.stack([p for (n, p) in results], axis=2))

    # Store translations as reference for model
    # len(angles) is the number of files that we have
    alignment                       = xp.zeros(n_frames, dtype=dar)
    alignment['AlignmentTransX']    = xp.array([x for (x,y,z) in cumulative_translations])
    alignment['AlignmentTransY']    = xp.array([y for (x,y,z) in cumulative_translations])
    alignment['Magnification']      = xp.repeat(1.0, n_frames)
    for i in range(n_frames):
        alignment['FileName'][i]    = os.path.join(save_path, 'projections', f'synthetic_{i+1}.mrc')

    # Write the alignment file as a text file
    filename_align                      = os.path.join(save_path, 'alignment_simulated.txt')
    savestar(filename_align, alignment, fmt=fmtAlignmentResults, header=HEADER_ALIGNMENT_RESULTS)
    return


def FSS(fvolume1, fvolume2, numberBands, verbose=False):
    """
    algorithm FSS = Fourier Shell Scaling
    Scale the values of fvolume1 to the values in fvolume2 per band in fourier space.
    """
    from pytom.tompy.correlation import meanUnderMask
    from pytom.simulation.support import bandpass_mask

    assert fvolume1.shape == fvolume2.shape, "volumes not of same size"
    assert len(set(fvolume1.shape)) == 1, "volumes are not perfect cubes"

    if verbose:
        print(f'shape of images is: {fvolume1.shape}')

    increment = int(fvolume1.shape[0] / 2 * 1 / numberBands)
    band = [-1, -1]

    output = xp.zeros(fvolume1.shape)

    for i in range(0, fvolume1.shape[0] // 2, increment):

        band[0] = i
        band[1] = i + increment

        if verbose:
            print('Band : ', band)

        if band[1] >= fvolume1.shape[0] // 2:
            bandpass = bandpass_mask(fvolume1.shape, 0, high=band[0])
            bandpass = (bandpass == 0) * 1
        else:
            bandpass = bandpass_mask(fvolume1.shape, band[0], band[1])

        if i == 0:
            # remove center point from mask
            c = bandpass.shape[0] // 2
            bandpass[c, c] = 0
            # scale center separately as center adjustment has large influence on the image
            output[c, c] = fvolume1[c, c] / (fvolume1[c, c] / fvolume2[c, c])

        n = bandpass.sum()
        # get mean amplitude of each band
        m1 = meanUnderMask(fvolume1, bandpass, p=n)
        m2 = meanUnderMask(fvolume2, bandpass, p=n)

        # scale the values inside the band of volume1
        outband = fvolume1 * bandpass / (m1 / m2)

        output += outband

    del outband, bandpass
    return output


def scale_image(volume1, volume2, numberBands):
    """
    Scale amplitudes in fourier space of volume1 to those of volume2 using fourier shells.
    Output is the real space result of scaled volume1.
    """
    assert volume1.shape == volume2.shape, "volumes not of same size"
    assert len(set(volume1.shape)) == 1, "volumes are not perfect cubes"

    # scale the amplitudes of 1 to those of 2 using fourier shells
    scaled_amp = FSS(xp.abs(xp.fft.fftshift(xp.fft.fftn(volume1))), xp.abs(xp.fft.fftshift(xp.fft.fftn(volume2))),
                     numberBands, verbose=False)
    # construct the output volume with the scaled amplitudes and phase infomation of volume 1
    fout = xp.fft.ifftn(xp.fft.ifftshift(scaled_amp) * xp.exp(1j * xp.angle(xp.fft.fftn(volume1))))

    return fout.real


def parallel_scale(number, projection, example, pixel_size, example_pixel_size, binning, make_even_factor):
    from pytom.tompy.transform import resize

    print(projection.shape, example.shape)

    print(f' -- scaling projection {number+1}')
    if pixel_size != (example_pixel_size * binning):
        # magnify or de-magnify if the pixel size does not match yet
        print('(de)magnifying pixel size')
        example = resize(example, (example_pixel_size * binning) / pixel_size, interpolation='Spline')

    # prevent issues later on with binning in case experimental and simulated pixel size do not match
    if example.shape[0] % (2*make_even_factor):
        example = example[:-(example.shape[0] % (2*make_even_factor)), :-(example.shape[0] % (2*make_even_factor))]

    if projection.shape != example.shape:
        # crop the largest
        if projection.shape[0] > example.shape[0]:
            cut = (projection.shape[0] - example.shape[0]) // 2
            projection = projection[cut:-cut, cut:-cut]
        else:
            cut = (example.shape[0] - projection.shape[0]) // 2
            example = example[cut:-cut, cut:-cut]

    print('using FSS to scale amplitudes')
    return number, scale_image(projection, example, projection.shape[0] // 4)


def scale_projections(save_path, pixel_size, example_folder, example_pixel_size, binning, nodes,
                      make_even_factor):
    from pytom.tompy.transform import resize
    from pytom.simulation.support import reduce_resolution
    from pytom.tompy.io import read_mrc, write
    from joblib import Parallel, delayed

    # generate list of all the possible example projections in the provided parameter example_folder
    files = [f for f in os.listdir(example_folder) if os.path.isfile(os.path.join(example_folder, f))]
    random_file = xp.random.randint(0, len(files))
    print(f'Selected {files[random_file]} for experimental amplitude scaling.')

    filename_pr = os.path.join(save_path, 'projections.mrc')
    filename_ex = os.path.join(example_folder, files[random_file])
    projections         = read_mrc(filename_pr)
    example_projections = read_mrc(filename_ex)

    # assert projections.shape[2] == example_projections.shape[2], 'not enough or too many example projections'
    sim_size = projections.shape
    exp_size = example_projections.shape
    if sim_size[2] > exp_size[2]:  # make sure number of poejection angles match
        diff = sim_size[2] - exp_size[2]
        new = xp.zeros((exp_size[0], exp_size[1], sim_size[2]))
        new[..., diff//2: - (diff//2+diff%2)] = example_projections
        new[..., :diff//2] = example_projections[..., :diff//2]
        new[..., -(diff//2+diff%2):] = example_projections[..., -(diff//2+diff%2):]
        example_projections = new
    elif sim_size[2] < exp_size[2]:
        diff = exp_size[2] - sim_size[2]
        example_projections = example_projections[..., diff//2: -(diff//2 + diff%2)]

    # joblib automatically memory maps a numpy array to child processes
    print(f'Scaling projections with {nodes} processes')

    verbosity = 55  # set to 55 for debugging, 11 to see progress, 0 to turn off output
    results = Parallel(n_jobs=nodes, verbose=verbosity, prefer="threads") \
        (delayed(parallel_scale)(i, projections[:, :, i].squeeze(),
                                 resize(reduce_resolution(example_projections[:, :, i], 1, 2*binning),
                                        1/binning, interpolation='Spline').squeeze(),
                                 pixel_size, example_pixel_size, binning, make_even_factor)
         for i in range(projections.shape[2]))

    sys.stdout.flush()

    if results.count(None) == 0:
        print('All scaling processes finished successfully')
    else:
        print(f'{results.count(None)} scaling processes did not finish successfully')

    new_projections = xp.dstack(tuple([r for (i,r) in results]))

    filename_scaled = os.path.join(save_path, 'projections_scaled.mrc')
    write(filename_scaled, new_projections)

    return


def reconstruct_tomogram(save_path, weighting=-1, reconstruction_bin=1,
                         filter_projections=False, use_scaled_projections=False):
    """
    Reconstruction of simulated tilt series into a tomogram. First creates an alignemnt file, then uses weighted back
    projection to make a reconstruction.

    @param prefix: file name prefix
    @type prefix: string
    @param suffix: file extension
    @type suffix: basestring
    @param start_idx: starting index number of first tilt image in file names
    @type start_idx: int
    @param end_idx: ending index of last tilt image in file names
    @type end_idx: int
    @param vol_size: size of the reconstruction volume
    @type vol_size: [int, int, int]
    @param angles: list of angles of the tilt projections
    @type angles: [float, float, ...]
    @param output_folder: output directory path
    @type output_folder: string
    @param weighting: weighting factor, either -1, 0, or 1
    @type weighting: int

    @return: -
    @rtype: None

    @author: Gijs van der Schot, Marten Chaillet
    """
    from pytom.reconstruction.reconstructionStructures import Projection, ProjectionList
    from pytom.tompy.transform import resize
    from pytom.simulation.support import reduce_resolution
    from pytom.tompy.io import read_mrc, write

    # create folder for individual projections
    projection_folder = os.path.join(save_path, 'projections')
    if not os.path.exists(projection_folder):
        os.mkdir(projection_folder)

    # Possibly apply low pass filter at this point
    if use_scaled_projections:
        filename_sc = os.path.join(save_path, 'projections_scaled.mrc')
        projections = read_mrc(filename_sc)
    else:
        filename_pr = os.path.join(save_path, 'projections.mrc')
        projections = read_mrc(filename_pr)

    if filter_projections:
        for i in range(projections.shape[2]):
            # 2.3 corresponds to circular filter with width 0.9 of half of the image
            projection_scaled = reduce_resolution(projections[:,:,i].squeeze(), 1.0, 2.3 * reconstruction_bin)
            filename_single = os.path.join(save_path, 'projections', f'synthetic_{i+1}.mrc')
            write(filename_single, projection_scaled)
    else:
        for i in range(projections.shape[2]):
            filename_single = os.path.join(save_path, 'projections', f'synthetic_{i+1}.mrc')
            write(filename_single, projections[:,:,i])

    # size and volume shape of the reconstruction
    size_reconstruction = projections.shape[0] // reconstruction_bin
    vol_size = [size_reconstruction, ] * 3
    # alignment file name for reconstruction
    filename_align = os.path.join(save_path, 'alignment_simulated.txt')

    if reconstruction_bin == 1:
        filename_output = os.path.join(save_path, 'reconstruction.em')
    else:
        filename_output = os.path.join(save_path, f'reconstruction_bin{reconstruction_bin}.em')
    # IF EM alignment file provided, filters applied and reconstruction will be identical.
    projections = ProjectionList()
    vol = projections.reconstructVolume(dims=vol_size, reconstructionPosition=[0,0,0], binning=reconstruction_bin,
                                        applyWeighting=weighting, alignResultFile=filename_align)
    vol.write(filename_output)
    os.system(f'em2mrc.py -f {filename_output} -t {os.path.dirname(filename_output)}')
    os.system(f'rm {filename_output}')

    # Adjust ground truth data to match cuts after scaling or after binning for reconstruction
    filename_gm = os.path.join(save_path, 'grandmodel.mrc')
    filename_cm = os.path.join(save_path, 'class_mask.mrc')
    filename_om = os.path.join(save_path, 'occupancy_mask.mrc')
    filename_gm_bin = os.path.join(save_path, f'grandmodel_bin{reconstruction_bin}.mrc')
    filename_cm_bin = os.path.join(save_path, f'class_mask_bin{reconstruction_bin}.mrc')
    filename_om_bin = os.path.join(save_path, f'occupancy_mask_bin{reconstruction_bin}.mrc')

    cell = read_mrc(filename_gm)
    # find cropping indices
    lind = (cell.shape[0] - reconstruction_bin * size_reconstruction)//2
    rind = cell.shape[0] - lind
    cell = resize(reduce_resolution(cell[lind:rind,lind:rind,:], 1, 2*reconstruction_bin), 1/reconstruction_bin,
                  interpolation='Spline')
    write(filename_gm_bin, cell)
    # bin class mask and bbox needed for training
    cell = read_mrc(filename_cm)
    write(filename_cm_bin, downscale_class_mask(cell[lind:rind,lind:rind,:], reconstruction_bin))
    # bin occupancy mask as well
    cell = read_mrc(filename_om)
    write(filename_om_bin, downscale_class_mask(cell[lind:rind,lind:rind,:], reconstruction_bin))

    # create particle locations bin2 file
    adjusted_ground_truth = ''
    filename_loc = os.path.join(save_path,'particle_locations.txt')
    filename_loc_bin = os.path.join(save_path, f'particle_locations_bin{reconstruction_bin}.txt')
    with open(filename_loc, 'r') as fin:
        line = fin.readline()
        while line:
            data = line.split()
            data[1] = int(data[1]) // reconstruction_bin - lind
            data[2] = int(data[2]) // reconstruction_bin - lind
            data[3] = int(data[3]) // reconstruction_bin - lind
            # only add the location back if the center of particle is still inside the box after binning
            if 0 <= data[1] < rind and 0 <= data[2] < rind and 0 <= data[3] < rind:
                data[1] = str(data[1])
                data[2] = str(data[2])
                data[3] = str(data[3])
                adjusted_ground_truth += ' '.join(data) + '\n'
            line = fin.readline()
    with open(filename_loc_bin, 'w') as fout:
        fout.write(adjusted_ground_truth)

    return


if __name__ == '__main__':
    from pytom.gui.guiFunctions import loadstar, datatype
    from ast import literal_eval
    # Use tracemalloc to record the peak memory usage of the program
    tracemalloc.start()

    # Read config
    config = configparser.ConfigParser()
    try:
        if len(sys.argv) > 1:
            config_given = sys.argv[1]
            if config_given and os.path.exists(config_given):
                print(f'\nLoading a given configuration file: {config_given}')
                config.read_file(open(config_given))
        else:
            print(f'\nLoading default configuration file: pytom/simulation/simulation.conf')
            config.read_file(open('simulation.conf'))
    except Exception as e:
        print(e)
        raise Exception('Could not open config file.')

    print('Configuration sections:', config.sections())

    # Set simulation parameters
    try:
        output_folder           = config['General']['OutputFolder']
        simulator_mode          = config['General']['Mode']
        device                  = config['General']['Device']
        nodes                   = config['General'].getint('Nodes')
        model_ID                = config['General'].getint('ModelID')
        seed                    = config['General'].getint('Seed')
        pixel_size              = config['General'].getfloat('PixelSize') * 1E-10 # pixel_size in nm
        binning                 = config['General'].getint('Binning')
        solvent_potential       = config['General'].getfloat('SolventConstant')
        absorption_contrast     = config['General'].getboolean('AbsorptionContrast')
        voltage                 = config['General'].getfloat('Voltage') * 1E3  # voltage in keV
        # voltage and pixelsize are needed for model generation and projection, thus general parameters

        # adjust pixel size with binning factor
        pixel_size *= binning

        # ensure simulator mode and device are valid options
        if (simulator_mode in ['TiltSeries', 'FrameSeries']) or (device in ['CPU', 'GPU']):
            print(f'Generating model {model_ID} on {device} in folder {output_folder}')
        else:
            print('Invalid entry for simulator mode or device in config.')
            sys.exit(0)
    except Exception as e:
        print(e)
        raise Exception('Missing general parameters in config file.')

    if 'GenerateModel' in config.sections():
        try:
            # We assume the particle models are in the desired voxel spacing for the pixel size of the simulation!
            particle_folder     = config['GenerateModel']['ParticleFolder']
            listpdbs            = literal_eval(config['GenerateModel']['Models'])
            listmembranes       = literal_eval(config['GenerateModel']['MembraneModels'])
            size                = config['GenerateModel'].getint('Size')
            placement_size      = config['GenerateModel'].getint('PlacementSize')
            # parse range of ice thickness, provided in nm
            thickness           = draw_range(literal_eval(config['GenerateModel']['Thickness']), float, 'Thickness') * 1E-9
            thickness_voxels    = int(thickness / pixel_size) # calculate thickness in number of voxels!
            # make even number to solve tomogram reconstruction mismatch bug
            thickness_voxels -= (thickness_voxels % 4)
            # gold markers
            number_of_markers   = draw_range(literal_eval(config['GenerateModel']['NumberOfMarkers']), int,
                                           'NumberOfMarkers')
            # parse range of number of particles
            number_of_particles = draw_range(literal_eval(config['GenerateModel']['NumberOfParticles']), int,
                                             'NumberOfParticles')
            number_of_membranes = draw_range(literal_eval(config['GenerateModel']['NumberOfMembranes']), int,
                                             'NumberOfMembranes')
            # pass motion blur in angstrom units
            # motion_blur         = config['GenerateModel'].get_float('MotionBlurResolution')
            motion_blur         = draw_range(literal_eval(config['GenerateModel']['MotionBlurResolution']), float,
                                             'MotionBlurResolution')
        except Exception as e:
            print(e)
            raise Exception('Missing generate model parameters in config file.')

    if 'Microscope' in config.sections():
        try:
            camera                  = config['Microscope']['Camera']
            camera_folder           = config['Microscope']['CameraFolder']
            defocus                 = draw_range(literal_eval(config['Microscope']['Defocus']), float, 'Defocus') * 1E-6
            electron_dose           = draw_range(literal_eval(config['Microscope']['ElectronDose']), float, 'ElectronDose')
            spherical_aberration    = config['Microscope'].getfloat('SphericalAberration') * 1E-3
            chromatic_aberration    = config['Microscope'].getfloat('ChromaticAberration') * 1E-3
            energy_spread           = config['Microscope'].getfloat('EnergySpread')
            illumination_aperture   = config['Microscope'].getfloat('IlluminationAperture') * 1E-3
            objective_diameter      = config['Microscope'].getfloat('ObjectiveDiameter') * 1E-6
            focus_length            = config['Microscope'].getfloat('FocalDistance') * 1E-3
            astigmatism             = config['Microscope'].getfloat('Astigmatism') * 1E-9
            astigmatism_angle       = config['Microscope'].getfloat('AstigmatismAngle')
        except Exception as e:
            print(e)
            raise Exception('Missing microscope parameters in config file.')

    if simulator_mode in config.sections():
        try:
            # first read common parameters between tilt and frame series
            image_size      = config[simulator_mode].getint('ImageSize')
            msdz            = config[simulator_mode].getfloat('MultisliceStep') * 1E-9
            beam_damage     = config[simulator_mode].getfloat('BeamDamage')
            # mode specific parameters
            if simulator_mode == 'TiltSeries':
                metadata            = loadstar(config['TiltSeries']['MetaFile'], dtype=datatype)
                angles              = metadata['TiltAngle'] # in degrees
            elif simulator_mode == 'FrameSeries':
                number_of_frames    = config['FrameSeries'].getint('NumberOfFrames')
        except Exception as e:
            print(e)
            raise Exception(f'Missing {simulator_mode} parameters in config file.')

    if 'ScaleProjections' in config.sections():
        try:
            example_folder      = config['ScaleProjections']['ExampleFolder']
            example_pixel_size  = config['ScaleProjections'].getfloat('ExamplePixelSize')
            # If experimental and simulated projections have different size, we need to crop. This should be done with
            # care if the option binning is set for reconstructions, because in that case the ground truth data needs
            # to be binned and cropped as well. Uneven size of the volume means the ground truth data will be shifted
            # by half a pixel compared to the reconstruction. This options makes sure that this not happen.
            make_even_factor    = config['ScaleProjections'].getint('EvenSizeFactor')
        except Exception as e:
            print(e)
            raise Exception('Missing experimental projection scaling parameters.')

    if 'TomogramReconstruction' in config.sections():
        try:
            weighting               = config['TomogramReconstruction'].getint('Weighting')
            reconstruction_bin      = config['TomogramReconstruction'].getint('Binning')
            filter_projections      = config['TomogramReconstruction'].getboolean('FilterProjections')
            use_scaled_projections  = config['TomogramReconstruction'].getboolean('UseScaledProjections')
        except Exception as e:
            print(e)
            raise Exception('Missing tomogram reconstruction parameters in config file.')

    # Create directories and logger
    save_path = os.path.join(output_folder, f'model_{model_ID}')
    if not os.path.exists(save_path):
        os.mkdir(save_path)

    logging.basicConfig(filename='{}/simulator-{date:%Y-%m-%d_%H:%M:%S}.log'.format(save_path,
                                                                date=datetime.datetime.now()), level=logging.INFO)
    config_logger = ConfigLogger(logging)
    config_logger(config)

    logging.info('Values of parameters that we randomly vary per simulation (only for generate model and generate projections):')
    if 'GenerateModel' in config.sections():
        logging.info(f'model thickness = {thickness_voxels*pixel_size*1E9:.2f}nm (adjusted to be an even number of voxels)')
        logging.info(f'# of particles = {number_of_particles}')
        logging.info(f'# of markers = {number_of_markers}')
        logging.info(f'# of membranes = {number_of_membranes}')
    if 'Microscope' in config.sections():
        logging.info(f'defocus = {defocus*1E6:.2f}um')
        logging.info(f'electron dose = {electron_dose} e-/A^2')

    # Generate a grand model
    if 'GenerateModel' in config.sections():
        # set seed for random number generation
        xp.random.seed(seed)
        random.seed(seed)

        print('\n- Generating grand model')
        generate_model(particle_folder, save_path, listpdbs, listmembranes,
                       pixel_size           =pixel_size,
                       size                 =size,
                       thickness            =thickness_voxels,
                       placement_size       =placement_size,
                       solvent_potential    =solvent_potential,
                       number_of_particles  =number_of_particles,
                       number_of_markers    =number_of_markers,
                       absorption_contrast  =absorption_contrast,
                       voltage              =voltage,
                       number_of_membranes  =number_of_membranes,
                       sigma_motion         =motion_blur)

    if simulator_mode in config.sections() and simulator_mode == 'TiltSeries':
        # set seed for random number generation
        xp.random.seed(seed)
        random.seed(seed)
        # Grab the ice thickness from the initial model in case program is only executed for projections
        print('\n- Generating projections')
        if device == 'CPU':
            generate_tilt_series_cpu(save_path, angles,
                                      nodes                 =nodes,
                                      image_size            =image_size,
                                      rotation_box_height   =None, # will automatically calculate fitting size if None
                                      pixel_size            =pixel_size,
                                      binning               =binning,
                                      dose                  =electron_dose,
                                      voltage               =voltage,
                                      spherical_aberration  =spherical_aberration,
                                      chromatic_aberration  =chromatic_aberration,
                                      energy_spread         =energy_spread,
                                      illumination_aperture =illumination_aperture,
                                      objective_diameter    =objective_diameter,
                                      focus_length          =focus_length,
                                      astigmatism           =astigmatism,
                                      astigmatism_angle     =astigmatism_angle,
                                      msdz                  =msdz,
                                      defocus               =defocus,
                                      sigma_damage          =beam_damage,
                                      camera_type           =camera,
                                      camera_folder         =camera_folder,
                                      solvent_potential     =solvent_potential,
                                      absorption_contrast   =absorption_contrast)
        elif device == 'GPU':
            print('This option needs to be implemented.')
            sys.exit(0)
        else:
            print('Invalid device type.')
            sys.exit(0)
    elif simulator_mode in config.sections() and simulator_mode == 'FrameSeries':
        xp.random.seed(seed)
        random.seed(seed)
        print('\n- Generate frame series projections')
        if device == 'CPU':
            generate_frame_series_cpu(save_path,
                                      n_frames              =number_of_frames,
                                      nodes                 =nodes,
                                      image_size            =image_size,
                                      pixel_size            =pixel_size,
                                      binning               =binning,
                                      dose                  =electron_dose,
                                      voltage               =voltage,
                                      spherical_aberration  =spherical_aberration,
                                      chromatic_aberration  =chromatic_aberration,
                                      energy_spread         =energy_spread,
                                      illumination_aperture =illumination_aperture,
                                      objective_diameter    =objective_diameter,
                                      focus_length          =focus_length,
                                      astigmatism           =astigmatism,
                                      astigmatism_angle     =astigmatism_angle,
                                      msdz                  =msdz,
                                      defocus               =defocus,
                                      sigma_damage          =beam_damage,
                                      camera_type           =camera,
                                      camera_folder         =camera_folder,
                                      solvent_potential     =solvent_potential,
                                      absorption_contrast   =absorption_contrast)
        elif device == 'GPU':
            print('This option needs to be implemented.')
            sys.exit(0)
        else:
            print('Invalid device type.')
            sys.exit(0)

    if 'ScaleProjections' in config.sections():
        # set seed for random number generation
        xp.random.seed(seed)
        random.seed(seed)
        print('\n- Scaling projections with experimental data')
        scale_projections(save_path, pixel_size * 1E10, example_folder,
                                            example_pixel_size, binning, nodes, make_even_factor)

    if 'TomogramReconstruction' in config.sections():
        print('\n- Reconstructing tomogram')
        reconstruct_tomogram(save_path,
                             weighting=weighting,
                             reconstruction_bin=reconstruction_bin,
                             filter_projections=filter_projections,
                             use_scaled_projections=use_scaled_projections)

    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
    tracemalloc.stop()
