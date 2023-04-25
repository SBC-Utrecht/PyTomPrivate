#!/usr/bin/env python

import os
import sys
import pytom.simulation.physics as physics
from pytom.agnostic.io import write
from pytom.simulation.template import generate_template_from_pdb, generate_template_from_map


if __name__ == '__main__':
    # todo template generation could be updated with an option for spherical mask generation

    # parameters: file_path, destination, spacing, binning (optional, default is 1), solvent_correction (optional),
    # solvent_density (optional, default 0.93),
    # apply_ctf_correction (optional), defocus (optional, default is 3 um, negative is overfocus),
    # amplitude contrast (optional, default 0.07), voltage (optional, default is 300 keV),
    # Cs (optional, default is 2.7 mm), ctf_decay (optional, default is 0.4), display_ctf (optional),
    # apply_low_pass (optional), resolution_filter (optional, default is 2*spacing*binning), box_size (optional)

    # IN ORDER TO FUNCTION, SCRIPT REQUIRES INSTALLATION OF PYTOM (and dependencies), CHIMERA

    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2

    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Generate a template for template matching!',
        authors='Marten Chaillet',
        options=[
            ScriptOption2(['-f', '--file'], 'Protein structure file, either pdb/cif or mrc/em.', 'file', 'required'),
            ScriptOption2(['-d', '--destination'], 'Folder where output should be stored.', 'directory', 'required'),
            ScriptOption2(['-o', '--output-name'], 'Name of file to write as output, with extension (mrc or em).',
                          'string', 'optional'),
            ScriptOption2(['-s', '--spacing'], 'The pixel spacing of original projections of the dataset in A,'
                                             ' e.g. 2.62', 'float', 'required'),
            ScriptOption2(['--map-spacing'], 'Additional spacing of the EM map, provide in addition to spacing '
                                             'argument.', 'float', 'optional'),
            ScriptOption2(['-b', '--binning'], 'Number of times to bin the template. Default is 1 (no binning). If '
                                             'set to 2 with a spacing of 2.62 the resulting voxel size will '
                                             'be 5.24', 'int', 'optional', 1),
            ScriptOption2(['--modify-structure'], 'Activate to call Chimera for adding hydrogen, symmetry and'
                                                      'removing water molecules from the structure.', 'no arguments',
                          'optional'),
            ScriptOption2(['--solvent-correction'], 'Whether to exclude solvent around each atom as a '
                                                        'correction of the potential, options: gaussian, masking.',
                          'string', 'optional'),
            ScriptOption2(['-r', '--solvent-density'], 'Density of solvent, value should not be higher than 1.35 as'
                                                     ' that is the density of proteins. Default is 0.93 g/cm^3.',
                          'float', 'optional', physics.AMORPHOUS_ICE_DENSITY),
            ScriptOption2(['-c', '--ctf-correction'], 'Correct the volume by applying a CTF. Default parameters are '
                                                  'defocus 3 um, amplitude contrast 0.07, voltage 300 keV, '
                                                  'spherical abberation (Cs) 2.7 mm, sigma of gaussian decay 0.4, '
                                                  'optionally plot the CTF to inspect.', 'no arguments', 'optional'),
            ScriptOption2(['-z', '--defocus'], 'Defocus in um (negative value is overfocus).', 'float', 'optional', 3.0),
            ScriptOption2(['-a', '--amplitude-contrast'], 'Amplitude contrast fraction.', 'float', 'optional', 0.07),
            ScriptOption2(['-v', '--voltage'], 'Acceleration voltage in keV', 'float', 'optional', 300),
            ScriptOption2(['--Cs'], 'Spherical abberration in mm.', 'float', 'optional', 2.7),
            ScriptOption2(['--decay'], 'Sigma of gaussian CTF decay function, 0.4 default.', 'float', 'optional', 0.4),
            ScriptOption2(['--flip-phase'], 'Apply phase flipped CTF, corresponding to tomogram which has been CTF '
                                            'corrected through phase-flipping.', 'no arguments', 'optional'),
            ScriptOption2(['--cut-first-zero'], 'Cut ctf after first zero crossing', 'no arguments', 'optional'),
            ScriptOption2(['--plot'], 'Give this option for plotting the CTF for visual inspection.', 'no arguments',
                          'optional'),
            ScriptOption2(['-l', '--lpf-resolution'], 'Specify the resolution of the low pass filter that is applied.'
                                                    'The default value is 2 x spacing x binning (in angstrom), a '
                                                    'smaller resolution than this cannot be selected.', 'float',
                          'optional'),
            ScriptOption2(['-x', '--xyz'], 'Specify a desired size for the output box of the template in number of '
                                           'pixels. By default the molecule is placed in a box with 30A overhang. '
                                           'This usually does not offer enough room to apply a spherical mask.',
                          'int', 'optional'),
            ScriptOption2(['-i', '--invert'], 'Multiplies template by -1. WARNING not needed if ctf with defocus is '
                                              'already applied!', 'no arguments', 'optional'),
            ScriptOption2(['-m', '--mirror'], 'Produce a mirrored and non-mirrored version.', 'no arguments',
                          'optional'),
            ScriptOption2(['--cores'], 'Number of cores to run template generation on.', 'int', 'optional', 1),
            ScriptOption2(['-g', '--gpuID'], 'GPU index to run the program on.', 'int', 'optional')])

    options = parse_script_options2(sys.argv[1:], helper)

    filepath, output_folder, output_name, spacing, original_spacing, binning, modify_structure, solvent_correction, \
        solvent_density, ctf_correction, defocus, amplitude_contrast, voltage, Cs, sigma_decay, phase_flipped, \
        cut_zero, display_ctf, resolution, box_size, invert, mirror, cores, gpuID = options

    if resolution is None:
        resolution = 2 * spacing * binning

    if cut_zero is not None:
        zero_cut = 0
    else:
        zero_cut = -1

    ext = os.path.splitext(filepath)[1]
    if ext == '.pdb' or ext == '.cif':
        template = generate_template_from_pdb(filepath, spacing,
                                     binning=binning,
                                     modify_structure=modify_structure,
                                     solvent_correction=solvent_correction,
                                     solvent_density=solvent_density,
                                     apply_ctf_correction=ctf_correction,
                                     defocus=defocus * 1e-6,
                                     amplitude_contrast=amplitude_contrast,
                                     voltage=voltage * 1e3,
                                     Cs=Cs * 1e-3,
                                     ctf_decay=sigma_decay,
                                     phase_flip=phase_flipped,
                                     zero_cut=zero_cut,
                                     display_ctf=display_ctf,
                                     resolution=resolution,
                                     box_size=box_size,
                                     cores=cores,
                                     gpu_id=gpuID)

    elif ext == '.mrc' or ext == '.em':
        template = generate_template_from_map(filepath, spacing,
                                              original_spacing=original_spacing,
                                              binning=binning,
                                              apply_ctf_correction=ctf_correction,
                                              defocus=defocus * 1e-6,
                                              amplitude_contrast=amplitude_contrast,
                                              voltage=voltage * 1e3,
                                              Cs=Cs * 1e-3,
                                              ctf_decay=sigma_decay,
                                              phase_flip=phase_flipped,
                                              zero_cut=zero_cut,
                                              display_ctf=display_ctf,
                                              resolution=resolution,
                                              box_size=box_size)

    else:
        print('Invalid input file provided, should be either pdb, cif, mrc, or em.')
        sys.exit(0)

    if invert:
        template *= -1

    # output structure
    _, file = os.path.split(filepath)
    id, _ = os.path.splitext(file)
    if output_name:
        output_filepath = os.path.join(output_folder, output_name)
    else:
        output_filepath = os.path.join(output_folder, f'template_{id}_{spacing*binning:.2f}A_{template.shape[0]}px.mrc')

    print(f'Writing template as {output_filepath}')
    write(output_filepath, template)

    from pytom.gpu.initialize import xp, device

    if mirror:
        output_filepath_mirror = os.path.splitext(output_filepath)[0] + '_mirror' + os.path.splitext(output_filepath)[1]
        print(f'Writing template as {output_filepath_mirror}')
        write(output_filepath_mirror, xp.flip(template, 0))
