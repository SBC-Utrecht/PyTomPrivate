#!/usr/bin/env pytom

import sys
import os
import multiprocessing
from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
from pytom.tools.parse_script_options import parse_script_options2
from pytom.agnostic.io import read, write, read_pixelsize
from pytom.basic.structures import ParticleList
from pytom.agnostic.transform import fourier_reduced2full, fourier_full2reduced
from functools import partial


def convert_subtomos(pl, invert_contrast=False):
    for i, p in enumerate(pl):
        if i > 0 and i % 1000 == 0:
            print(i, ' done')

        if invert_contrast:  # invert contrast if needed
            subtomo_file = p.getFilename()
            folder, filename = os.path.split(subtomo_file)
            basename, ext = os.path.splitext(filename)
            new_subtomo_file = os.path.join(folder, basename + '_pytom' + ext)
            pixel_size = read_pixelsize(subtomo_file)
            subtomo = read(subtomo_file)
            write(new_subtomo_file, -1 * subtomo, pixel_size=pixel_size)
            # p.setFilename(new_subtomo_file)

        # change fourier reduced dimension of 3d ctf
        ctf_file = p.getWedge().getFilename()
        folder, filename = os.path.split(ctf_file)
        basename, ext = os.path.splitext(filename)
        new_ctf_file = os.path.join(folder, basename + '_pytom' + ext)
        pixel_size = read_pixelsize(ctf_file)
        ctf = read(ctf_file)
        # ctf[ctf < 0] = 0
        # ctf[ctf > 1] = 1
        full_ctf = fourier_reduced2full(ctf, reduced_axis=0, isodd=ctf.shape[1] % 2)
        red_ctf = fourier_full2reduced(full_ctf, reduced_axis=2)
        write(new_ctf_file, red_ctf, pixel_size=pixel_size)
        # p.getWedge.setFilename(new_ctf_file)


def import_subtomos(star_file, xml_file, invert_contrast=False, n_cores=8):
    # read star only handles relion3.0/3.1 currently
    target, file_name = os.path.split(xml_file)
    os.system(f"convert.py -f {star_file} -o xml -t "
              f"{target if target != '' else './'} --outname {file_name.strip('.xml')}")
    particle_list = ParticleList()
    particle_list.fromXMLFile(xml_file)

    print('each process will print some output about its progress')
    print('total number of particles to loop over: ', len(particle_list))
    # map cores to convert the files onto the split particle list
    with multiprocessing.Pool(processes=n_cores) as pool:
        pool.map(partial(convert_subtomos, invert_contrast=invert_contrast), particle_list.splitNSublists(n_cores))

    # update filenames in unsplit particle list
    for p in particle_list:
        if invert_contrast:  # invert contrast if needed
            subtomo_file = p.getFilename()
            folder, filename = os.path.split(subtomo_file)
            basename, ext = os.path.splitext(filename)
            new_subtomo_file = os.path.join(folder, basename + '_pytom' + ext)
            p.setFilename(new_subtomo_file)

        # change fourier reduced dimension of 3d ctf
        ctf_file = p.getWedge().getFilename()
        folder, filename = os.path.split(ctf_file)
        basename, ext = os.path.splitext(filename)
        new_ctf_file = os.path.join(folder, basename + '_pytom' + ext)
        p.getWedge().setFilename(new_ctf_file)

    particle_list.toXMLFile(xml_file)


if __name__ == '__main__':
    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Import subtomograms from other software. Currently supported: WarpM',
        authors='Marten Chaillet',
        options=[
            ScriptOption2(['-i', '--input-star-file'], 'Input particle list in .star format',
                          'file', 'required'),
            ScriptOption2(['-o', '--output-xml-file'], 'PyTom xml file that will store the particle information',
                          'string', 'optional'),
            ScriptOption2(['--invert'], 'Invert contrast of subtomograms. PyTom wants particles to be negative, '
                                        'while Relion wants positive particles.',
                          'no arguments', 'optional'),
            ScriptOption2(['-c', '--cores'], 'number of cpu cores to split processes on', 'int', 'optional', 1)])

    options = parse_script_options2(sys.argv[1:], helper)
    input_star, output_xml, invert_contrast_flag, n_cores = options

    # make boolean
    invert_contrast_flag = False if invert_contrast_flag is None else True

    # set some paths for the output
    if output_xml is None:
        folder, filename = os.path.split(input_star)
        file_basename = os.path.splitext(filename)[0]
        output_xml = os.path.join(folder, file_basename + '.xml')

    import_subtomos(input_star, output_xml, invert_contrast=invert_contrast_flag, n_cores=n_cores)
