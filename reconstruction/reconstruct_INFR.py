#!/usr/bin/env python

'''
Created on Jun 6, 2013

@author: yuxiangchen
'''

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    helper = ScriptHelper(sys.argv[0].split('/')[-1], # script name
                          description='Reconstruct the whole tomogram from projections using iterative NFFT method.',
                          authors='Yuxiang Chen',
                          options=[ScriptOption('-d', 'Unbinned projection directory. Note the projections should be aligned but not weighted!', 'string', 'required'),
                                   ScriptOption('-o', 'Output filename.', 'string', 'required'),
                                   ScriptOption('-i', 'Number of iterations to run.', 'uint', 'optional', 10)])

    proj_dir, output_filename, iter = parse_script_options(sys.argv[1:], helper)

    # start reconstruction
    from tompy.io import read, write
    from nufft.reconstruction import fourier_2d1d_iter_reconstruct
    from pytom.reconstruction.reconstructionStructures import ProjectionList
    projections = ProjectionList()
    projections.loadDirectory(proj_dir)
    projections.sort()
    
    projs = []
    tilt_angles = []
    for p in projections:
        projs.append(read(p.getFilename()))
        tilt_angles.append(p.getTiltAngle())
    
    v = fourier_2d1d_iter_reconstruct(projs, tilt_angles, iter)
    write(output_filename, v)

