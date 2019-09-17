#!/usr/bin/env pytom

"""
Created on Sep 02, 2019

@author: dschulte
"""

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.tompy.mpi import MPI

    helper = ScriptHelper(sys.argv[0].split('/')[-1],  # script name
                          description='Reconstruct a local alignment of particles based on a global alignment.',
                          authors='Douwe Schulte',
                          options=[ScriptOption(['-p', '--particleList'],
                                                'Particle list (filename).', True, False),
                                   ScriptOption(['-d', '--projectionDirectory'],
                                                'Unbinned projection directory. Note the projections should be aligned '
                                                'but not weighted! (path)', True, False),
                                   ScriptOption(['-s', '--particleSize'],
                                                'Output particle size (int).', True, False),
                                   ScriptOption(['-b', '--binning'],
                                                'Binning factor of the particle list (int).', True, False),
                                   ScriptOption(['-o', '--cuttingOffset'],
                                                'Cutting offset of the particle list (int x, int y, int z).',
                                                True, False),
                                   ScriptOption(['-a', '--averagedSubtomogram'],
                                                'The path to an averaged subtomogram, to use instead of many '
                                                'subtomograms', True, True),
                                   ScriptOption(['-i', '--INFRIterations'],
                                                'Number of iterations to run, when running INFR, using this option '
                                                'automatically sets the reconstruction method to INFR.', True, True),
                                   ScriptOption(['-m', '--reconstructionMethod'],
                                                'The reconstruction method for creating the initial subtomograms, has '
                                                'to be INFR or WBP. If this is not specified no initial subtomograms '
                                                'are created and subtomograms are expected to be made available in the '
                                                'normal folder.', True, True),
                                   ScriptOption(['-g', '--createGraphics'],
                                                'Flag to turn on to create graphical reports of intermediate steps of '
                                                'the particle polishing', False, True),
                                   ScriptOption(['-n', '--numberOfParticles'],
                                                'To take a subset of the particlelist for debugging purposes (int)',
                                                True, True),
                                   ScriptOption(['--skipAlignment'],
                                                'Skips the alignment/particle polish phase, only does the '
                                                'reconstruction and FRM alignment.', False, True),
                                   ScriptOption(['-h', '--help'],
                                                'Help.', False, True)])
    if len(sys.argv) == 1:
        print helper
        sys.exit()
    try:
        pl_filename, proj_dir, vol_size, binning, offset, averaged_subtomogram, infr_iter, reconstruction_method, \
         create_graphics, number_of_particles, skip_alignment, b_help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print e
        sys.exit()
    if b_help is True:
        print helper
        sys.exit()

    # parse the arguments
    vol_size = int(vol_size)
    binning = int(binning)

    offset = [int(i) for i in offset.split(',')]
    if len(offset) != 3:
        raise Exception("The offset should be defined with three values. Structure: x,y,z. "
                        "In which x, y and z are integers.")

    if reconstruction_method == "WBP" or reconstruction_method == "INFR":
        # fine go on
        create_subtomograms = True
    elif reconstruction_method == "" or not reconstruction_method:
        create_subtomograms = False
    else:
        raise Exception((reconstruction_method if reconstruction_method else "[nothing]") +
                        " is not a known reconstruction method, use WBP or INFR")

    if not averaged_subtomogram:
        averaged_subtomogram = False

    if create_graphics:
        create_graphics = True
    else:
        create_graphics = False

    if skip_alignment:
        skip_alignment = True
    else:
        skip_alignment = False

    if infr_iter:
        infr_iter = int(infr_iter)
        reconstruction_method = "INFR"
    else:
        infr_iter = -1

    # force the user to specify even-sized volume
    assert vol_size % 2 == 0

    # pass everything to the function in reconstruction/reconstruct_local_alignment.py
    from pytom.reconstruction.reconstruct_local_alignment import local_alignment

    if number_of_particles:
        number_of_particles = int(number_of_particles)
    else:
        number_of_particles = -1

    # Split particlelist based on filename
    # Map filename to right projection directory

    print("Parsed arguments")

    names = [("particleList_TM_tomogram_010_WBP", "/data2/dschulte/BachelorThesis/Data/VPP2/03_Tomographic_Reconstruction/tomogram_000/alignment/marker_0001_-60.0,60.0/"),
             ("particleList_TM_tomogram_011_WBP", "/data2/dschulte/BachelorThesis/Data/VPP2/03_Tomographic_Reconstruction/tomogram_001/alignment/marker_0001_-60.0,60.0/"),
             ("particleList_TM_tomogram_012_WBP", "/data2/dschulte/BachelorThesis/Data/VPP2/03_Tomographic_Reconstruction/tomogram_002/alignment/marker_0001_-60.0,60.0/"),
             ("particleList_TM_tomogram_013_WBP", "/data2/dschulte/BachelorThesis/Data/VPP2/03_Tomographic_Reconstruction/tomogram_003/alignment/marker_0001_-60.0,60.0/")]
    for i, n in enumerate(names):
        print("Running "+n[0])
        proj_dir = n[1]

        # load projection list
        from pytom.reconstruction.reconstructionStructures import ProjectionList

        projections = ProjectionList()
        projections.loadDirectory(proj_dir)
        projections.sort()

        # create the list of projectionfilenames and tiltangles
        proj = []
        tilt_angles = []
        for p in projections:
            proj.append(p.getFilename())
            tilt_angles.append(p.getTiltAngle())

        local_alignment(proj, vol_size, binning, offset, tilt_angles, pl_filename, proj_dir, reconstruction_method,
                        infr_iter, create_graphics, create_subtomograms, averaged_subtomogram, number_of_particles,
                        skip_alignment, n[0], 1 if n == 3 else 0)
        print("Finished "+n[0])
