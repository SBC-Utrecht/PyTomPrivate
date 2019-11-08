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
    from time import gmtime, strftime
    from pytom.reconstruction.reconstruct_local_alignment import polish_particles

    helper = ScriptHelper(sys.argv[0].split('/')[-1],  # script name
                          description='Reconstruct a local alignment of particles based on a global alignment.',
                          authors='Douwe Schulte',
                          options=[ScriptOption(['-d', '--data'],
                                                'A list (delimited by semicolons) which '
                                                'consists of a path to a particlelist followed by the path to the '
                                                'projection directory. Ex path/to/tomo_000.xml;path/to/projection;path/'
                                                'to/tomo_001.xml;path/to/projection', 'has arguments', 'required'),
                                   ScriptOption(['-s', '--particleSize'],
                                                'Output particle size.', 'uint', 'required'),
                                   ScriptOption(['-b', '--binning'],
                                                'Binning factor of the particle list.', 'uint', 'required'),
                                   ScriptOption(['-o', '--cuttingOffset'],
                                                'Cutting offset of the particle list (int x, int y, int z).',
                                                'int,int,int', 'required'),
                                   ScriptOption(['-a', '--averagedSubtomogram'],
                                                'The path to an averaged subtomogram, to use instead of many '
                                                'subtomograms', 'string', 'required'),
                                   ScriptOption(['-g', '--createGraphics'],
                                                'Flag to turn on to create graphical reports of intermediate steps of '
                                                'the particle polishing', 'no arguments', 'optional', False),
                                   ScriptOption(['-n', '--numberOfParticles'],
                                                'To take a subset of the particlelist for debugging purposes',
                                                'uint', 'optional', -1),
                                   ScriptOption(['--skipAlignment'],
                                                'Skips the alignment/particle polish phase, only does the '
                                                'reconstruction and FRM alignment.', 'no arguments', 'optional', False),
                                   ScriptOption(['-f', '--FSCPath'], "The path to an FSC file (.dat) to use as a filter"
                                                " for the cutouts.", 'has arguments', 'optional', ''),
                                   ScriptOption(['--GJobName'], 'The name of the GLocal job', 'string',
                                                'optional', strftime("glocaljob-%D-%m-%Y", gmtime())),
                                   ScriptOption(['--GNodes'], 'The amount of nodes GLocal can use', 'uint', 'optional', 5),
                                   ScriptOption(['--Gparticlelist'], 'The particlelist to be used by GLocal', 'string', 'optional'),
                                   ScriptOption(['--dimZ'], 'The dimension on the Z axis unbinned, default is the same as dimension X', 'uint', 'optional'),
                                   ScriptOption(['-std', '--shiftStandardDeviation'],
                                                'Used to constrain the maximum shift found by the peak finding'
                                                ' algorithm, defined as the standard deviation excepted', 'uint', 'optional', 5),
                                   ScriptOption(['--peakBorderSize'],
                                                'Used to constrain the maximum shift found by the peak finding'
                                                ' algorithm, defined as the maximum number of standard deviations'
                                                ' a peak can be found from the center.', 'uint', 'optional', 5)
                                   ])

    proj_dir, vol_size, binning, offset, averaged_subtomogram, create_graphics, number_of_particles, skip_alignment, \
    fsc_path, glocal_jobname, glocal_nodes, glocal_particlelist, dimz, std, std_num = parse_script_options(sys.argv[1:], helper)

    if vol_size % 2 != 0:
        raise ValueError("The particle size has to be an even number.")

    peak_border = vol_size / 2 - std * std_num
    if peak_border < 0:
        raise ValueError("The given shiftStandardDeviation and peakBorderSize result in a maximal shift "
                         "bigger than the given volume. Please use a bigger volume or a smaller maximal shift.")

    names = []
    if ':' in proj_dir:
        temp = proj_dir.split(':')
        if len(temp) % 2 == 1:
            raise Exception("The data given is not valid (invalid number of items, should be even)")
        for i in range(len(temp) / 2):
            names.append((temp[i*2], temp[i*2+1]))
    else:
        raise Exception("The data given is not valid (invalid number of items, should be at least two)")

    print("Parsed arguments")

    from pytom.tompy.mpi import MPI

    mpi = MPI()
    mpi.begin()

    for i, n in enumerate(names):
        print("Running " + n[0])
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

        print(tilt_angles)

        polish_particles(proj, vol_size, binning, offset, tilt_angles, n[0], proj_dir, mpi,
                        averaged_subtomogram, number_of_particles, skip_alignment, (True if i == 3 else False),
                        fsc_path, glocal_jobname, glocal_nodes, glocal_particlelist, dimz, peak_border, create_graphics)
        print("Finished "+n[0])

    mpi.end()
