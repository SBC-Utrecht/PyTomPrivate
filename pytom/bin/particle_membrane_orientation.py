#!/usr/bin/env pytom

"""
Determine orientation of membrane-bound particles in relation to the segmentation map of a membrane.

Author: Marten Chaillet
"""
import numpy as np
import sys
import os
from scipy.spatial import distance
from pytom.simulation.membrane import Vector
from pytom.basic.structures import ParticleList


def convert_to_mesh(volume, cutoff=0.2, mesh_detail=2, display=False):
    from skimage import measure

    verts, faces, normals, values = measure.marching_cubes(volume, level=cutoff, step_size=mesh_detail)

    if display:
        try:
            import matplotlib.pyplot as plt
        except Exception as e:
            import matplotlib
            matplotlib.use('Qt5Agg')
            import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        # Display resulting triangular mesh using Matplotlib. This can also be done
        # with mayavi (see skimage.measure.marching_cubes_lewiner docstring).
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Fancy indexing: `verts[faces]` to generate a collection of triangles
        mesh = Poly3DCollection(verts[faces])
        mesh.set_edgecolor('k')
        ax.add_collection3d(mesh)

        ax.set_xlabel(f"x-axis: {volume.shape[0]}")
        ax.set_ylabel(f"y-axis: {volume.shape[1]}")
        ax.set_zlabel(f"z-axis: {volume.shape[2]}")

        ax.set_xlim(0, volume.shape[0])
        ax.set_ylim(0, volume.shape[1])
        ax.set_zlim(0, volume.shape[2])

        plt.tight_layout()
        plt.show()

    # invert normals because they point inwards
    return verts, faces, normals, values


def visualize(volume, cutoff, mesh_detail, particle_vecs, membrane_vecs, precalc=None):
    from skimage import measure
    if precalc is not None:
        verts, faces = precalc
    else:
        verts, faces, normals, values = measure.marching_cubes(volume, level=cutoff, step_size=mesh_detail)

    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        import matplotlib
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    # Display resulting triangular mesh using Matplotlib. This can also be done
    # with mayavi (see skimage.measure.marching_cubes_lewiner docstring).
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh = Poly3DCollection(verts[faces])
    mesh.set_edgecolor('k')
    ax.add_collection3d(mesh)

    pdata = np.array([list(c) + list(p) for c, p in particle_vecs])
    X, Y, Z, U, V, W = zip(*pdata)
    ax.quiver(X, Y, Z, [10 * u for u in U], [10 * u for u in V], [10 * u for u in W], color='red')

    mdata = np.array([list(c) + list(p) for c, p in membrane_vecs])
    X, Y, Z, U, V, W = zip(*mdata)
    ax.quiver(X, Y, Z, [10 * u for u in U], [10 * u for u in V], [10 * u for u in W], color='blue')

    ax.set_xlabel(f"x-axis: {volume.shape[0]}")
    ax.set_ylabel(f"y-axis: {volume.shape[1]}")
    ax.set_zlabel(f"z-axis: {volume.shape[2]}")

    ax.set_xlim(0, volume.shape[0])
    ax.set_ylim(0, volume.shape[1])
    ax.set_zlim(0, volume.shape[2])

    plt.tight_layout()
    plt.show()


def point_3d_in_triangle(point, v1, v2, v3):
    """

    Reference: W. Heidrich, Journal of Graphics, GPU, and Game Tools,Volume 10, Issue 3, 2005
    @param point:
    @type point: np.array(3)
    @param v1:
    @type v1: np.array(3)
    @param v2:
    @type v2: np.array(3)
    @param v3:
    @type v3: np.array(3)
    @return:
    @rtype: Bool
    """
    type_list = [type(point), type(v1), type(v2), type(v3)]
    assert len(set(type_list)) == 1, "Input is not of same type."
    if all([t is list for t in type_list]):
        point, v1, v2, v3 = map(np.array, [point, v1, v2, v3])

    uv = v2 - v1
    vv = v3 - v1
    nv = np.cross(uv, vv)
    ov = point - v1

    gamma = np.dot(np.cross(uv, ov), nv) / np.dot(nv, nv)
    beta = np.dot(np.cross(ov, vv), nv) / np.dot(nv, nv)
    alpha = 1 - gamma - beta

    projected_point = alpha * v1 + beta * v2 + gamma * v3
    is_in_triangle = (0 <= alpha <= 1) and (0 <= beta <= 1) and (0 <= gamma <= 1)

    return projected_point, is_in_triangle


def point_3d_in_line(point, a, b):
    """
    @param point:
    @type point: np.array(3)
    @param a: first point of line segment
    @type a: np.array(3)
    @param b: second point of line segment
    @type b: np.array(3)
    @rtype: Bool
    """
    ap = point - a
    ab = b - a
    projected_point = a + np.dot(ap, ab) / np.dot(ab, ab) * ab
    if distance.euclidean(a, b) == distance.euclidean(a, projected_point) + distance.euclidean(b, projected_point):
        return projected_point, True
    else:
        return projected_point, False


def faces_to_edges(faces, vert):
    index_dict = {0: [1, 2], 1: [0, 2], 2: [0, 1]}
    edges = []
    for face in faces:
        id = face.tolist().index(vert)
        rem = index_dict[id]
        e1 = [face[id], face[rem[0]]]
        e2 = [face[id], face[rem[1]]]
        if not any([True for e in edges if (e == e1 or e == e1[::-1])]):
            edges.append(e1)
        if not any([True for e in edges if (e == e2 or e == e2[::-1])]):
            edges.append(e2)
    return np.array(edges)


def get_edge_vector(edge, faces, verts):
    ids = (np.any(faces == edge[0], axis=1) == np.any(faces == edge[1], axis=1))
    faces_duo = faces[ids]
    assert len(faces_duo) == 2, "something is wrong in finding edge vector, could not find two connected faces"
    # find first normal
    v1, v2, v3 = verts[faces[0][0]], verts[faces[0][1]], verts[faces[0][2]]
    normal1 = Vector(np.cross(v2 - v1, v3 - v1))
    normal1.normalize()
    normal1.inverse()

    # find second normal
    v1, v2, v3 = verts[faces[1][0]], verts[faces[1][1]], verts[faces[1][2]]
    normal2 = Vector(np.cross(v2 - v1, v3 - v1))
    normal2.normalize()
    normal2.inverse()

    # calculate average
    average = normal1.average(normal2)
    average.normalize()
    return average


def find_variation(triangle_normals, vert_normal):
    diff_angles = []
    for n in triangle_normals:
        diff_angles.append(vert_normal.angle(n, degrees=True))
    diff_angles = np.array(diff_angles)
    return diff_angles.mean(), diff_angles.std(axis=0)


def find_orientations(plist, segmentation, cutoff, mesh_detail, reference_normal, verbose=False,
                      improved_precision=True):
    from pytom.voltools.utils import rotation_matrix

    # get the triangular mesh
    verts, faces, normals, values = convert_to_mesh(segmentation, cutoff, mesh_detail)

    # load reference unit vector and make sure it is actually a 'unit' vector
    unit_vector = Vector(reference_normal)
    unit_vector.normalize()

    distances, orientations, angle_stds = [], [], []
    particle_arrows, membrane_arrows = [], []
    exclusion_count = 0
    for p in plist:
        # get rotation matrix and convert to axis-angle
        rotation = p.getRotation().toVector()  # z1, z2, x
        matrix = rotation_matrix(rotation=(rotation[0], rotation[2], rotation[1]), rotation_order='rzxz')

        # rotate unit vector by rotation matrix...
        particle_normal = Vector(unit_vector.get())  # copy reference normal
        particle_normal.rotate(matrix[:3, :3])

        # get the coordinates of the particle
        coordinates = np.array(p.getPickPosition().toVector())

        # distance to each vertex in the triangle mesh
        distance_vector = np.sqrt(np.sum(np.subtract(verts, coordinates) ** 2, axis=1))
        # min_distance = np.min(distance)  # this is no longer needed?
        min_distance_idx = np.argmin(distance_vector)

        # try to find if the particle is directly above a triangle of the closest point
        face_idx = np.any(faces == min_distance_idx, axis=1)
        faces_subset = faces[face_idx]
        triangle_projection_points, point_on_triangle_surface, triangle_normals = [], [], []
        for face in faces_subset:
            # get triangle vertices
            v1, v2, v3 = verts[face[0]], verts[face[1]], verts[face[2]]

            # get projected point and bool telling whether its in the triangle
            triangle_projection, lies_in_triangle = point_3d_in_triangle(coordinates, v1, v2, v3)
            triangle_projection_points.append(triangle_projection)
            point_on_triangle_surface.append(lies_in_triangle)

            # calculate normal and append it
            triangle_norm = Vector(np.cross(v2 - v1, v3 - v1))
            triangle_norm.normalize()
            triangle_norm.inverse()
            triangle_normals.append(triangle_norm)

        # visualize(segmentation, cutoff, mesh_detail, [(c, v.get()) for c, v in zip(triangle_projection_points,
        #                                                                       triangle_normals)],
        #           [(coordinates, particle_normal.get())], precalc=(verts, faces_subset))

        # compare vertex normal with triangle normals
        # visualize(segmentation, cutoff, mesh_detail, [(c, v.get()) for c, v in zip(triangle_projection_points,
        #                                                                            triangle_normals)],
        #           [(verts[min_distance_idx], normals[min_distance_idx])], precalc=(verts, faces_subset))

        mean, std = find_variation(triangle_normals, particle_normal)
        if verbose:
            print('average and standard deviation of membrane normals: ', mean, std)

        # select triangle if applicable, otherwise select edge
        if improved_precision:
            if sum(point_on_triangle_surface) == 1:  # this is the easy case, just above a single triangle
                if verbose: print('particle above triangle')
                id = point_on_triangle_surface.index(True)
                membrane_normal, membrane_point = triangle_normals[id], triangle_projection_points[id]
            elif sum(point_on_triangle_surface) > 1:  # above two or more, select triangle with shortest distance
                if verbose: print('particle above two or more triangles')
                ids = [i for i, t in enumerate(point_on_triangle_surface) if t]
                dists = [distance.euclidean(triangle_projection_points[i], coordinates) for i in ids]
                id = ids[dists.index(min(dists))]
                membrane_normal, membrane_point = triangle_normals[id], triangle_projection_points[id]
            else:  # not above any triangle, it is either above an edge or above a point
                # first select all possible edges
                edges = faces_to_edges(faces_subset, min_distance_idx)
                line_projection_points, point_on_line, line_normals = [], [], []
                for edge in edges:
                    line_projection, lies_on_line = point_3d_in_line(coordinates, verts[edge[0]], verts[edge[1]])
                    line_projection_points.append(line_projection)
                    point_on_line.append(lies_on_line)
                    line_normals.append(get_edge_vector(edge, faces_subset, verts))

                # select edge if applicable, otherwise select vertex
                if sum(point_on_line) == 1:
                    if verbose: print('particle above edge')
                    id = point_on_line.index(True)
                    membrane_normal, membrane_point = line_normals[id], line_projection_points[id]
                elif sum(point_on_line) > 1:
                    if verbose: print('particle above two or more edges')
                    ids = [i for i, t in enumerate(point_on_line)]
                    dists = [distance.euclidean(line_projection_points[i], coordinates) for i in ids]
                    id = ids[dists.index(min(dists))]
                    membrane_normal, membrane_point = line_normals[id], line_projection_points[id]
                else:  # finally, if not above anything else, select the vertex
                    if verbose: print('particle above vertex')
                    membrane_normal = Vector(normals[min_distance_idx])
                    membrane_point = verts[min_distance_idx]
        else:
            membrane_normal = Vector(np.array([n.get() for n in triangle_normals]).sum(axis=0) / len(triangle_normals))
            membrane_point = verts[min_distance_idx]

        # get the difference angle
        difference_angle = particle_normal.angle(membrane_normal, degrees=True)
        orientations.append(difference_angle)

        # get the distance
        distances.append(distance.euclidean(membrane_point, coordinates))

        # append the variation
        angle_stds.append(std)

        # create arrows for bild file
        particle_arrows.append((coordinates, particle_normal.get()))
        membrane_arrows.append((membrane_point, membrane_normal.get()))

    return distances, orientations, angle_stds, particle_arrows, membrane_arrows


def write_arrow_bild(p_arrows, m_arrows, filename, outlier_filter=None):
    if outlier_filter is None:
        outlier_filter = [True] * len(p_arrows)
    with open(filename, 'w') as stream:
        for (arrow, normal), inlier in zip(p_arrows, outlier_filter):
            if inlier:
                stream.write('.color red\n')
                stream.write(f'.arrow {arrow[0]:.2f} {arrow[1]:.2f} {arrow[2]:.2f} '
                             f'{arrow[0] + 15 * normal[0]:.2f} '
                             f'{arrow[1] + 15 * normal[1]:.2f} '
                             f'{arrow[2] + 15 * normal[2]:.2f}\n')
        for (arrow, normal), inlier in zip(m_arrows, outlier_filter):
            if inlier:
                stream.write('.color blue\n')
                stream.write(f'.arrow {arrow[0]:.2f} {arrow[1]:.2f} {arrow[2]:.2f} '
                             f'{arrow[0] + 15 * normal[0]:.2f} '
                             f'{arrow[1] + 15 * normal[1]:.2f} '
                             f'{arrow[2] + 15 * normal[2]:.2f}\n')


def write_split_particle_lists(orientations_per_particle_list, distance_range, angle_outlier_filter,
                               split_lists_angle, voxel_size, output_name):

    for data in orientations_per_particle_list:
        d_filter = np.logical_and(np.array(data['distances']) * voxel_size >= distance_range[0],
                                  np.array(data['distances']) * voxel_size <= distance_range[1])
        o_filter = (np.array(data['stds']) < angle_outlier_filter)

        c_filter = np.logical_and(d_filter, o_filter)
        set_1_filter = np.logical_and(np.array(data['orientations']) <= split_lists_angle, c_filter)  # include cutoff
        set_2_filter = np.logical_and(np.array(data['orientations']) > split_lists_angle, c_filter)  # above without

        # apply filter to particle list and write to correct destination
        set_1, set_2 = ParticleList(), ParticleList()
        for s1, s2, part in zip(set_1_filter, set_2_filter, data['particle_list']):
            if s1:
                set_1.append(part)
            elif s2:
                set_2.append(part)

        print('particles in set1 and set2 after splitting angles: ', len(set_1), len(set_2))

        # output file names
        set_1_file_name = os.path.join(output_name,
                                       os.path.splitext(os.path.split(data['segmentation_file'])[1])[0] + \
                                       f'_<{split_lists_angle}deg.xml')
        set_2_file_name = os.path.join(output_name,
                                       os.path.splitext(os.path.split(data['segmentation_file'])[1])[0] + \
                                       f'_>{split_lists_angle}deg.xml')
        # write the sets
        set_1.toXMLFile(set_1_file_name)
        set_2.toXMLFile(set_2_file_name)


if __name__ == '__main__':
    from pytom.tools.script_helper import ScriptHelper2, ScriptOption2
    from pytom.tools.parse_script_options import parse_script_options2
    from pytom.agnostic.io import read

    helper = ScriptHelper2(
        sys.argv[0].split('/')[-1],  # script name
        description='Find the orientations of particles with respect to a membrane. The particles should be in the '
                    'pytom particle list (xml) format. Membrane volume should be a segmented map of the membrane in '
                    'the tomogram, take care sizes correspond.',
        authors='Marten Chaillet',
        options=[
            ScriptOption2(['-i', '--input_file'], 'Particle list (.xml) or linker file (.txt) listing matching '
                                                  'segmentation and particle lists.', 'file', 'required'),
            ScriptOption2(['-s', '--segmentation_file'], 'Membrane map in mrc/em/rec format.', 'file',
                          'optional'),
            ScriptOption2(['-o', '--output_name'], 'Folder where the output report is stored',
                          'string', 'optional'),
            ScriptOption2(['-v', '--voxel_size'], 'Voxel size of segmentation model.',
                          'float', 'optional'),
            ScriptOption2(['-c', '--cutoff'], 'Cutoff value for converting membrane model to a triangular mesh.',
                          'float', 'optional', 0.2),
            ScriptOption2(['-m', '--mesh_detail'], 'Detail of the mesh, i.e. how fine it should be sampled from the '
                                                   'volume.', 'int', 'optional', 2),
            ScriptOption2(['-u', '--template_unit_vector'], 'Direction of unit vector of template that the '
                                                            'orientations are relative to. It will not affect the '
                                                            'distribution of populations that you find, only their '
                                                            'angular difference with the template.',
                          'float,float,float', 'optional', [.0, .0, 1.]),
            ScriptOption2(['-f', '--filter_nstd'], 'Number of standard deviations to filter outliers based on '
                                                   'distance to membrane mesh and variation of angles under the '
                                                   'ribosome.', 'int',
                          'optional', 3),
            ScriptOption2(['-r', '--distance_range'], 'Min and max distance from closest membrane point to consider, '
                                                      'min distance can for example be set to the radius of the '
                                                      'particle (if the current coordinate is in the center',
                          'float,float', 'optional'),
            ScriptOption2(['-a', '--angle_std_cutoff'], 'Set a threshold for removing orientations that are likely '
                                                        'imprecise due to inaccurate patches of membrane. Cutoff is '
                                                        'relative to std of triangles directly under membrane. So a '
                                                        'value of 15 degrees implies the std cannot be larger than 15 '
                                                        'for the membrane patch.',
                          'float', 'optional'),
            ScriptOption2(['--axis_limit_distance'], 'Specify up to what value distances should be plotted. Distances '
                                                     'of particles can become huge, which decreases plot '
                                                     'visibility. Default 500A. ', 'float', 'optional', 500),
            ScriptOption2(['-n', '--bins'], 'Number of bins for histogram.',
                          'int', 'optional', 30),
            ScriptOption2(['--write_split_particle_lists'], 'Write .xml particle lists split on specified angle.',
                          'float', 'optional'),
            ScriptOption2(['--verbose'], 'Be verbose.', 'no arguments', 'optional')])

    options = parse_script_options2(sys.argv[1:], helper)

    input_file, segmentation_file, output_name, voxel_size, cutoff, mesh_detail, template_normal, nstd, \
        distance_range, angle_std_cutoff, axis_distance_limit, nbins, split_lists_angle, verbose = options

    if output_name is not None and not os.path.exists(output_name):
        os.mkdir(output_name)

    _, ext = os.path.splitext(input_file)
    if ext == '.xml':
        segmentation = read(segmentation_file)
        particle_list = ParticleList()
        particle_list.fromXMLFile(input_file)
        distances, orientations, stds, p_arrows, m_arrows = find_orientations(particle_list, segmentation, cutoff,
                                                                        mesh_detail,
                                                                        template_normal, verbose=verbose)
        # create dict with information per list
        orientations_per_particle_list = [{'particle_list': particle_list,
                                                   'segmentation_file': segmentation_file,
                                                   'particle_list_file': input_file,
                                                   'distances': distances,
                                                   'orientations': orientations,
                                                   'stds': stds}]
    elif ext == '.txt':
        linker_data = np.genfromtxt(input_file, dtype=[('segmentation', 'U1000'), ('particle_list', 'U1000')],
                                    skip_header=1)
        segmentation_files = linker_data['segmentation']
        particle_list_files = linker_data['particle_list']
        distances, orientations, stds = [], [], []
        orientations_per_particle_list = []
        for s, p in zip(segmentation_files, particle_list_files):
            print(f'Run on segmentation {s} with plist {p}')
            name = os.path.splitext(os.path.split(s)[1])[0]
            segmentation = read(s)
            particle_list = ParticleList()
            particle_list.fromXMLFile(p)
            d, o, std, p_arrows, m_arrows = find_orientations(particle_list, segmentation, cutoff, mesh_detail,
                                                         template_normal, verbose=verbose)
            distances += d
            orientations += o
            stds += std
            write_arrow_bild(p_arrows, m_arrows, name + '.bild')

            # add a dict that stores information per list
            orientations_per_particle_list.append({'particle_list': particle_list,
                                                   'segmentation_file': s,
                                                   'particle_list_file': p,
                                                   'distances': d,
                                                   'orientations': o,
                                                   'stds': std})
    else:
        print('Invalid input extension.')
        sys.exit(0)

    # total number of particles
    n_pre_filter = len(orientations)

    # conver to numpy arrays for calculation
    distances = np.array(distances)
    distances *= voxel_size if voxel_size is not None else 1
    orientations = np.array(orientations)
    stds = np.array(stds)

    # find distance cutoff  and angle cutoff if not specified
    distance_range = (distances.mean() - nstd * distances.std(), distances.mean() + nstd * distances.std()) if \
        distance_range is None else distance_range
    angle_std_cutoff = nstd * stds.std() + stds.mean() if angle_std_cutoff is None else angle_std_cutoff

    if split_lists_angle is not None and output_name is not None:
        write_split_particle_lists(orientations_per_particle_list, distance_range, angle_std_cutoff,
                                   split_lists_angle, voxel_size, output_name)

    # Filter the combined dataset to get better statistics of distributions for plot
    distance_outlier_filter = np.logical_and(distances >= distance_range[0], distances <= distance_range[1])
    angle_outlier_filter = (stds < angle_std_cutoff)
    outlier_filter = np.logical_and(distance_outlier_filter, angle_outlier_filter)
    distances_filtered = distances[outlier_filter]
    orientations = orientations[outlier_filter]
    n = len(orientations)

    # output number of particle excluded due to filtering
    total_outliers = (outlier_filter == False).sum()
    print(f'Filtered {total_outliers} from total of {n_pre_filter} particles')

    # ============================== create plots
    try:
        import matplotlib.pyplot as plt
    except Exception as e:
        import matplotlib
        matplotlib.use('Qt5Agg')
        import matplotlib.pyplot as plt

    font = {'size': 16}
    plt.rc('font', **font)

    # do some plotting
    fig, ax = plt.subplots(figsize=(5, 5))

    ax.hist(orientations, bins=nbins, color='black', histtype='stepfilled', alpha=0.5)
    ax.hist(orientations, bins=nbins, color='black', histtype='step')
    ax.set_ylabel('Number of particles')
    ax.set_xlabel('Angle (degrees)')

    # information about number of particles
    textstr = f'N={n}'
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
               verticalalignment='top', bbox=props)

    plt.tight_layout()

    if output_name is not None:
        plt.savefig(os.path.join(output_name, 'orientation_distribution.png'), dpi=300, format='png', transparent=True)
        print(f"wrote {os.path.join(output_name, 'orientation_distribution.png')}")
        if ext == '.xml':
            write_arrow_bild(p_arrows, m_arrows, os.path.join(output_name, 'vectors.bild'),
                             outlier_filter=outlier_filter)
            print(f"wrote {os.path.join(output_name, 'vectors.bild')}")
    else:
        plt.show()

    # ==================== create plot of distance vs. angles

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(distances_filtered, orientations, alpha=0.7, color='gray')
    if voxel_size is not None:
        ax.set_xlabel(r'Distance ($\AA$)')
    else:
        ax.set_xlabel('Distance (voxels)')
    ax.set_ylabel('Angle (degrees)')

    # information about number of particles
    textstr = f'N={n}'
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
               verticalalignment='top', bbox=props)

    plt.tight_layout()

    if output_name is not None:
        plt.savefig(os.path.join(output_name, 'scatter.png'), dpi=300, format='png', transparent=True)
        print(f"wrote {os.path.join(output_name, 'scatter.png')}")
    else:
        plt.show()

    # ==================== plot histogram of distances with cutoffs

    fig, ax = plt.subplots(figsize=(5, 5))
    bins, _, _ = ax.hist(distances[distances < axis_distance_limit],
                            bins=nbins, color='black', histtype='stepfilled', alpha=0.5)
    bins, _, _ = ax.hist(distances[distances < axis_distance_limit],
                            bins=nbins, color='black', histtype='step')
    if voxel_size is not None:
        ax.set_xlabel(r'Distance ($\AA$)')
    else:
        ax.set_xlabel('Distance (voxels)')
    ax.set_ylabel('Number of particles')
    ax.vlines(distance_range[0], 0, max(bins), label='cutoff', color='orange', linestyles='dashed', linewidth=3)
    ax.vlines(distance_range[1], 0, max(bins), color='orange', linestyles='dashed', linewidth=3)
    ax.legend()

    plt.tight_layout()

    if output_name is not None:
        plt.savefig(os.path.join(output_name, 'distances.png'), dpi=300, format='png', transparent=True)
        print(f"wrote {os.path.join(output_name, 'distances.png')}")
    else:
        plt.show()

    # ==================== plot histogram of angle stds with cutoffs

    fig, ax = plt.subplots(figsize=(5, 5))
    bins, _, _ = ax.hist(stds, bins=nbins, color='black', histtype='stepfilled', alpha=0.5)
    bins, _, _ = ax.hist(stds, bins=nbins, color='black', histtype='step')
    if angle_std_cutoff is not None:
        ax.vlines(angle_std_cutoff, 0, max(bins), label='cutoff', color='orange', linestyles='dashed', linewidth=3)
        ax.legend()
    ax.set_xlabel('Std (degrees)')
    ax.set_ylabel('Number of particles')

    # information about number of particles
    textstr = f'N={n_pre_filter}'
    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
    # place a text box in upper left in axes coords
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)

    plt.tight_layout()

    if output_name is not None:
        plt.savefig(os.path.join(output_name, 'std.png'), dpi=300, format='png', transparent=True)
        print(f"wrote {os.path.join(output_name, 'std.png')}")
    else:
        plt.show()
