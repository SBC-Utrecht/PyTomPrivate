#!/usr/bin/env pytom

"""
Determine orientation of membrane-bound particles in relation to the segmentation map of a membrane.

Author: Marten Chaillet
"""
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import sklearn.cluster as cluster
import hdbscan
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import distance
from pytom.simulation.membrane import Vector
from pytom.basic.structures import ParticleList
from skimage import measure
from pytom.voltools.utils import rotation_matrix
from pytom.angles.angleFnc import matToZXZ
from numba import njit, prange

epsilon = 0.0000005


# ==================================CONTAINED IN MEMBRANE MESH CLASS ===================================================
def convert_to_mesh(volume, cutoff=0.2, mesh_detail=2, display=False):

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


def find_membrane_surface(coordinates, verts, faces, normals):
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

        # get particle coordinate projected on triangle face and bool telling whether its in the triangle
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

    # TODO this is unnecessary

    std = find_std(triangle_normals)

    # consider situations of coordinate above triangle, line, or vertex
    if sum(point_on_triangle_surface) == 1:  # this is the easy case, just above a single triangle
        id = point_on_triangle_surface.index(True)
        membrane_normal, membrane_point = triangle_normals[id], triangle_projection_points[id]

    elif sum(point_on_triangle_surface) > 1:  # above two or more, select triangle with shortest distance
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
            id = point_on_line.index(True)
            membrane_normal, membrane_point = line_normals[id], line_projection_points[id]

        elif sum(point_on_line) > 1:
            ids = [i for i, t in enumerate(point_on_line)]
            dists = [distance.euclidean(line_projection_points[i], coordinates) for i in ids]
            id = ids[dists.index(min(dists))]
            membrane_normal, membrane_point = line_normals[id], line_projection_points[id]

        else:  # finally, if not above anything else, select the vertex
            membrane_normal = Vector(normals[min_distance_idx])
            membrane_point = verts[min_distance_idx]

    return membrane_point, membrane_normal, std


def find_orientations(plist, segmentation, cutoff, mesh_detail, reference_normal, verbose=False):
    from pytom.voltools.utils import rotation_matrix
    from pytom.tools.maths import Matrix
    from pytom.angles.angleFnc import matToZXZ, differenceAngleOfTwoRotations, matToAxisAngle

    # get the triangular mesh
    verts, faces, normals, values = convert_to_mesh(segmentation, cutoff, mesh_detail)

    # load reference unit vector and make sure it is actually a 'unit' vector
    unit_vector = Vector([0, 0, 1])
    unit_vector.normalize()

    # template rotation
    reference_basis = Vector(reference_normal)
    reference_basis.normalize()

    # basis rotation -> rotation of unit vector onto reference basis
    # rm_to_ref_bas = reference_basis.get_rotation(unit_vector)
    rm_to_ref_bas = unit_vector.get_rotation(reference_basis)

    # initialize some lists
    distances, orientations, orientations_around_axis, angle_stds, final_zxz = [], [], [], [], []
    particle_arrows, membrane_arrows = [], []
    # exclusion_count = 0
    for p in plist:

        # get the coordinates of the particle
        coordinates = np.array(p.getPickPosition().toVector())

        # find closest point on membrane
        membrane_point, membrane_normal, angular_variation = find_membrane_surface(coordinates,
                                                                                   verts, faces, normals)

        if verbose:
            print('average and standard deviation of membrane normals: ', angular_variation, std)

        # get rotation matrix and convert to axis-angle
        zxz_par_to_ref = p.getRotation().toVector(convention='zxz')  # z1, z2, x
        rot_ref_to_par = rotation_matrix(rotation=zxz_par_to_ref,
                                         rotation_order='rzxz')[:3, :3].T  # transpose for ref to par
        # par_rot_basis = np.dot(rm_to_ref_bas.T, par_rot)  # matrix multiply

        rot_mem_to_ref = reference_basis.get_rotation(membrane_normal)
        # 'subtract' reference rotation from the particle rotation
        rot_par_to_mem = np.dot(rot_mem_to_ref, rot_ref_to_par)

        # convert to zxz, but invert the matrix to accomodate to pytom convention
        zxz_in_mem = matToZXZ(np.linalg.inv(rot_par_to_mem)).toVector(convention='zxz')

        # rotate unit vector by rotation matrix...
        # particle_normal = Vector(unit_vector.get())  # copy reference normal
        # particle_normal.rotate(par_rot_basis)

        # ===> check that above is the same as this >
        # test_normal = Vector(reference_basis.get())
        # test_normal.rotate(par_rot[:3, :3])
        # if sum([abs(a - b) for a, b  in zip(particle_normal.get(), test_normal.get())]) > 0.001:
        #     print(particle_normal.get(), test_normal.get())

        # # rotation matrix of reference vector onto membrane normal
        # # rm_ref_to_membrane = membrane_normal.get_rotation(unit_vector)  # this is the rotation to membrane basis
        # rm_ref_to_membrane = membrane_normal.get_rotation(unit_vector)
        # # remove the base rotation to compare rotation in membrane frame
        # par_rot_memb_basis = np.dot(par_rot_basis, rm_ref_to_membrane.T)
        # # particle_normal.rotate(rm_ref_to_membrane.T) => same as recalculating particle normal as below
        # # print(particle_normal.get())
        #
        # # membrane_normal.rotate(par_rot_basis.T)
        #
        # # rotate unit vector with relative rotation matrix
        # particle_normal = Vector(unit_vector.get())
        # particle_normal.rotate(par_rot_memb_basis)
        # # we are now in the coordinates system of the [0, 0, 1] vector
        # difference_angle = particle_normal.angle(unit_vector, degrees=True)
        #
        # # get rotation around axis
        # rm_unit_to_par = particle_normal.get_rotation(unit_vector)
        # # rot_around_axis = differenceAngleOfTwoRotations(matToZXZ(rm_unit_to_par), matToZXZ(par_rot_memb_basis))
        #
        # remainder = np.dot(rm_unit_to_par, par_rot_memb_basis.T)
        # matrix_test = Matrix(3, 3)
        # matrix_test[0, 0], matrix_test[0, 1], matrix_test[0, 2] = remainder[0]
        # matrix_test[1, 0], matrix_test[1, 1], matrix_test[1, 2] = remainder[1]
        # matrix_test[2, 0], matrix_test[2, 1], matrix_test[2, 2] = remainder[2]
        # angle, axis = matToAxisAngle(matrix_test)
        # sign = axis[2]
        # # print(axis, angle)
        #
        # if sign < 0:
        #     rot_around_axis = angle * -1 + 360
        # else:
        #     rot_around_axis = angle
        #
        # # get the difference angle
        # # difference_angle = particle_normal.angle(membrane_normal, degrees=True)
        # orientations.append(difference_angle)  # probability fit should correct for random angular distribution
        #
        # # get axis rotation
        # # matrix_basic = particle_normal.get_rotation(unit_vector)
        # # rot_around_axis = differenceAngleOfTwoRotations(matToZXZ(matrix_basic), matToZXZ(matrix))
        # orientations_around_axis.append(rot_around_axis)
        # # difference of this angle with the plane???

        # get the distance
        distances.append(distance.euclidean(membrane_point, coordinates))

        # append the variation
        angle_stds.append(angular_variation)

        # create arrows for bild file
        orientations.append(0)
        orientations_around_axis.append(0)
        vector_for_bild = unit_vector.copy()
        # vector_for_bild.rotate(par_rot_basis)
        particle_arrows.append((coordinates, vector_for_bild.get()))
        membrane_arrows.append((membrane_point, membrane_normal.get()))
        final_zxz.append(zxz_in_mem)

    return distances, orientations, orientations_around_axis, angle_stds, particle_arrows, membrane_arrows, final_zxz


# ===================================NUMBA ACCELERATED ANGULAR DISTANCE MATRIX==========================================
@njit
def check_eps(value):
    if abs(value) < epsilon:
        return 0.0
    else:
        return float(value)


@njit
def fill_zmat(z):
    mat = np.identity(3)
    sinz = check_eps(np.sin(z))
    cosz = check_eps(np.cos(z))
    mat[0, 0] = cosz
    mat[0, 1] = -sinz
    mat[1, 0] = sinz
    mat[1, 1] = cosz
    return mat


@njit
def fill_xmat(x):
    mat = np.identity(3)
    sinx = check_eps(np.sin(x))
    cosx = check_eps(np.cos(x))
    mat[1, 1] = cosx
    mat[1, 2] = -sinx
    mat[2, 1] = sinx
    mat[2, 2] = cosx
    return mat


@njit
def diff(z1_1, x_1, z2_1, z1_2, x_2, z2_2):
    """
    Input and output of angular diff in radians.
    """
    z1mat_1 = fill_zmat(z1_1)
    xmat_1 = fill_xmat(x_1)
    z2mat_1 = fill_zmat(z2_1)
    rot_1 = np.dot(np.dot(z1mat_1, xmat_1), z2mat_1)

    z1mat_2 = fill_zmat(z1_2)
    xmat_2 = fill_xmat(x_2)
    z2mat_2 = fill_zmat(z2_2)
    rot_2 = np.linalg.inv(np.dot(np.dot(z1mat_2, xmat_2), z2mat_2))  # invert for diff angle

    rot = np.dot(rot_1, rot_2)

    trace = rot[0, 0] + rot[1, 1] + rot[2, 2]
    cosAng = .5 * (trace - 1)
    if cosAng > 1:
        cosAng = 1.
    if cosAng < -1.:
        cosAng = -1.
    the = np.arccos(cosAng)

    return the


@njit(parallel=True)
def angular_distance(mat):
    newmat = np.zeros((mat.shape[0], mat.shape[1]), dtype=np.float32)  # output matrix shape

    for i in prange(mat.shape[0]):
        for j in range(mat.shape[1]):
            newmat[i, j] = diff(mat[i, j, 0], mat[i, j, 1], mat[i, j, 2],
                                mat[i, j, 3], mat[i, j, 4], mat[i, j, 5])

    for i in prange(mat.shape[0]):
        newmat[i, i] = .0

    return newmat


def angular_distance_matrix(rotations):
    """
    Calculate a angular distance matrix with the 'distance' bewteen all the rotations in the input.
    """
    # construct matrix of (n, n, 6) where is the number zxz of rotations in the input
    ang1 = np.tile(rotations[:, np.newaxis, :], (1, rotations.shape[0], 1))
    ang2 = np.tile(rotations[np.newaxis, :, :], (rotations.shape[0], 1, 1))
    mat = np.deg2rad(np.concatenate((ang1, ang2), axis=2))  # combine to matrix shape

    # pass to calculation of angular distance which returns a (n,n,1) matrix
    return np.rad2deg(angular_distance(mat)) / 180.


# ==============================================MEMBRANE MESH CLASS AND SUPPORT=========================================
def point_3d_in_triangle(point, v1, v2, v3):
    """
    Reference: W. Heidrich, Journal of Graphics, GPU, and Game Tools,Volume 10, Issue 3, 2005
    @param point: coordinate in 3D
    @type point: np.array(3)
    @param v1: vertices defined by xyz
    @type v1: np.array(3)
    @type v2: np.array(3)
    @type v3: np.array(3)
    @return: (projection of point on triangle plane; whether the projection is inside the triangle)
    @rtype: np.array(3), bool
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

    pp = alpha * v1 + beta * v2 + gamma * v3  # the projected point (pp)
    is_in_triangle = (0 <= alpha <= 1) and (0 <= beta <= 1) and (0 <= gamma <= 1)

    return pp, is_in_triangle


def point_3d_in_line(point, a, b):
    """
    @param point: coordinate in 3D
    @type point: np.array(3)
    @param a: first point of line segment
    @type a: np.array(3)
    @param b: second point of line segment
    @type b: np.array(3)
    @return: (projection of point on line (i.e. perpendicular); whether projection is within a, b)
    @rtype: np.array(3), bool
    """
    ap = point - a
    ab = b - a

    # projected point (pp)
    pp = a + np.dot(ap, ab) / np.dot(ab, ab) * ab

    # only if distance of ( p -> a ) + ( p -> b) == (a -> b) the point is in between a and b
    if distance.euclidean(a, b) == \
            distance.euclidean(a, pp) + \
            distance.euclidean(b, pp):
        return pp, True

    return pp, False  # else


def faces_to_edges(faces, vert):
    """
    Convert faces to list of edges connected to vert, where each edge holds the indices to its adjacent faces.
    """
    # dict that maps vertex to the remaining vertices of the face
    index_dict = {0: [1, 2], 1: [0, 2], 2: [0, 1]}

    edges = []

    for face in faces:
        id = face.tolist().index(vert)  # each face can have two possible edges connected to the vert
        rem = index_dict[id]  # an edge is defined by two vertices
        e1 = [face[id], face[rem[0]]]  # => these are the two possible combinations
        e2 = [face[id], face[rem[1]]]

        # append only if not yet recorded
        if not any([True for e in edges if (e == e1 or e == e1[::-1])]):
            edges.append(e1)
        if not any([True for e in edges if (e == e2 or e == e2[::-1])]):
            edges.append(e2)

    return np.array(edges)


def get_edge_vector(edge, faces, verts):
    """
    Find the normal on an edge by averaging the normals of the two adjacent faces.
    """
    # get the indices of the faces next to the edge
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


def find_std(face_normals):
    """
    Calculate the variation of a set of triangle normals.
    """
    # get the mean vector of all the normals: sum all vectors and divice by number of vectors
    mean_vec = Vector(np.array([n.get() for n in face_normals]).sum(axis=0) / len(face_normals))

    # calculate difference of each one to the mean
    diff_angles = np.array(list(map(lambda x: mean_vec.angle(x, degrees=True), face_normals)))

    return diff_angles.std()  # return their sigma


class MembraneMesh:
    def __init__(self, volume, cutoff=0.3, mesh_detail=2, ref_vector=[0, 0, 1]):
        # ensure volume is normalized between 0 and 1
        self.volume = (volume - volume.min()) / (volume.max() - volume.min())
        self.verts, self.faces, self.normals, self.values = \
            measure.marching_cubes(self.volume, level=cutoff, step_size=mesh_detail)
        if ref_vector == [0, 0, 1]:
            self.reference_unit_vector = None
        else:
            self.reference_unit_vector = Vector(ref_vector, normalize=True)
        self.z_axis_unit_vector = Vector([0, 0, 1], normalize=True)

    def display(self):
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Fancy indexing: `verts[faces]` to generate a collection of triangles
        mesh = Poly3DCollection(self.verts[self.faces])
        mesh.set_edgecolor('k')
        ax.add_collection3d(mesh)

        ax.set_xlabel(f"x-axis: {self.volume.shape[0]}")
        ax.set_ylabel(f"y-axis: {self.volume.shape[1]}")
        ax.set_zlabel(f"z-axis: {self.volume.shape[2]}")

        ax.set_xlim(0, self.volume.shape[0])
        ax.set_ylim(0, self.volume.shape[1])
        ax.set_zlim(0, self.volume.shape[2])

        plt.tight_layout()
        plt.show()

    def visualize_vectors(self, particle_vecs, membrane_vecs):
        # Display resulting triangular mesh using Matplotlib. This can also be done
        # with mayavi (see skimage.measure.marching_cubes_lewiner docstring).
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Fancy indexing: `verts[faces]` to generate a collection of triangles
        mesh = Poly3DCollection(self.verts[self.faces])
        mesh.set_edgecolor('k')
        ax.add_collection3d(mesh)

        pdata = np.array([list(c) + list(p) for c, p in particle_vecs])
        X, Y, Z, U, V, W = zip(*pdata)
        ax.quiver(X, Y, Z, [10 * u for u in U], [10 * u for u in V], [10 * u for u in W], color='red')

        mdata = np.array([list(c) + list(p) for c, p in membrane_vecs])
        X, Y, Z, U, V, W = zip(*mdata)
        ax.quiver(X, Y, Z, [10 * u for u in U], [10 * u for u in V], [10 * u for u in W], color='blue')

        ax.set_xlabel(f"x-axis: {self.volume.shape[0]}")
        ax.set_ylabel(f"y-axis: {self.volume.shape[1]}")
        ax.set_zlabel(f"z-axis: {self.volume.shape[2]}")

        ax.set_xlim(0, self.volume.shape[0])
        ax.set_ylim(0, self.volume.shape[1])
        ax.set_zlim(0, self.volume.shape[2])

        plt.tight_layout()
        plt.show()

    def find_membrane_surface(self, coordinate):
        """
        Find the membrane surface closest to a coordinate. Returns closest point on the membrane, the normal of
        the membrane at that position, and the variation of the normals around that position.
        @type coordinate: np.array(3)
        @return: membrane point, membrane normal, standard deviation
        @rtype: np.array(3), pytom.simulation.membrane.Vector, float
        """
        # distance to each vertex in the triangle mesh
        distance_vector = np.sqrt(np.sum(np.subtract(self.verts, coordinate) ** 2, axis=1))

        # select faces containing the closest vertex
        closest_faces = self.faces[np.any(self.faces == np.argmin(distance_vector), axis=1)]

        # point of coordinate of plane formed by triangle, boolean whether this point is inside the triangle,
        # and the normal of the triangle
        face_projections, face_normals, point_on_face = [None, ] * len(closest_faces), \
                                                        [None, ] * len(closest_faces), \
                                                        [None, ] * len(closest_faces)

        for i, face in enumerate(closest_faces):
            # get triangle vertices
            v1, v2, v3 = self.verts[face[0]], self.verts[face[1]], self.verts[face[2]]

            # calculate triangle normal
            normal = Vector(np.cross(v2 - v1, v3 - v1), normalize=True).inverse()  # normal needs to point outwards

            # get particle coordinate projected on triangle face and bool telling whether its
            face_projections[i], point_on_face[i], face_normals[i] = \
                *point_3d_in_triangle(coordinate, v1, v2, v3), normal

        sigma = find_std(face_normals)  # standard deviation of normals

        # consider cases of where the coordinate is relative to the surface
        if sum(point_on_face) == 1:  # this is the easy case, just above a single triangle

            membrane_normal, membrane_point = face_normals[point_on_face.index(True)], \
                                              face_projections[point_on_face.index(True)]

        elif sum(point_on_face) > 1:  # above two or more, select triangle with shortest distance

            # get the id of the closest face
            id, _ = min([(i, distance.euclidean(p, coordinate)) for i, (p, t)
                         in enumerate(zip(face_projections, point_on_face)) if t], key=lambda x: x[1])
            membrane_normal, membrane_point = face_normals[id], face_projections[id]

        else:  # not above any triangle, it is either above an edge or above a vert

            # get all edges connected to closest vert
            edges = faces_to_edges(closest_faces, np.argmin(distance_vector))

            # get edge projections, edge normals, and whether the projected coordinate is on the line segment
            # (see above, similar to faces)
            edge_projections, edge_normals, point_on_edge = [None, ] * len(edges), \
                                                            [None, ] * len(edges), \
                                                            [None, ] * len(edges)

            for i, edge in enumerate(edges):
                edge_projections[i], point_on_edge[i], edge_normals[i] = \
                    *point_3d_in_line(coordinate, self.verts[edge[0]], self.verts[edge[1]]), \
                    get_edge_vector(edge, closest_faces, self.verts)

            if sum(point_on_edge) == 1:  # easy case, above single edge

                membrane_normal, membrane_point = edge_normals[point_on_edge.index(True)], \
                                                  edge_projections[point_on_edge.index(True)]

            elif sum(point_on_edge) > 1:  # above two or more edges

                id, _ = min([(i, distance.euclidean(t, coordinate)) for i, (p, t)
                             in enumerate(zip(edge_projections, point_on_edge)) if t], key=lambda x: x[1])
                membrane_normal, membrane_point = edge_normals[id], edge_projections[id]

            else:  # finally, if not above anything else, select the vertex

                membrane_normal = Vector(self.normals[np.argmin(distance_vector)])
                membrane_point = self.verts[np.argmin(distance_vector)]

        return membrane_point, membrane_normal, sigma

    def find_particle_orientations(self, particle_list, verbose=False):

        # initialize some lists
        distances, angle_stds, final_zxz, particle_arrows, membrane_arrows = \
            [None, ] * len(particle_list), [None, ] * len(particle_list), [None, ] * len(particle_list), \
            [None, ] * len(particle_list), [None, ] * len(particle_list)

        # exclusion_count = 0
        for i, p in enumerate(particle_list):

            # get the coordinates of the particle
            coordinate = np.array(p.getPickPosition().toVector())

            # find closest point on membrane
            membrane_point, membrane_normal, angular_variation = self.find_membrane_surface(coordinate)

            if verbose:
                print('angular variation of membrane normals: ', angular_variation)

            # get rotation matrix and convert to axis-angle
            zxz_par_to_ref = p.getRotation().toVector(convention='zxz')  # z1, z2, x
            rot_ref_to_par = rotation_matrix(rotation=zxz_par_to_ref,
                                             rotation_order='rzxz')[:3, :3].T  # transpose for ref to par

            if self.reference_unit_vector is not None:
                pass  # add additional rotation to particle

            # get rot of membrane surface to z axis unit vec
            rot_mem_to_ref = self.reference_unit_vector.get_rotation(membrane_normal)

            # 'subtract' reference rotation from the particle rotation
            rot_par_to_mem = np.dot(rot_mem_to_ref, rot_ref_to_par)

            # convert to zxz, but invert the matrix to accomodate to pytom convention
            zxz_in_mem = matToZXZ(np.linalg.inv(rot_par_to_mem)).toVector(convention='zxz')

            # create arrows for bild file
            particle_arrows[i] = (coordinate, self.reference_unit_vector.rotate(rot_ref_to_par).get())
            membrane_arrows[i] = (membrane_point, membrane_normal.get())
            distances[i], angle_stds[i], final_zxz[i] = \
                distance.euclidean(membrane_point, coordinate), angular_variation, zxz_in_mem

        return distances, angle_stds, particle_arrows, membrane_arrows, final_zxz


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


def write_clusters(particle_list, clusters):
    for i in np.unique(clusters):
        if i == -1:  # write a list for each cluster that is not 0
            continue
        else:
            cluster_list = ParticleList()
            for c, p in zip(clusters, particle_list):
                if i == c:
                    cluster_list.append(p)
            cluster_list.toXMLFile(f'particles_c{i}.xml')
            print(f'wrote particles_c{i}.xml with {len(cluster_list)} particles')


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
            ScriptOption2(['-c', '--cutoff'], 'Cutoff value for converting membrane model to a triangular mesh. '
                                              'IMPORTANT: make sure to check segmentation for a threshold that does '
                                              'not include artifacts.',
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
            ScriptOption2(['--interactive-clustering'], 'Whether to do clustering interactively to find optimal '
                                                        'parameters.',
                          'no arguments', 'optional'),
            ScriptOption2(['--verbose'], 'Be verbose.', 'no arguments', 'optional')])

    options = parse_script_options2(sys.argv[1:], helper)

    input_file, segmentation_file, output_name, voxel_size, cutoff, mesh_detail, template_normal, nstd, \
        distance_range, angle_std_cutoff, axis_distance_limit, nbins, split_lists_angle, interactive_clustering,\
        verbose = options

    interactive_clustering = False if interactive_clustering is None else True

    if output_name is not None and not os.path.exists(output_name):
        os.mkdir(output_name)

    _, ext = os.path.splitext(input_file)

    if ext == '.xml':

        # load particles
        particle_list = ParticleList()
        particle_list.fromXMLFile(input_file)

        # create membrane mesh
        membrane = MembraneMesh(read(segmentation_file), cutoff=cutoff, mesh_detail=mesh_detail,
                                ref_vector=template_normal)

        dists, stds, p_arrows, m_arrows, zxzs = membrane.find_particle_orientations(particle_list, verbose=verbose)

    elif ext == '.txt':

        linker_data = np.genfromtxt(input_file, dtype=[('segmentation', 'U1000'), ('particle_list', 'U1000')],
                                    skip_header=1)
        segmentation_files = linker_data['segmentation']
        particle_list_files = linker_data['particle_list']

        distances, stds, zxzs = [], [], []
        combined_particle_list = ParticleList()

        for s, p in zip(segmentation_files, particle_list_files):
            print(f'Run on segmentation {s} with plist {p}')
            name = os.path.splitext(os.path.split(s)[1])[0]

            # read particles
            particle_list = ParticleList()
            particle_list.fromXMLFile(p)

            # create membrane mesh
            membrane = MembraneMesh(read(s), cutoff=cutoff, mesh_detail=mesh_detail, ref_vector=template_normal)

            d, std, p_arrows, m_arrows, zxz = membrane.find_particle_orientations(particle_list, verbose=verbose)
            distances += d
            stds += std
            zxzs += zxz

            # add to the full length particle list
            combined_particle_list += particle_list

            # membrane_arrow_test += m_arrows
            write_arrow_bild(p_arrows, m_arrows, os.path.join(output_name, name + '.bild'))

    else:
        print('Invalid input extension.')
        sys.exit(0)

    # total number of particles
    # n_pre_filter = len(orientations)

    # convert to numpy arrays for calculation
    distances = np.array(distances)
    distances *= voxel_size if voxel_size is not None else 1
    # orientations = np.array(orientations)
    # orientations_around_axis = np.array(orientations_around_axis)
    stds = np.array(stds)
    zxzs = np.array(zxzs)
    z1_angles = zxzs[:, 0]
    x_angles = zxzs[:, 1]
    z2_angles = zxzs[:, 2]

    n_pre_filter = zxzs.shape[0]

    # find distance cutoff  and angle cutoff if not specified
    distance_range = (distances.mean() - nstd * distances.std(), distances.mean() + nstd * distances.std()) if \
        distance_range is None else distance_range
    angle_std_cutoff = nstd * stds.std() + stds.mean() if angle_std_cutoff is None else angle_std_cutoff

    # if split_lists_angle is not None and output_name is not None:
    #     write_split_particle_lists(orientations_per_particle_list, distance_range, angle_std_cutoff,
    #                                split_lists_angle, voxel_size, output_name)

    # Filter the combined dataset to get better statistics of distributions for plot
    distance_outlier_filter = np.logical_and(distances >= distance_range[0], distances <= distance_range[1])
    angle_outlier_filter = (stds < angle_std_cutoff)
    outlier_filter = np.logical_and(distance_outlier_filter, angle_outlier_filter)
    distances_filtered = distances[outlier_filter]
    # orientations = orientations[outlier_filter]
    # orientations_around_axis = orientations_around_axis[outlier_filter]
    z1_angles, x_angles, z2_angles = z1_angles[outlier_filter], x_angles[outlier_filter], z2_angles[outlier_filter]
    n = len(distances_filtered)

    # output number of particle excluded due to filtering
    total_outliers = (outlier_filter == False).sum()
    print(f'Filtered {total_outliers} from total of {zxzs.shape[0]} particles')

    plt.hist(z1_angles)  # => z1 angles are completely random, as they should be:
    plt.show()

    rot_filtered = zxzs[outlier_filter]
    rot_filtered[:, 0] = 0.
    ang_matrix = angular_distance_matrix(rot_filtered)

    # use dbscan, eps is the distance of points to be considered in the neighbourhood
    # model = cluster.OPTICS(metric='precomputed', max_eps=0.04, min_cluster_size=0.1, n_jobs=8).fit(ang_matrix)
    model = hdbscan.HDBSCAN(metric='precomputed',
                            min_cluster_size=int(0.067 * ang_matrix.shape[0]),
                            min_samples=int(0.013 * ang_matrix.shape[0]),
                            cluster_selection_epsilon=0.)
    model.fit(ang_matrix.astype(np.float64))

    while interactive_clustering:
        min_cluster = int(input('Minimum cluster size: '))
        min_samples = int(input('Minimum samples size: '))
        eps = float(input('Epsilon: '))
        # model = cluster.OPTICS(metric='precomputed', max_eps=eps, min_cluster_size=min_cluster, min_samples=min_samples,
        #                        n_jobs=8, algorithm='brute')
        model = hdbscan.HDBSCAN(metric='precomputed',
                                min_cluster_size=min_cluster,  # int(0.05 * ang_matrix.shape[0]),
                                min_samples=min_samples,  # int(0.02 * ang_matrix.shape[0]))
                                cluster_selection_epsilon=eps)
        model.fit(ang_matrix.astype(np.float64))
        print(np.unique(model.labels_))

        fig, ax = plt.subplots(figsize=(4, 8))
        ax.scatter(x_angles, z2_angles, alpha=0.7, s=4, c=model.labels_.astype(float))
        ax.set_xlim(0, 180)
        ax.set_xlabel('euler X')
        ax.set_ylim(0, 360)
        ax.set_ylabel('euler Z2')
        plt.show()

        interactive_clustering = int(input('Continue with interactive clustering (0 or 1): '))
        interactive_clustering = True if interactive_clustering == 1 else False

    # write clusters as particle lists
    new_list = ParticleList()
    for i, (f, p) in enumerate(zip(outlier_filter, combined_particle_list)):
        new_list.append(p)
    write_clusters(new_list, model.labels_)

    # ============================== create plots
    # try:
    #     import matplotlib.pyplot as plt
    # except Exception as e:
    #     import matplotlib
    #     matplotlib.use('Qt5Agg')
    #     import matplotlib.pyplot as plt
    #
    # # fig = plt.figure(figsize=(5, 5))
    # # ax = fig.add_subplot(111, projection='3d')
    # #
    # # points = np.array([p for o, p in membrane_arrow_test]).T
    # # ax.scatter(points[0], points[1], points[2], s=1)
    # # # ax.quiver(0, 0, 0, *norm, length=100, color='red')
    # #
    # # ax.set_xlabel('x')
    # # # ax.set_xlim(-300, 300)
    # # ax.set_ylabel('y')
    # # # ax.set_ylim(-300, 300)
    # # ax.set_zlabel('z')
    # # # ax.set_zlim(-300, 300)
    # # plt.show()
    #
    # font = {'size': 16}
    # plt.rc('font', **font)
    #
    #
    #
    # # do some plotting
    # fig, ax = plt.subplots(figsize=(5, 5))
    #
    # ax.hist(orientations, bins=nbins, color='black', histtype='stepfilled', alpha=0.5)
    # h = ax.hist(orientations, bins=nbins, color='black', histtype='step')
    # values, bins = h[0], h[1]
    # ax.set_ylabel('Number of particles')
    # ax.set_xlabel('Angle (degrees)')
    #
    # # information about number of particles
    # textstr = f'N={n}'
    # # these are matplotlib.patch.Patch properties
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
    # # place a text box in upper left in axes coords
    # ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    #            verticalalignment='top', bbox=props)
    #
    # plt.tight_layout()
    #
    # if output_name is not None:
    #     plt.savefig(os.path.join(output_name, 'orientation_distribution.png'), dpi=300, format='png', transparent=True)
    #     plt.close()
    #     print(f"wrote {os.path.join(output_name, 'orientation_distribution.png')}")
    #     if ext == '.xml':
    #         write_arrow_bild(p_arrows, m_arrows, os.path.join(output_name, 'vectors.bild'),
    #                          outlier_filter=outlier_filter)
    #         print(f"wrote {os.path.join(output_name, 'vectors.bild')}")
    # else:
    #     plt.show()
    #
    # from scipy.stats import kstest # compare distribution to scipy.stats.uniform ??
    # # => this is the kolmogorov smirnov test
    #
    # # statistical test whether angles come from sin distribution
    # random_z = np.random.uniform(-1, 1, size=len(orientations))
    # random_angles = np.rad2deg(np.arccos(random_z))
    # ks_result = kstest(orientations, random_angles)
    # print(ks_result)
    #
    # # ks plot
    # sort_rand = sorted(random_angles)
    # sort_orie = sorted(orientations)
    # p = 1. * np.arange(len(sort_rand)) / (len(sort_rand) - 1)
    # fig, ax = plt.subplots()
    # ax.plot(p, sort_rand, label='random')
    # ax.plot(p, sort_orie, label='data')
    # ax.legend()
    # ax.set_title(f'D({len(orientations)})={ks_result.statistic:.2f}, p={ks_result.pvalue}')
    #
    # plt.tight_layout()
    #
    # if output_name is not None:
    #     plt.savefig(os.path.join(output_name, 'kolmogorov-smirnov.png'), dpi=300, format='png',
    #                 transparent=True)
    #     plt.close()
    #     print(f"wrote {os.path.join(output_name, 'kolmogorov-smirnov.png')}")
    # else:
    #     plt.show()
    #
    # # do some plotting
    # fig, ax = plt.subplots(figsize=(5, 5))
    #
    # angles = (bins[:-1] + bins[1:]) / 2
    # # print(angles)
    # random_prob = np.sin(np.deg2rad(angles)) # * -1 + 1
    #
    # # ax.plot(angles, random_prob / random_prob)
    # # ax.plot(angles, values / random_prob)
    # ax.hist(angles, weights=values / random_prob, color='black', alpha=0.5, histtype='stepfilled', bins=nbins)
    # ax.hist(angles, weights=values / random_prob, color='black', histtype='step', bins=nbins)
    # # ax.hlines(0.5, xmin=-1, xmax=1, linestyles='dashed', color='black', label='random dist')
    # # ax.legend()
    # ax.set_ylabel('Number of particles')
    # ax.set_xlabel('Angle (degrees)')
    # ax.set_xlim(0, 180)
    #
    # # information about number of particles
    # textstr = f'N={n}'
    # # these are matplotlib.patch.Patch properties
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
    # # place a text box in upper left in axes coords
    # ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    #         verticalalignment='top', bbox=props)
    #
    # plt.tight_layout()
    #
    # if output_name is not None:
    #     plt.savefig(os.path.join(output_name, 'orientation_distribution_normalized.png'), dpi=300, format='png',
    #                 transparent=True)
    #     plt.close()
    #     print(f"wrote {os.path.join(output_name, 'orientation_distribution_normalized.png')}")
    # else:
    #     plt.show()

    # ==================== create plot of distance vs. angles

    # fig, ax = plt.subplots(figsize=(5, 5))
    # ax.scatter(distances_filtered, orientations, alpha=0.7, color='gray')
    # if voxel_size is not None:
    #     ax.set_xlabel(r'Distance ($\AA$)')
    # else:
    #     ax.set_xlabel('Distance (voxels)')
    # ax.set_ylabel('Angle (degrees)')
    #
    # # information about number of particles
    # textstr = f'N={n}'
    # # these are matplotlib.patch.Patch properties
    # props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
    # # place a text box in upper left in axes coords
    # ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    #            verticalalignment='top', bbox=props)
    #
    # plt.tight_layout()
    #
    # if output_name is not None:
    #     plt.savefig(os.path.join(output_name, 'scatter.png'), dpi=300, format='png', transparent=True)
    #     plt.close()
    #     print(f"wrote {os.path.join(output_name, 'scatter.png')}")
    # else:
    #     plt.show()

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
        plt.close()
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
        plt.close()
        print(f"wrote {os.path.join(output_name, 'std.png')}")
    else:
        plt.show()

    # ============ test angle vs angle
    # => scipy.stats.gaussian_kde
    # could do a 2d kde on the data and then correct the density for random angular distribution of the vector angle

    # plt.hist(x_angles)
    # plt.show()
    # plt.hist(z2_angles)
    # plt.show()

    fig, ax = plt.subplots(figsize=(4.5, 7))
    h = ax.hist2d(x_angles, z2_angles, bins=(nbins//2, nbins))
    plt.colorbar(h[3], ax=ax)
    ax.set_xlabel('euler X')
    ax.set_ylabel('euler Z2')
    # plt.tight_layout()
    if output_name is not None:
        plt.savefig(os.path.join(output_name, '2d_rot_density.png'), dpi=300, format='png',
                    bbox_inches='tight')
        plt.close()
        print(f"wrote {os.path.join(output_name, '2d_rot_density.png')}")
    else:
        plt.show()

    # print(x_angles.min(), x_angles.max())
    # print(z2_angles.min(), z2_angles.max())

    fig, ax = plt.subplots(figsize=(4, 8))
    ax.scatter(x_angles, z2_angles, alpha=0.7, s=4, c=model.labels_.astype(float))
    ax.set_xlim(0, 180)
    ax.set_xlabel('euler X')
    ax.set_ylim(0, 360)
    ax.set_ylabel('euler Z2')
    # plt.tight_layout()
    if output_name is not None:
        plt.savefig(os.path.join(output_name, '2d_rot_scatter.png'), dpi=300, format='png',
                    bbox_inches='tight')
        plt.close()
        print(f"wrote {os.path.join(output_name, '2d_rot_scatter.png')}")
    else:
        plt.show()
