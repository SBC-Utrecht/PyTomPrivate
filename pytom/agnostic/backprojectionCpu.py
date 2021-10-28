from numba import jit, njit, prange
import numpy as np
import numpy as cp

@njit(parallel=True)
def singlewbp(reconstruction, src, dims, dims_src, center_recon, tr11, tr21, tr31, tr41, tr12, tr22, tr32, tr42):

    # tr11, tr21, tr31, tr41 = mtx[:,0]
    # tr12, tr22, tr32, tr42 = mtx[:,1]


    for x in prange(dims[0]):
        for y in prange(dims[1]):
            for z in prange(dims[2]):
                # x = i  // (dims[1] * dims[2])
                # y = (i // dims[2]) % dims[1]
                # z = (i % dims[2])

                #if np.sqrt((x-center_recon[0])**2+(z-center_recon[2])**2) > max(dims[0], dims[2])//2: continue

                # Unexpected +1 needed for interpolation_positions -- is this a BUG inhereted from c implentation!
                interpolation_position_x = tr11 * (x - center_recon[0] ) + tr21 * (y - center_recon[1] ) + tr31 * (z - center_recon[2] ) + tr41 #+ center_proj[0] + 1
                interpolation_position_y = tr12 * (x - center_recon[0] ) + tr22 * (y - center_recon[1] ) + tr32 * (z - center_recon[2] ) + tr42 #+ center_proj[0] + 1

                # interpolation_position_x = tr11 * x + tr31 * z + center_proj[n, 0]
                # interpolation_position_y = tr22 * y + center_proj[n, 1]

                px = int(interpolation_position_x)
                py = int(interpolation_position_y)

                ix = float(interpolation_position_x - float(px))
                iy = float(interpolation_position_y - float(py))

                if (px >= 1) and (px <= dims_src[0]) and (py >= 1) and (py <= dims_src[1]):
                    mipx = px - 1
                    mipy = py - 1
                    value_x1 = src[mipx,mipy] if px == dims_src[0] else src[mipx,mipy]*(1-ix) + src[px,mipy]*ix
                    if py < dims_src[1]:
                        value_x2 = src[mipx, py] if px==dims_src[0] else src[mipx, py]*(1-ix) + src[px,py]*ix
                        #value_x1 += (src[proj_position_x, proj_position_y - 1] - src[proj_position_x - 1, proj_position_y - 1]) * interpolation_offset_x

                    # if proj_position_y == dims_src[1]:
                    #     value_x2 = 0.
                    #     interpolation_offset_y = 0.
                    # elif proj_position_x==dims_src[0]:
                    #     value_x2 = src[proj_position_x - 1, proj_position_y]
                    # else:
                    #     value_x2 = src[proj_position_x - 1, proj_position_y] + \
                    #                (src[proj_position_x, proj_position_y] -src[proj_position_x - 1, proj_position_y]) * interpolation_offset_x

                    value_y = value_x1 if py == dims_src[1] else value_x1 + (value_x2 - value_x1) * iy
                    reconstruction[x,y,z] = reconstruction[x,y,z] + value_y




def backProjectCPU(projections, vol_bp, phi_angles, proj_angles, recPosVol=None, vol_offsetProjections=None,
                   interpolation='', particleObject=None, alignment_params=None, alignment_order=(2,1,0)):
    from pytom.agnostic.io import read_size

    if recPosVol.squeeze().shape[0] == 3:
        recPosVol = recPosVol.squeeze().T
    else:
        recPosVol = recPosVol.squeeze()

    if vol_offsetProjections.squeeze().shape[0] == 2:
        vol_offsetProjections = vol_offsetProjections.squeeze().T
    else:
        vol_offsetProjections = vol_offsetProjections.squeeze()


    proj_angles = proj_angles.squeeze()

    # preparation
    reconstruction = cp.array(vol_bp, dtype=cp.float32)
    psi_angles   = cp.zeros((proj_angles.shape[0]), dtype=cp.float32)
    phi_angles   = cp.deg2rad(cp.array(phi_angles.squeeze()), dtype=cp.float32)
    theta_angles = cp.deg2rad(cp.array(proj_angles), dtype=cp.float32)

    dims = cp.array(reconstruction.shape, dtype=cp.int32)
    dims_src = cp.array(projections.shape, dtype=cp.int32)

    assert len(theta_angles) == projections.shape[2]  # 'Number of angles and projections should match'

    center_recon = cp.zeros((3), dtype=cp.int32)
    center_recon[0] = int(dims[0] // 2 - recPosVol[0, 0] )
    center_recon[1] = int(dims[1] // 2 - recPosVol[0, 1] )
    center_recon[2] = int(dims[2] // 2 - recPosVol[0, 2] )

    dims_proj = cp.array(projections.shape, dtype=cp.int32)

    center_proj = cp.zeros_like(vol_offsetProjections)
    center_proj[:, 0] += dims_proj[0] // 2 + vol_offsetProjections[:, 0]
    center_proj[:, 1] += dims_proj[1] // 2 + vol_offsetProjections[:, 1]


    mtx = cp.identity(4, dtype=cp.float32)
    from numpy import cos, sin
    from pytom.voltools.utils.matrices import rotation_matrix, translation_matrix, scale_matrix

    subtomoAlignment = cp.identity(4, dtype=cp.float32)
    imageAlignment = cp.identity(4, dtype=cp.float32)


    if not particleObject is None:
        t = translation_matrix(translation=particleObject.getShift().toVector())
        r = rotation_matrix(rotation=particleObject.getRotation().toVector())
        subtomoAlignment = r.dot(t)

    if not alignment_params is None:
        t = translation_matrix(translation=(-alignment_params['AlignmentTransX'], -alignment_params['AlignmentTransY'], 0))
        r = rotation_matrix(rotation=(-alignment_params['InPlaneRotation'],0,0), rotation_order='rzxz')
        s = scale_matrix(coefficients=(alignment_params['Magnification'], alignment_params['Magnification'], 1))

        sx,sy,sz = read_size(alignment_params['FileName'])

        postC = translation_matrix(translation=(-sx//2, -sy//2, -sz//2))
        preC = translation_matrix(translation=(sx//2, sy//2, sz//2))

        s = postC.dot(s.dot(preC))
        r = postC.dot(r.dot(preC))

        a = [r,t,s]

        for i in range(alignment_order):
            imageAlignment = a[i].dot(imageAlignment)


    for n in range(projections.shape[2]):
        p = translation_matrix(translation=(-center_proj[n,0]-1, -center_proj[n,1]-1,0))

        # psi = psi_angles[n]
        # phi = phi_angles[n]
        # theta = theta_angles[n]
        Z1 = Z2 = 0.0
        Y = theta_angles[n]

        # tr11 = cos(Y) * cos(Z1) * cos(Z2) - sin(Z1) * sin(Z2)
        # tr21 = cos(Y) * sin(Z1) * cos(Z2) + cos(Z1) * sin(Z2)
        # tr31 = -sin(Y) * cos(Z2)
        # tr41 = 0.
        # tr12 = -cos(Y) * cos(Z1) * sin(Z2) - sin(Z1) * cos(Z2)
        # tr22 = -cos(Y) * sin(Z1) * sin(Z2) + cos(Z1) * cos(Z2)
        # tr32 = sin(Y) * sin(Z2)
        # tr42 = 0.

        tiltAxisMtx = rotation_matrix(rotation=(Z1,Y,Z2), rotation_order='rzyz', rotation_units='rad').astype(cp.float32)
        mtx = subtomoAlignment.dot(p.dot(tiltAxisMtx.dot(imageAlignment)))

        # print(mtx)
        # print(tr11, tr21, tr31, tr12, tr22, tr32)
        # # print(mtx)
        # #raise Exception('pause')
        #
        # # wbp_single_image(reconstruction, projections[:,:,n], tr11, tr21, tr31, tr12, tr22, tr32,
        # #                  dims=dims, center_recon=center_recon, center_proj=center_proj[n,:])
        #
        tr11, tr21, tr31, tr41 = mtx[0, :]
        tr12, tr22, tr32, tr42 = mtx[1, :]

        singlewbp(reconstruction, cp.array(projections[:,:,n],dtype=cp.float32), dims, dims_src, center_recon,
                  tr11, tr21, tr31, tr41, tr12, tr22, tr32, tr42)

    rec = reconstruction

    return rec

