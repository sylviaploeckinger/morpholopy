import numpy as np
from scipy.spatial.transform import Rotation as R

def apply_rotation(pos,momentum):
    
    # Apply rotation
    Lz = momentum[2]
    L = np.sqrt(momentum[0]**2+momentum[1]**2+momentum[2]**2)
    Lyz = np.sqrt(momentum[1]**2+momentum[2]**2)
    phi = np.arccos(Lz/Lyz)
    if momentum[1]<0: phi = 2*np.pi-phi
    
    pstars = pos[:,0:3].copy()
    rotation = R.from_euler('x', [phi], degrees=False)
    pstars_rot = rotation.apply(pstars)
    momentum_rot = rotation.apply(momentum)
    
    Lz = momentum_rot[0][2]
    Lxz = np.sqrt(momentum_rot[0][0]**2+momentum_rot[0][2]**2)
    beta = np.arccos(Lz/Lxz)
    if momentum[1]<0: beta = 2*np.pi-beta
    rotation = R.from_euler('y', [beta], degrees=False)
    pstars_rot = rotation.apply(pstars_rot)
    momentum_rot = rotation.apply(momentum_rot)
    
    Lx = momentum_rot[0][0]
    Lxy = np.sqrt(momentum_rot[0][0]**2+momentum_rot[0][1]**2)
    gamma = np.arccos(Lx/Lxy)
    rotation = R.from_euler('z', [gamma], degrees=False)
    pstars_rot = rotation.apply(pstars_rot)
    momentum_rot = rotation.apply(momentum_rot)
    
    return pstars_rot
