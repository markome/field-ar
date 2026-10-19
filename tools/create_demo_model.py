"""Create a physically motivated demo of field lines around an ideal circular current loop.

Requires: numpy, scipy, trimesh
The sample is for demonstrating AR registration, not measurement-grade simulation.
Coordinates are in 'marker-width' units: the physical HIRO marker black square is
represented as 1 unit wide. Its plane is XZ, and Y points out of the paper.
"""
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipk, ellipe
import trimesh

RADIUS = 0.55
COIL_HEIGHT = 0.60
OUT = Path(__file__).resolve().parents[1] / 'models' / 'demo_fieldlines.glb'


def field_rz(v):
    """Axisymmetric Biot-Savart field components (relative scale only)."""
    r, z = float(v[0]), float(v[1])
    d = (RADIUS - r)**2 + z*z
    s = (RADIUS + r)**2 + z*z
    m = np.clip(4 * RADIUS * r / s, 0., 1. - 1e-10)
    K, E = ellipk(m), ellipe(m)
    if r < 1e-8:
        return np.array([0., 2 * np.pi * RADIUS**2 / (RADIUS**2 + z*z)**1.5])
    br = z / (r*np.sqrt(s)) * (-K + (RADIUS**2+r*r+z*z)/d * E)
    bz = (K + (RADIUS**2-r*r-z*z)/d * E) / np.sqrt(s)
    return np.array([br, bz])


def streamline(seed_r):
    """Follow unit-B until the streamline closes around the conductor."""
    start = np.array([seed_r, 0.])
    def derivative(_, p):
        B = field_rz(p)
        return B / np.linalg.norm(B)
    sol = solve_ivp(derivative, (0., 15.), start, max_step=.008,
                    rtol=1e-7, atol=1e-8)
    d = np.linalg.norm(sol.y.T - start, axis=1)
    selected = np.where((sol.t > .45) & (d < .007))[0]
    if not len(selected):
        raise RuntimeError(f'Failed to close field line for seed {seed_r}')
    end_index = selected[0]
    path = sol.y.T[:end_index+1]
    return path


def tube(points, radius, color, sides=8):
    """Construct a closed low-poly tube mesh, avoiding WebGL 1-pixel line primitives."""
    P = np.asarray(points, dtype=np.float64)
    T = np.roll(P, -1, axis=0) - np.roll(P, 1, axis=0)
    T = T / np.linalg.norm(T, axis=1)[:, None]
    # Start with Y as reference unless path tangent is almost parallel to Y.
    ref = np.broadcast_to([0., 1., 0.], P.shape).copy()
    parallel = np.abs(np.einsum('ij,ij->i', T, ref)) > .93
    ref[parallel] = [1., 0., 0.]
    N = np.cross(T, ref)
    N = N / np.linalg.norm(N, axis=1)[:, None]
    B = np.cross(T, N)
    theta = np.linspace(0, 2*np.pi, sides, endpoint=False)
    v = (P[:, None, :] + radius * (
        np.cos(theta)[None,:,None]*N[:,None,:] +
        np.sin(theta)[None,:,None]*B[:,None,:]
    )).reshape(-1,3)
    f = []
    for i in range(len(P)):
        j = (i+1) % len(P)
        for k in range(sides):
            l = (k+1)%sides
            a, b, c, d = i*sides+k, i*sides+l, j*sides+k, j*sides+l
            f.append((a,c,b)); f.append((b,c,d))
    mesh = trimesh.Trimesh(vertices=v, faces=np.asarray(f), process=False)
    mesh.visual.vertex_colors = np.tile(np.array(color, dtype=np.uint8), (len(v),1))
    return mesh


def make_model():
    scene=trimesh.Scene()
    # The coil itself is a calibration reference in the sample; optional for real overlay.
    phi=np.linspace(0,2*np.pi,192,endpoint=False)
    coil=np.column_stack([RADIUS*np.cos(phi), np.full_like(phi,COIL_HEIGHT), RADIUS*np.sin(phi)])
    scene.add_geometry(tube(coil,.018,[245,156,63,255]), geom_name='reference_conductor')

    seeds=[.36,.45,.75,1.05]
    colors=[[39,216,234,255],[63,175,245,255],[109,135,244,255],[181,142,255,255]]
    for i,r in enumerate(seeds):
        trace=streamline(r)
        # Uniform arc-length resampling and periodic closure adjustment.
        dist=np.r_[0., np.cumsum(np.linalg.norm(np.diff(trace,axis=0),axis=1))]
        total=dist[-1]
        s=np.linspace(0,total,180,endpoint=False)
        rz=np.column_stack([np.interp(s,dist,trace[:,0]), np.interp(s,dist,trace[:,1])])
        # Equal-angle meridional planes; the circular loop has rotational symmetry.
        for j,alpha in enumerate([0., np.pi/4, np.pi/2, 3*np.pi/4]):
            P=np.column_stack([rz[:,0]*np.cos(alpha), rz[:,1]+COIL_HEIGHT, rz[:,0]*np.sin(alpha)])
            scene.add_geometry(tube(P,.009,colors[i]), geom_name=f'field_seed{i}_azimuth{j}')
    OUT.write_bytes(scene.export(file_type='glb'))
    print(f'Created {OUT}: {OUT.stat().st_size:,} bytes, {len(scene.geometry)} meshes')

if __name__=='__main__':
    make_model()
