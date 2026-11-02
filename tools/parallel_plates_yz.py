"""3D BEM model dveh tankih kovinskih plošč v ravnini YZ.
Zahteva: numpy, matplotlib, scikit-image.
"""
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from numpy.polynomial.legendre import leggauss
from skimage.measure import marching_cubes

EPS0 = 8.8541878128e-12
PLATE_Y = 0.15              # m
PLATE_Z = 0.06              # m
GAP = 0.02                  # m (spremenljiv)
VOLTAGE = 1000.0            # V
N_EQUIPOTENTIALS = 7        # spremenljivo, >= 0
NY, NZ = 30, 12             # 5-mm ploščice
QUAD_ORDER = 3
SEED_Y, SEED_Z = 9, 5
STEP = 0.0005               # 0.5 mm
GRID_SHAPE = (38, 42, 64)   # x, y, z; manjše = hitreje
OUT_JSON = 'field_lines_yz.json'

K = 1 / (4 * np.pi * EPS0)
dy, dz = PLATE_Y / NY, PLATE_Z / NZ
ys = np.linspace(-PLATE_Y/2 + dy/2, PLATE_Y/2 - dy/2, NY)
zs = np.linspace(-PLATE_Z/2 + dz/2, PLATE_Z/2 - dz/2, NZ)
Y, Z = np.meshgrid(ys, zs, indexing='ij')
yz = np.column_stack((Y.ravel(), Z.ravel()))
nplate = len(yz)
positive = np.column_stack((np.full(nplate, -GAP/2), yz))
negative = np.column_stack((np.full(nplate, +GAP/2), yz))
centres = np.vstack((positive, negative))
N = len(centres)
Vtarget = np.r_[np.full(nplate, +VOLTAGE/2), np.full(nplate, -VOLTAGE/2)]

# Integracija konstantne povrsinske gostote naboja po vsakem panelu
nodes, w1 = leggauss(QUAD_ORDER)
offsets = np.array([(0, a*dy/2, b*dz/2) for a in nodes for b in nodes])
weights = np.array([wa*wb*dy*dz/4 for wa in w1 for wb in w1])

P = np.zeros((N, N))
for offset, weight in zip(offsets, weights):
    source = centres + offset
    distance = np.linalg.norm(centres[:, None, :] - source[None, :, :], axis=2)
    # diagonalo zamenjamo z analiticnim self-termom
    np.fill_diagonal(distance, np.inf)
    P += weight / distance

a, b = dy/2, dz/2
self_integral = 4*(a*np.arcsinh(b/a) + b*np.arcsinh(a/b))
np.fill_diagonal(P, self_integral)
sigma = np.linalg.solve(K*P, Vtarget)
charge = sigma[:nplate].sum()*dy*dz
print(f'Kapacitivnost: {charge / VOLTAGE * 1e12:.3f} pF')

source_points = (centres[:, None, :] + offsets[None, :, :]).reshape(-1, 3)
source_charges = (sigma[:, None] * weights[None, :]).reshape(-1)

# Enaka BEM resitev za polje in potencial
# Chunking varuje porabo RAM, predvsem pri izracunu 3D mreze.
def evaluate(points, want_field=False, chunk=120):
    pts = np.atleast_2d(np.asarray(points, float))
    result = np.empty((len(pts), 3) if want_field else len(pts))
    for i in range(0, len(pts), chunk):
        r = pts[i:i+chunk, None, :] - source_points[None, :, :]
        r2 = np.einsum('ijk,ijk->ij', r, r)
        r2 = np.maximum(r2, 1e-20)
        if want_field:
            factors = K*source_charges[None, :] / (r2*np.sqrt(r2))
            result[i:i+chunk] = np.einsum('ij,ijk->ik', factors, r)
        else:
            result[i:i+chunk] = K*np.sum(source_charges[None, :]/np.sqrt(r2), axis=1)
    return result

def E(point):
    return evaluate(np.asarray(point), want_field=True)[0]

def direction(p):
    v = E(p)
    return v / max(np.linalg.norm(v), 1e-30)

def rk4(p, h):
    k1 = direction(p)
    k2 = direction(p + h*k1/2)
    k3 = direction(p + h*k2/2)
    k4 = direction(p + h*k3)
    return p + h*(k1 + 2*k2 + 2*k3 + k4)/6

# Silnice potujejo s pozitivne proti negativni plosci (v smeri +X).
# Za zarise od 0.7 mm do povrsine je izbran kratek zacetni odmik.
start_eps = 0.0007
lines = []
for y0 in np.linspace(-.92*PLATE_Y/2, .92*PLATE_Y/2, SEED_Y):
    for z0 in np.linspace(-.94*PLATE_Z/2, .94*PLATE_Z/2, SEED_Z):
        p = np.array([-GAP/2 + start_eps, y0, z0], float)
        path = [p.copy()]
        for _ in range(250):
            pnew = rk4(p, STEP)
            # Preprecimo nadaljevanje v prevodnik; silnice, ki
            # zaidejo izven roba plosce, se zaradi prikaza omejijo.
            if pnew[0] >= GAP/2 - start_eps:
                if pnew[0] != p[0]:
                    t = (GAP/2 - start_eps - p[0])/(pnew[0]-p[0])
                    pnew = p + np.clip(t,0,1)*(pnew-p)
                path.append(pnew.copy())
                break
            if abs(pnew[1]) > 2*PLATE_Y or abs(pnew[2]) > 2*PLATE_Z:
                break
            path.append(pnew.copy())
            p = pnew
        lines.append(np.asarray(path))

# Omejitev ekvipotencialnih ploskev na prikazano polje:
# enak prostorski obseg YZ kot silnice; X znotraj obeh plosc.
all_line_points = np.vstack(lines)
pad_y, pad_z = 0.002, 0.002
ylim = (min(all_line_points[:,1].min()-pad_y, -PLATE_Y/2),
        max(all_line_points[:,1].max()+pad_y, +PLATE_Y/2))
zlim = (min(all_line_points[:,2].min()-pad_z, -PLATE_Z/2),
        max(all_line_points[:,2].max()+pad_z, +PLATE_Z/2))
xgrid = np.linspace(-GAP/2 + start_eps, GAP/2 - start_eps, GRID_SHAPE[0])
ygrid = np.linspace(*ylim, GRID_SHAPE[1])
zgrid = np.linspace(*zlim, GRID_SHAPE[2])

# Povrsine, ki predstavljajo fiksne potenciale, iz 3D skalarnega polja.
if N_EQUIPOTENTIALS > 0:
    gx, gy, gz = np.meshgrid(xgrid, ygrid, zgrid, indexing='ij')
    gridpts = np.column_stack((gx.ravel(), gy.ravel(), gz.ravel()))
    potential = evaluate(gridpts).reshape(GRID_SHAPE)
    # Izpustimo sami plosci (+/- V/2); ravni segajo od -0.8 do +0.8.
    levels = np.linspace(-0.8*VOLTAGE/2, 0.8*VOLTAGE/2,
                         N_EQUIPOTENTIALS)

fig = plt.figure(figsize=(13, 8))
ax = fig.add_subplot(111, projection='3d')

# Obe plosci v YZ, narisani kot pravokotnika
YY, ZZ = np.meshgrid([-PLATE_Y/2, PLATE_Y/2],
                     [-PLATE_Z/2, PLATE_Z/2])
for xpos, color, lab in [(-GAP/2, 'tomato', '+V/2'),
                          (+GAP/2, 'royalblue', '-V/2')]:
    ax.plot_surface(np.full_like(YY, xpos), YY, ZZ,
                    color=color, alpha=.36, shade=False)
    ax.text(xpos, PLATE_Y/2, PLATE_Z/2, lab, fontsize=11)

if N_EQUIPOTENTIALS > 0:
    cmap = plt.get_cmap('coolwarm')
    norm = colors.Normalize(vmin=-VOLTAGE/2, vmax=VOLTAGE/2)
    spacing = (xgrid[1]-xgrid[0], ygrid[1]-ygrid[0], zgrid[1]-zgrid[0])
    for level in levels:
        if not (potential.min() < level < potential.max()):
            continue
        verts, faces, _, _ = marching_cubes(potential, level=level,
                                             spacing=spacing)
        verts += np.array([xgrid[0], ygrid[0], zgrid[0]])
        ax.plot_trisurf(verts[:,0], verts[:,1], verts[:,2],
                        triangles=faces, color=cmap(norm(level)),
                        alpha=.12, linewidth=0, shade=False)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, shrink=.57, pad=.09, label='Potencial [V]')

for line in lines:
    ax.plot(line[:,0], line[:,1], line[:,2], color='black', linewidth=1.0)

ax.set(xlabel='X [m] (razmik)', ylabel='Y [m]', zlabel='Z [m]',
       title=f'Elektricne silnice in {N_EQUIPOTENTIALS} ekvipotencialnih ploskev')
ax.set_xlim(-.018, .018)
ax.set_ylim(ylim)
ax.set_zlim(zlim)
ax.set_box_aspect((.036, ylim[1]-ylim[0], zlim[1]-zlim[0]))
ax.view_init(elev=23, azim=-55)
plt.tight_layout()
plt.savefig('parallel_plates_yz.png', dpi=160)
plt.show()

# Za AR: tocke v metrih, potenciali ekvipotencialnih ploskev
payload = {
    'units': {'position':'m', 'electric_field':'V/m', 'potential':'V'},
    'geometry': {'normal_axis':'x', 'plate_y':PLATE_Y,
                 'plate_z':PLATE_Z, 'gap':GAP},
    'voltage': VOLTAGE,
    'equipotential_levels': levels.tolist() if N_EQUIPOTENTIALS else [],
    'field_lines': [[{'x':float(p[0]),'y':float(p[1]),'z':float(p[2])}
                      for p in line] for line in lines]
}
with open(OUT_JSON, 'w', encoding='utf-8') as f:
    json.dump(payload, f, indent=2)
print('Datoteki: parallel_plates_yz.png,', OUT_JSON)
print(f'Polje v sredini: {np.linalg.norm(E([0,0,0])):.1f} V/m; idealno U/d: {VOLTAGE/GAP:.1f} V/m')
