import numpy as np
import networkx as nx
import trimesh
from utils import *
import scipy.sparse
import time

# --------------------------------- Load Mesh -------------------------------- #
mesh = trimesh.load('./meshes/bunny.ply', process=False)
export_fn = './exports/bunny_%s.ply' % (time.strftime("%Y%m%d%H%M%S", time.localtime()) )

# ---------------------------- Editing Parameters ---------------------------- #
handles = {
    24092 : [-0.01134  ,  0.151374 , -0.0242688]
}
boundary_ctrl_points = [15617, 24120, 30216, 11236, 6973]


# -------------------------------- Build Graph ------------------------------- #
g = mesh_to_graph(mesh)
boundary = get_boundary(g, boundary_ctrl_points)
editable_verts = get_editable_vertices(g, boundary, list(handles.keys())[0])


# --------------------- Subgraph of the Editiable Region --------------------- #
subgraph = g.subgraph(boundary + editable_verts)
g2l = {}
for n in subgraph.nodes:
    g2l[n] = len(g2l)
l2g = list(subgraph.nodes)
def get_local_neighbor(subgraph, n, l2g, g2l):
    nb = []
    for i in subgraph.neighbors(l2g[n]):
        nb.append(g2l[i])
    return nb


# -------------------------- Build the Linear System ------------------------- #
L = rw_laplacian_matrix(subgraph).todense()
V = np.matrix([subgraph.nodes[i]['pos'] for i in subgraph.nodes])
Delta = L.dot(V)
n = L.shape[0]

LS = np.zeros([3*n, 3*n])
LS[0*n:1*n, 0*n:1*n] = (-1) * L
LS[1*n:2*n, 1*n:2*n] = (-1) * L
LS[2*n:3*n, 2*n:3*n] = (-1) * L

for i in range(n):
    nb_idx = get_local_neighbor(subgraph, i, l2g, g2l)
    ring = np.array([i] + nb_idx)
    V_ring = V[ring]
    n_ring = V_ring.shape[0]
    
    A = np.zeros([n_ring * 3, 7])
    for j in range(n_ring):
        A[j]          = [V_ring[j,0],           0 ,   V_ring[j,2], -V_ring[j,1], 1, 0, 0]
        A[j+n_ring]   = [V_ring[j,1], -V_ring[j,2],            0 ,  V_ring[j,0], 0, 1, 0]
        A[j+2*n_ring] = [V_ring[j,2],  V_ring[j,1], -V_ring[j, 0],           0 , 0, 0, 1]
        
    # Moore-Penrose Inversion
    A_pinv = np.linalg.pinv(A)
    s = A_pinv[0]
    h = A_pinv[1:4]
    t = A_pinv[4:7]
    

    T_delta = np.vstack([
         Delta[i,0]*s    - Delta[i,1]*h[2] + Delta[i,2]*h[1],
         Delta[i,0]*h[2] + Delta[i,1]*s    - Delta[i,2]*h[0],
        -Delta[i,0]*h[1] + Delta[i,1]*h[0] + Delta[i,2]*s   ,
    ])
        
    LS[i, np.hstack([ring, ring+n, ring+2*n])] += T_delta[0]
    LS[i+n, np.hstack([ring, ring+n, ring+2*n])] += T_delta[1]
    LS[i+2*n, np.hstack([ring, ring+n, ring+2*n])] += T_delta[2]


# ------------------- Add Constraints to the Linear System ------------------- #
constraint_coef = []
constraint_b = []

# Boundary constraints
boundary_idx = [g2l[i] for i in boundary_ctrl_points]
for idx in boundary_idx:
    constraint_coef.append(np.arange(3*n) == idx)
    constraint_coef.append(np.arange(3*n) == idx + n)
    constraint_coef.append(np.arange(3*n) == idx + 2*n)
    constraint_b.append(V[idx, 0])
    constraint_b.append(V[idx, 1])
    constraint_b.append(V[idx, 2])
    
# Handle constraints
for gid, pos in handles.items():
    idx = g2l[gid]
    constraint_coef.append(np.arange(3*n) == idx)
    constraint_coef.append(np.arange(3*n) == idx + n)
    constraint_coef.append(np.arange(3*n) == idx + 2*n)
    constraint_b.append(pos[0])
    constraint_b.append(pos[1])
    constraint_b.append(pos[2])
    
constraint_coef = np.matrix(constraint_coef)
constraint_b = np.array(constraint_b)


# -------------------------- Solve the Linear System ------------------------- #
A = np.vstack([LS, constraint_coef])
b = np.hstack([np.zeros(3*n), constraint_b])
spA = scipy.sparse.coo_matrix(A)

V_prime = scipy.sparse.linalg.lsqr(spA, b)


# -------------------------- Output the Edited Mesh -------------------------- #
new_pnts = []
for i in range(n):
    new_pnts.append([V_prime[0][i], V_prime[0][i+n], V_prime[0][i+2*n]])
    
new_mesh = mesh.copy()
for idx, pnt in enumerate(new_pnts):
    gid = l2g[idx]
    new_mesh.vertices[gid] = pnt

new_mesh.export(export_fn)