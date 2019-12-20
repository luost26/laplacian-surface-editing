import numpy as np
import networkx as nx
import trimesh
import scipy.sparse


def mesh_to_graph(mesh):
    edges = mesh.edges_unique
    verts = mesh.vertices
    g = nx.Graph()
    for i, pos in enumerate(verts):
        g.add_node(i, pos=pos)
    for edge in edges:
        g.add_edge(*edge)
    return g


def get_boundary(g, ctrl_points):
    boundary = []
    for i in range(len(ctrl_points)):
        p1 = ctrl_points[i]
        p2 = ctrl_points[(i+1)%len(ctrl_points)]
        path = nx.shortest_path(g, p1, p2)
        boundary += path[0:-1]
    return boundary


def rw_laplacian_matrix(g):
    lap = nx.laplacian_matrix(g)
    adj = nx.adjacency_matrix(g)
    inv_deg = scipy.sparse.diags(1/adj.dot(np.ones([adj.shape[0]])))
    return inv_deg.dot(lap)


def get_editable_vertices(g, boundary, manip_handle):
    g_part = g.copy()
    g_part.remove_nodes_from(boundary)
    return list(nx.node_connected_component(g_part, manip_handle))


def get_mesh_scene(mesh, boundary=None, editable_vertices=None):
    scene = [mesh]
    if boundary:
        scene.append(trimesh.load_path(mesh.vertices[boundary + [boundary[0]]]))
    if editable_vertices:
        scene.append(trimesh.points.PointCloud(mesh.vertices[editable_vertices + boundary]))
    
    return trimesh.Scene(scene)
    