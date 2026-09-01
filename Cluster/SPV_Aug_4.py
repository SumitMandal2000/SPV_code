#1. Import libraries:
import numpy as np
import time
from scipy.spatial import Voronoi, Delaunay, voronoi_plot_2d
from matplotlib import cm
from collections import defaultdict
from tqdm import trange
import random
import matplotlib.pyplot as plt
import os

#2. Define Functions:
#2.1: Random sequantial Addition to produce initial cells:
def random_sequential_addition(n_points, domain_size, min_distance):
    points = []
    iterations = 0
    max_iterations=1000*n_points
    while len(points) < n_points and iterations < max_iterations:
        candidate = np.random.rand(2) * domain_size
        if all(np.linalg.norm(candidate - p) > min_distance for p in points):
            points.append(candidate)
        iterations += 1

    success = len(points) == n_points
    if not success:
        print("Warning: Reached maximum iterations (",max_iterations,") before generating all points.")
    else:
        print("The number of points generated is ",len(points))

    return np.array(points)


# 2.2.Periodic boundary cpoints: To creat extra 8 copies of original cells coordinate
def peridodic_points(pts, L):
  Lx, Ly = L
  n = pts.shape[0]
  values = np.array([-1, 0, 1])
  offsets = np.array([[x, y] for x in values for y in values if not (x == 0 and y == 0)])
  offsets[:, 0] *= Lx
  offsets[:, 1] *= Ly
  new_pts = np.empty((n*9, pts.shape[1]))
  new_pts[:n] = pts
  for i in range(8):
    new_pts[n * (i + 1):n * (i + 2)] = pts + offsets[i]
  return new_pts


#2.3. function to get vertex coodinates for each triangles whlie preserving periodicity
def get_vertices(cent):
  ri, rj, rk = cent.transpose(1, 2, 0)
  ax, ay = ri
  bx, by = rj
  cx, cy = rk
  d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
  ux = ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (
          ay - by)) / d
  uy = ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (
          bx - ax)) / d
  vs = np.empty((ux.shape[0], 2))
  vs[:, 0] = ux
  vs[:, 1] = uy
  return vs

#2.4. Cell to triangles: find the relation between cells and triangles
def cell_to_triangles(tri, n_c):
    cell_triangles = defaultdict(list)

    # Loop through each triangle and each vertex in the triangle
    for idx, triangle in enumerate(tri):
        for vertex in triangle:
            if vertex < n_c:  # Only consider cells below n_c
                cell_triangles[vertex].append(idx)

    return cell_triangles

#2.5. find the vertices for each cells and order them in CCW (both indices and coordinates)
def cell_to_vertices(vertices, cell_triangles, pts, n_c):
    cell_vertices = []
    for i in range(n_c):
        v = vertices[cell_triangles[i]]
        center = pts[i]
        angles = np.arctan2(v[:, 1] - center[1], v[:, 0] - center[0])
        sorted_indices = np.argsort(angles)
        cell_vertices.append(v[sorted_indices])
        cell_triangles[i] = np.array(cell_triangles[i])[sorted_indices]
    return cell_vertices

#2.6. def a function to get all the needed triangles (sort the needed triangle to inside the box + neighbors)
def get_triangles(tri, cell_triangles, f_tri, f_nei, n_c):
  #create an arrray that contains all trianges in the needed order
  triangles = []
  neighbors = []
  for i in range(n_c):
      arr1 = f_tri[cell_triangles[i]]
      arr2 = f_nei[cell_triangles[i]]
      ids = []
      for row1, row2 in zip(arr1, arr2):
        idx = np.where(row1 == i)[0][0]
        ids.append(row2[idx])
        #neighbors.append(tri[row2[idx]])
        if idx != 0 :
          row1[:] = np.roll(row1, -idx)
      triangles.append(np.array(arr1))
      neighbors.append(tri[ids])
  return triangles, neighbors

#2.7. function to calculate area and perimeter for each cells
def calculate_area_perimeter(cell_vertices, n_c):
    area_array = []
    perimeter_array = []

    for i in range(n_c):
      # Perimeter Calculation
      perimeter = np.sum(np.linalg.norm(cell_vertices[i] - np.roll(cell_vertices[i], -1, axis=0), axis=1))
      perimeter_array.append(perimeter)

      # Area Calculation
      area = 0.5 * np.abs(np.dot(cell_vertices[i][:, 0], np.roll(cell_vertices[i][:, 1], -1)) -
                          np.dot(cell_vertices[i][:, 1], np.roll(cell_vertices[i][:, 0], -1)))
      area_array.append(area)

    return np.array(area_array), np.array(perimeter_array)
#2.8. Dot product
def dot(a, b):
    return np.sum(a * b, axis=1)

# 2.9. Vertex derivative:
#function to get vertex derivatives wrt cell center for all cells
def get_dHdr(triangles, p, f_tri, n_c):
  tr = [item for sublist in triangles for item in sublist]
  rijk = p[tr]
  ri = rijk[:, 0]
  rj = rijk[:, 1]
  rk = rijk[:, 2]
  # print(rj[:6])
  r_ik = rk - ri
  r_ij = rj - ri
  r_jk = rk - rj
  c = r_ij[:, 0]*r_jk[:, 1] - r_ij[:, 1]*r_jk[:, 0]
  d = 2*c**2
  bd = - (dot(r_ik, r_ik)) * dot(r_ij, r_jk)
  cd = (dot(r_ij, r_ij)) * dot(r_ik, r_jk)
  z = bd[:, np.newaxis]*r_ij + cd[:, np.newaxis]*r_ik
  del_bd = 2*(dot(r_ij, r_jk)[:, np.newaxis])*r_ik + r_jk*(dot(r_ik, r_ik)[:, np.newaxis])
  del_cd = -2*(dot(r_ik, r_jk)[:, np.newaxis])*r_ij - r_jk*(dot(r_ij, r_ij)[:, np.newaxis])
  n = ri.shape[0]
  del_d = np.column_stack([-(2/c) * r_jk[:, 1], (2/c) * r_jk[:, 0]])
  In = np.array([np.eye(2) for _ in range(n)])
  del_h_ijk = In + (np.einsum('ij, ik -> ijk', r_ij, del_bd) + np.einsum('ij, ik -> ijk', r_ik, del_cd) - In * (cd + bd)[:, np.newaxis, np.newaxis] - np.einsum('ij, ik -> ijk', z, del_d)) / d[:, np.newaxis, np.newaxis]

  return np.transpose(del_h_ijk, (0,2,1))

#define a dunction to calculate energy derivative wrt to vertex for all cells
def self_dEdh(v1, v2, v3, A, P, A0, p0, kp, cell_triangles, cell_vertices, f_tri, n_c, n_vertices_per_cell, cumulative_vertices):
  a = np.repeat(A/A0, n_vertices_per_cell)
  p = np.repeat(P/np.sqrt(A0), n_vertices_per_cell)

  # Adjust indices to account for cell boundaries
  for start, end in zip(cumulative_vertices[:-1], cumulative_vertices[1:]):
      v2[end - 1] = v1[start]  # Ensure v2 wraps around within the same cell
      v3[start] = v1[end - 1]  # Ensure v3 wraps around within the same cell
  #print(f'v1, v2, v3={v1, v2, v3}')
  d1 = np.linalg.norm(v1 - v3, axis=1)
  d2 = np.linalg.norm(v1 - v2, axis=1)

  Ax = (a - 1) * (v2[:, 1] - v3[:, 1]) / np.sqrt(A0)
  Px = 2 * kp * (p - p0) * ((v1[:, 0] - v3[:, 0]) / d1 + (v1[:, 0] - v2[:, 0]) / d2)
  Fx = Ax + Px

  Ay = (a - 1) * (v3[:, 0] - v2[:, 0]) / np.sqrt(A0)
  Py = 2 * kp * (p - p0) * ((v1[:, 1] - v3[:, 1]) / d1 + (v1[:, 1] - v2[:, 1]) / d2)
  Fy = Ay + Py

  return np.column_stack((Fx, Fy))


#function to find self force acting on each cell
def self_force(dEdh, dHdr, n_c, cumulative_vertices):

    # Compute dot products for all cells in a vectorized manner
    Fx = np.sum(dEdh * dHdr[:, 0], axis=1)  # Force in x direction
    Fy = np.sum(dEdh * dHdr[:, 1], axis=1)  # Force in y direction

    # Use slicing and np.sum to compute the sum of forces for each cell
    fx = -np.array([np.sum(Fx[cumulative_vertices[i]:cumulative_vertices[i + 1]]) for i in range(n_c)])
    fy = -np.array([np.sum(Fy[cumulative_vertices[i]:cumulative_vertices[i + 1]]) for i in range(n_c)])

    # Combine the force components into a single array
    return np.column_stack((fx, fy))

#function to get the necessary third vertex in counter clockwise direction
def get_n_vertices(neighbors, p, n_c):
    n_cell_vertices = []
    for i in range(n_c):
        n_cell_vertices.append(get_vertices(p[neighbors[i]]))
    return n_cell_vertices

#define a dunction to calculate energy derivative wrt to vertex
def neighbor_dEdh(v1, v2, v3, triangles, A, P, A0, kp, p0, f_tri, tri, cell_vertices, n_cell_vertices, n_c, cumulative_vertices):

  nv = np.concatenate(n_cell_vertices)
  # Adjust indices to account for cell boundaries
  for start, end in zip(cumulative_vertices[:-1], cumulative_vertices[1:]):
      v2[end - 1] = v1[start]  # Ensure v2 wraps around within the same cell
      v3[start] = v1[end - 1]  # Ensure v3 wraps around within the same cell
  triangles = [item for sublist in triangles for item in sublist]
  # print(f'ri : {np.array(triangles)[:, 0]}')
  # print(f'rj : {np.array(triangles)[:, 1]}')
  # print(f'rk : {np.array(triangles)[:, 2]}')
  rj = np.array(triangles)[:, 1] % n_c
  # print('Area:',A[rj][:6])
  # print('Perimeter:',P[rj][:6])
  #E1
  d1f = np.linalg.norm(v1 - v2, axis=1)
  d2f = np.linalg.norm(v1 - nv, axis=1)
  cf = np.array(triangles)[:, 2].ravel()  %n_c
  af = A[cf]/A0
  pf = P[cf]/np.sqrt(A0)
  Axf = (af - 1) * (nv[:,1] - v2[:,1]) / np.sqrt(A0)
  Pxf = 2 * kp * (pf - p0) * ((v1[:,0] - v2[:,0])/d1f + ((v1[:,0] - nv[:,0])/d2f))
  Fxf = Axf + Pxf
  Ayf = (af - 1) * (v2[:,0] - nv[:,0]) / np.sqrt(A0)
  Pyf = 2 * kp * (pf - p0) * ((v1[:,1] - v2[:,1])/d1f + ((v1[:,1] - nv[:,1])/d2f))
  Fyf = Ayf + Pyf
  E1 = np.concatenate((Fxf[:, np.newaxis], Fyf[:, np.newaxis]), axis=1)
  #for i in range(6):
  #    print(f'v1 = {v1[i]}\tv2 = {nv[i]}\tv3 = {v2[i]}')
  #E2
  d1b = np.linalg.norm(v1 - nv, axis=1)
  d2b = np.linalg.norm(v1 - v3, axis=1)
  cb = np.array(triangles)[:, 1].ravel()  %n_c
  ab = A[cb]/A0
  pb = P[cb]/np.sqrt(A0)
  Axb = (ab - 1) * (v3[:,1] - nv[:,1]) / np.sqrt(A0)
  Pxb = 2 * kp * (pb - p0) * ((v1[:,0] - nv[:,0])/d1b + ((v1[:,0] - v3[:,0])/d2b))
  Fxb = Axb + Pxb
  Ayb = (ab - 1) * (nv[:,0] - v3[:,0]) / np.sqrt(A0)
  Pyb = 2 * kp * (pb - p0) * ((v1[:,1] - nv[:,1])/d1b + ((v1[:,1] - v3[:,1])/d2b))
  Fyb = Ayb + Pyb
  E2 = np.concatenate((Fxb[:, np.newaxis], Fyb[:, np.newaxis]), axis=1)
  #for i in range(6):
  #    print(f'v1 = {v1[i]}\tv2 = {v3[i]}\tv3 = {nv[i]}')
  return np.array(E1 + E2)

#function to find neighbor force acting on each cell
def neighbor_force(dEdh, dHdr, n_c, cumulative_vertices):
    # Compute dot products for all cells in a vectorized manner
    Fx = np.sum(dEdh * dHdr[:, 0], axis=1)  # Force in x direction
    Fy = np.sum(dEdh * dHdr[:, 1], axis=1)  # Force in y direction

    # Use slicing and np.sum to compute the sum of forces for each cell
    fx = -np.array([np.sum(Fx[cumulative_vertices[i]:cumulative_vertices[i + 1]]) for i in range(n_c)])
    fy = -np.array([np.sum(Fy[cumulative_vertices[i]:cumulative_vertices[i + 1]]) for i in range(n_c)])

    # Combine the force components into a single array
    return np.column_stack((fx, fy))

def generate_square_lattice_with_noise(rows, cols, spacing=1.0, noise_amplitude=0.1):
    np.random.seed(10) # Need to remove it
    # Generate the square lattice
    x_coords = np.arange(0, cols * spacing, spacing)
    y_coords = np.arange(0, rows * spacing, spacing)
    x_grid, y_grid = np.meshgrid(x_coords, y_coords)
    lattice_points = np.column_stack([x_grid.ravel(), y_grid.ravel()])

    # Add random noise
    noise = np.random.uniform(-noise_amplitude, noise_amplitude, size=lattice_points.shape)
    displaced_points = lattice_points + noise

    return displaced_points

def generate_hexagonal_lattice_inside_box(rows, cols, spacing=1.0, noise_amplitude=0.01):
    # Horizontal and vertical offsets for the hexagonal grid
    dx = spacing  # Horizontal distance between points
    dy = spacing  # Vertical distance between points

    # Generate lattice points inside the box
    x_coords = []
    y_coords = []

    # To ensure points are within the bounding box, we need to adjust for the staggered grid
    for row in range(rows):
        for col in range(cols):
            x = col * dx
            y = row * dy
            if col % 2 == 1:  # Offset every other column to create hexagonal tiling
                y += dy / 2
            # Add coordinates only if they are inside the bounding box
            if x < cols * dx and y < rows * dy:
                x_coords.append(x)
                y_coords.append(y)

    lattice_points = np.column_stack([x_coords, y_coords])

    # Add random noise
    noise = np.random.uniform(-noise_amplitude, noise_amplitude, size=lattice_points.shape)
    displaced_points = lattice_points + noise

    return displaced_points



def Q_calculation(data_t, data_0, L, l_dash=0.64):
    dr = data_t - data_0
    dr_mag = np.linalg.norm(dr, axis=1)
    wi = (dr_mag <= l_dash).astype(float)
    Q_t = np.mean(wi)

    return Q_t

def uniform_pin_indices(pts, L, N_pin):
    '''
    Greedy Farthest Point Sampling (FPS):
    Every new pinned cell is chosen to be as far away as possible from all previously selected pinned
    cells.
    '''
    if N_pin==0:
        return np.array([], dtype=int)

    N = len(pts)
    pin_index = [np.random.randint(N)]                    # Choose the first pinned cell randomly.
    min_dist = np.full(N, np.inf)                         # min_dist[i] stores the distance of cell i to its nearest pinned cell.
                                                          # Initially, no distances are known, so set them to infinity.
    for _ in range(N_pin - 1):        
        last = pts[pin_index[-1]]                         # Coordinates of the most recently selected pinned cell.
        dr = pts - last                                   # Displacement from the last pinned cell to every cell.      
        dr -= L * np.round(dr / L)                        # Apply the minimum-image convention for periodic boundaries.
        dist = np.linalg.norm(dr, axis=1)                 # Distance from the last pinned cell to every cell.
        min_dist = np.minimum(min_dist, dist)             # It will choose the minimum distance between min_dis and dist
        min_dist[pin_index] = -1                          # Prevent already selected cells from being chosen again.
        pin_index.append(np.argmax(min_dist))             # Select the cell whose nearest pinned neighbour is farthest away.
    return np.array(pin_index)


#define a function for simulation
def Simulation(sample_pts, p0, v0, Dr, pinning_percent, L, h, init_steps,initial_steps_pin,  n_steps,msd_step, set_i, save_coordinate=False):
    kp = 1
    A0 = 1
    l_dash=0.60
    N = sample_pts.shape[0]
    N_pin= int((pinning_percent/100)*N)
    T = h*n_steps
    
    duplicate_pts = sample_pts.copy()
    theta = np.random.uniform(0, 2 * np.pi, N)
    times = np.linspace(0, T, n_steps + 1)

    stop_exp=int(np.log10(n_steps+0.1*n_steps))
    log_steps = np.logspace(0, stop_exp, num=400, base=10)
    save_steps = np.unique(np.concatenate(([0], np.round(log_steps).astype(int))))
    
    start_time = time.time()
    Time= np.array(save_steps)*h
    Coordinates = []
    Velocity = []
    Theta =[]
    t_Force_Equilibrium = []
    Force_Equilibrium=[]
    folder = "Data"
    save_F_motility=True
    os.makedirs(folder, exist_ok=True)

    ## 1. Equilibriation with zero motility:
    for step in range(init_steps+1):
        pts = peridodic_points(duplicate_pts, L)
        Tri = Delaunay(pts)
        tri = Tri.simplices
        nei = Tri.neighbors
        mask = np.any(tri < N, axis=1)
        f_tri = tri[mask]
        f_nei = nei[mask]

        cent = pts[f_tri]
        vertices = get_vertices(cent)
        cell_triangles = cell_to_triangles(f_tri, N)
        cell_vertices = cell_to_vertices(vertices, cell_triangles, pts, N)

        A, P = calculate_area_perimeter(cell_vertices, N)
        
        triangles, neighbors = get_triangles(tri, cell_triangles, f_tri, f_nei, N)
        n_vertices_per_cell = [len(v) for v in cell_vertices]
        cumulative_vertices = np.cumsum(np.concatenate(([0],n_vertices_per_cell)))
        v1 = np.concatenate(cell_vertices, axis=0)
        v2 = np.roll(v1, -1, axis=0)
        v3 = np.roll(v1, 1, axis=0)
        dHdr = get_dHdr(triangles, pts, f_tri, N)
        dEdh_s = self_dEdh(v1, v2, v3, A, P, A0, p0, kp, cell_triangles, cell_vertices, f_tri, N, n_vertices_per_cell, cumulative_vertices)
        F_self = self_force(dEdh_s, dHdr, N, cumulative_vertices)
        n_cell_vertices = get_n_vertices(neighbors, pts, N)
        dEdh_n = neighbor_dEdh(v1, v2, v3, triangles, A, P, A0, kp, p0, f_tri, tri, cell_vertices, n_cell_vertices, N, cumulative_vertices)
        F_nei = neighbor_force(dEdh_n, dHdr, N, cumulative_vertices)
        
        F = F_self + F_nei
        F-=np.mean(F, axis=0)
        duplicate_pts += h * (F)
        duplicate_pts = np.mod(duplicate_pts, L[0])


    if save_F_motility:
        Data_F_no_motility={}
        Data_F_no_motility={
            'time': step*dt,
            'Force': F
        }
        np.save("%s/F_no_motility_p0=%.2f_v0=%0.2f_pin=%d_set_%d.npy" % (folder, p0, v0, pinning_percent, set_i), Data_F_no_motility)

    ## 2. Equilibriation with zero motility with pin:
    pinned_cells = uniform_pin_indices(duplicate_pts, L, N_pin)
    unpinned_cells = np.setdiff1d(np.arange(N), pinned_cells)
    print('pinned_cells=', pinned_cells, 'number=', len(pinned_cells))
    
    
    for step in trange(initial_steps_pin+1):
        pts = peridodic_points(duplicate_pts, L)
        Tri = Delaunay(pts)
        tri = Tri.simplices
        nei = Tri.neighbors
        mask = np.any(tri < N, axis=1)
        f_tri = tri[mask]
        f_nei = nei[mask]

        cent = pts[f_tri]
        vertices = get_vertices(cent)
        cell_triangles = cell_to_triangles(f_tri, N)
        cell_vertices = cell_to_vertices(vertices, cell_triangles, pts, N)

        A, P = calculate_area_perimeter(cell_vertices, N)
        
        triangles, neighbors = get_triangles(tri, cell_triangles, f_tri, f_nei, N)
        n_vertices_per_cell = [len(v) for v in cell_vertices]
        cumulative_vertices = np.cumsum(np.concatenate(([0],n_vertices_per_cell)))
        v1 = np.concatenate(cell_vertices, axis=0)
        v2 = np.roll(v1, -1, axis=0)
        v3 = np.roll(v1, 1, axis=0)
        dHdr = get_dHdr(triangles, pts, f_tri, N)
        dEdh_s = self_dEdh(v1, v2, v3, A, P, A0, p0, kp, cell_triangles, cell_vertices, f_tri, N, n_vertices_per_cell, cumulative_vertices)
        F_self = self_force(dEdh_s, dHdr, N, cumulative_vertices)
        n_cell_vertices = get_n_vertices(neighbors, pts, N)
        dEdh_n = neighbor_dEdh(v1, v2, v3, triangles, A, P, A0, kp, p0, f_tri, tri, cell_vertices, n_cell_vertices, N, cumulative_vertices)
        F_nei = neighbor_force(dEdh_n, dHdr, N, cumulative_vertices)
        
        F_total = F_self + F_nei
        F_total[pinned_cells]=0
        F_total[unpinned_cells]-=np.mean(F_total[unpinned_cells], axis=0)
        
        duplicate_pts += h * F_total
        
        duplicate_pts = np.mod(duplicate_pts, L[0])

    if save_F_motility:
        Data_F_no_motility={}
        Data_F_no_motility={
            'time': step*dt,
            'Force': F_total
        }
        np.save("%s/F_no_motility_pin_p0=%.2f_v0=%0.2f_pin=%d_set_%d.npy" % (folder, p0, v0, pinning_percent, set_i), Data_F_no_motility)


    
    ## 2. Steady state with motility:
    msd_pts = duplicate_pts.copy()
    pts_zero_motility = duplicate_pts.copy()
    datai0_cm=duplicate_pts.copy()


    file = open('%s/%d_cells_tf=%0.1f_p0=%0.2f_v0=%0.2f_Dr=%0.3f_dt=%0.3f_p=%d_%d.txt'%(folder, N, T, p0, v0, Dr, h,pinning_percent, set_i), 'w')
    file.write('Time \t MSD \t\t Qt \t Energy \t q \n')
    for step in trange(n_steps+msd_step+1):
        pts = peridodic_points(duplicate_pts, L)
        Tri = Delaunay(pts)
        tri = Tri.simplices
        nei = Tri.neighbors
        mask = np.any(tri < N, axis=1)
        f_tri = tri[mask]
        f_nei = nei[mask]

        cent = pts[f_tri]
        vertices = get_vertices(cent)
        cell_triangles = cell_to_triangles(f_tri, N)
        cell_vertices = cell_to_vertices(vertices, cell_triangles, pts, N)

        A, P = calculate_area_perimeter(cell_vertices, N)
        A0 = 1
        triangles, neighbors = get_triangles(tri, cell_triangles, f_tri, f_nei, N)
        n_vertices_per_cell = [len(v) for v in cell_vertices]
        cumulative_vertices = np.cumsum(np.concatenate(([0],n_vertices_per_cell)))
        v1 = np.concatenate(cell_vertices, axis=0)
        v2 = np.roll(v1, -1, axis=0)
        v3 = np.roll(v1, 1, axis=0)
        dHdr = get_dHdr(triangles, pts, f_tri, N)
        dEdh_s = self_dEdh(v1, v2, v3, A, P, A0, p0, kp, cell_triangles, cell_vertices, f_tri, N, n_vertices_per_cell, cumulative_vertices)
        F_self = self_force(dEdh_s, dHdr, N, cumulative_vertices)
        n_cell_vertices = get_n_vertices(neighbors, pts, N)
        dEdh_n = neighbor_dEdh(v1, v2, v3, triangles, A, P, A0, kp, p0, f_tri, tri, cell_vertices, n_cell_vertices, N, cumulative_vertices)
        F_nei = neighbor_force(dEdh_n, dHdr, N, cumulative_vertices)


        F = F_self + F_nei
        n = np.zeros_like(F)
        eta=np.random.normal(loc=0.0, scale=1.0, size=N)
        theta += np.sqrt(2*Dr* h)* eta
        n[:, 0] = np.cos(theta)
        n[:, 1] = np.sin(theta)
        F_total=(F + v0 * n )
        velocity_t = F_total
        step1=step-msd_step

        if step1==0:
            pts_zero_motility=duplicate_pts.copy()
            msd_pts = duplicate_pts.copy()

        F_total[pinned_cells]=0
        F_total[unpinned_cells]-=np.mean(F_total[unpinned_cells], axis=0)
        duplicate_pts += h * F_total
        msd_pts +=h*F_total
        

        #data_inside_0=pts_zero_motility[inside_radius_indices]
        data_0=pts_zero_motility[unpinned_cells]
        
        if step1 in save_steps:
            
            data_t=msd_pts[unpinned_cells]
            q = np.sum(P[unpinned_cells]/np.sqrt(A[unpinned_cells]))
            energy = 0
            for a, p in zip(A[unpinned_cells], P[unpinned_cells]):
                energy +=  (a - 1)**2 + kp * (p - p0)**2

            msd=np.mean(np.sum((data_t-data_0)**2, axis=1))
            Qt = Q_calculation(data_t, data_0, L, l_dash=l_dash)

            ## Storing the data:
            file.write(f'{(step1)*h : 0.3f} \t {msd:0.6f} \t {Qt : 0.3f}\t {energy:0.4f} \t {q/N : 0.3f}\t\n')

        if save_coordinate:
            if step1 in save_steps:
                Coordinates.append(msd_pts.copy())
                Theta.append(theta.copy())
                Velocity.append(velocity_t.copy())

        duplicate_pts = np.mod(duplicate_pts, L[0])
    
    if save_coordinate:
        np.save("%s/pinned_cells_p0=%.2f_v0=%0.2f_pin=%d_set_%d.npy" % (folder, p0, v0, pinning_percent, set_i), pinned_cells)
        np.save("%s/coordinate_p0=%.2f_v0=%0.2f_pin=%d_set_%d.npy" % (folder, p0, v0, pinning_percent, set_i), Coordinates)
        np.save("%s/Theta_p0=%.2f_v0=%0.2f_pin=%d_set_%d.npy" % (folder, p0, v0, pinning_percent, set_i), Theta)
        np.save("%s/Velocity_p0=%.2f_v0=%0.2f_pin=%d_set_%d.npy" % (folder, p0, v0, pinning_percent, set_i), Velocity)
        np.save("%s/time.npy"%(folder), Time)
    print(len(save_steps))
    
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Total time taken: {total_time} seconds")

    return None

#-----------------------------------------------------------------------------------------------------------
l=32
L = np.array([l, l])
pts= random_sequential_addition(L[0]**2, L[0], 0.5)
p0=3.8
v0=0.1
Dr=1
pinning_percent=5
dt = 0.10
initial_steps=100000
initial_steps_pin=10000
n_steps =100000
msd_step=100000
set_i=1
Simulation(pts, p0, v0, Dr,  pinning_percent, L, dt,initial_steps, initial_steps_pin, n_steps,msd_step,  set_i, save_coordinate=True)
