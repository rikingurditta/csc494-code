import numpy as np


def make_disc(resolution):
    angles = np.arange(0., 2. * np.pi, 2. * np.pi / resolution)
    x_outer = np.cos(angles)
    y_outer = np.sin(angles)
    x_inner = 0.5 * np.cos(-angles)
    y_inner = 0.5 * np.sin(-angles)
    vertices = []
    edges = []
    for j in range(resolution):
        vertices.append(f'{x_outer[j]},{y_outer[j]}')
        edges.append(f'{j},{(j + 1) % resolution}')
    for j in range(resolution):
        vertices.append(f'{x_inner[j]},{y_inner[j]}')
        edges.append(f'{resolution + j},{resolution + (j + 1) % resolution}')
    return vertices, edges


def make_two_hole_disc(resolution):
    angles = np.arange(0., 2. * np.pi, 2. * np.pi / resolution)
    x_outer = np.cos(angles)
    y_outer = np.sin(angles)
    x_inner_1 = 0.5 * np.cos(-angles)
    y_inner_1 = 0.5 * np.sin(-angles)
    x_inner_2 = 0.2 * np.cos(-angles) + 0.75
    y_inner_2 = 0.2 * np.sin(-angles)
    vertices = []
    edges = []
    for j in range(resolution):
        vertices.append(f'{x_outer[j]},{y_outer[j]}')
        edges.append(f'{j},{(j + 1) % resolution}')
    for j in range(resolution):
        vertices.append(f'{x_inner_1[j]},{y_inner_1[j]}')
        edges.append(f'{resolution + j},{resolution + (j + 1) % resolution}')
    for j in range(resolution):
        vertices.append(f'{x_inner_2[j]},{y_inner_2[j]}')
        edges.append(f'{resolution + resolution + j},{resolution + resolution + (j + 1) % resolution}')
    return vertices, edges


def write_v_e_files(basename, resolution, vertices, edges):
    with open(f'{basename}_{resolution}_V.csv', 'w') as v_file:
        v_file.write('\n'.join(vertices))
    with open(f'{basename}_{resolution}_E.csv', 'w') as e_file:
        e_file.write('\n'.join(edges))


if __name__ == '__main__':
    basename = 'disc'
    resolutions = [256, 192, 128, 96, 64, 48, 32, 16, 8, 5, 3]
    for res in resolutions:
        vs, es = make_disc(res)
        write_v_e_files(basename, res, vs, es)
    vs, es = make_two_hole_disc(200)
    write_v_e_files('two_hole_right_small_disc', 200, vs, es)
