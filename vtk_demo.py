import os
import tempfile
import numpy as np
from vtk_export import vtk_export

ftype = 'ASCII' # Can be BINARY or ASCII

tmpdir = tempfile.gettempdir()

if True:
    '''
    # vtk DataFile Version 3.0
    Unstructured Grid - Pyramid
    ASCII
    DATASET UNSTRUCTURED_GRID
    POINTS 5 float
    0.0 0.0 0.0
    1.0 0.0 0.0
    1.0 1.0 0.0
    0.0 1.0 0.0
    0.5 0.2 1.0
    '''

    out_filename = os.path.join(tmpdir,'unstructured_grid_demo.vtk')
    points = np.array( [[0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [1.0, 1.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.5, 0.2, 1.0]])

    vtk_export(out_filename, points,
                        title='Unstructured Grid - Pyramid',
                        ftype=ftype)


if True:
    '''
    # vtk DataFile Version 3.0
    Unstructured Grid - Pyramid
    ASCII
    DATASET UNSTRUCTURED_GRID
    POINTS 5 float
    0.0 0.0 0.0
    1.0 0.0 0.0
    1.0 1.0 0.0
    0.0 1.0 0.0
    0.5 0.2 1.0
    CELLS 1 6
    5 0 1 2 3 4
    CELL_TYPES 1
    14
    '''

    out_filename = os.path.join(tmpdir,'unstructured_grid_cells_demo.vtk')

    points = np.array( [[0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [1.0, 1.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.5, 0.2, 1.0]])

    connectivity = {}

    vtk_export(out_filename, points,
                        dataset = 'UNSTRUCTURED_GRID',
                        connectivity = connectivity,
                        point_data = None,
                        point_data_name = 'point data',
                        title='Title',
                        ftype=ftype)

if True:
    '''
    # vtk DataFile Version 3.0
    Structured Grid Example C
    ASCII
    DATASET STRUCTURED_GRID
    DIMENSIONS 3 2 1
    POINTS 6 float
    0 0 0
    1 0 0
    2 0 0
    0 1 0
    1 0.5 0
    2 0.2 0
    '''

    out_filename = os.path.join(tmpdir,'structured_grid_demo.vtk')
    points = np.array( [[0, 0, 0],
                        [1, 0, 0],
                        [2, 0, 0],
                        [0, 1, 0],
                        [1, 0.5, 0],
                        [2, 0.2, 0]] )
    Nx, Ny, Nz = 3, 2, 1
    vtk_export(out_filename, points,
                        dataset = 'STRUCTURED_GRID',
                        connectivity = {'DIMENSIONS' : (Nx,Ny,Nz)},
                        title='Structured Grid Example C',
                        ftype=ftype)

if True:
    '''
    # vtk DataFile Version 3.0
    Dataset attribute example
    ASCII
    DATASET STRUCTURED_GRID
    DIMENSIONS 2 2 1
    POINTS 4 float
    0 0 0
    1 0 0
    0 1 0
    1 1 0
    POINT_DATA 4
    SCALARS point_scalars int 1
    LOOKUP_TABLE default
    0
    0
    3
    3
    '''

    out_filename = os.path.join(tmpdir,'structured_grid_scalars_demo.vtk')
    points = np.array( [[0, 0, 0],
                        [1, 0, 0],
                        [0, 1, 0],
                        [1, 1, 0]] )
    point_data = np.array([0, 0, 3, 3])
    Nx, Ny, Nz = 2, 2, 1
    vtk_export(out_filename, points,
                    dataset = 'STRUCTURED_GRID',
                    connectivity = {'DIMENSIONS' : (Nx,Ny,Nz)},
                    point_data = point_data,
                    texture = 'SCALARS',
                    point_data_name = 'point_scalars',
                    title='Dataset attribute example',
                    ftype=ftype)

if True:
    '''
    # vtk DataFile Version 3.0
    Title
    ASCII
    DATASET POLYDATA
    POINTS 13 float
    0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00
    1.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00
    2.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00
    2.000000000000000000e+00 1.000000000000000000e+00 0.000000000000000000e+00
    3.000000000000000000e+00 2.000000000000000000e+00 0.000000000000000000e+00
    0.000000000000000000e+00 0.000000000000000000e+00 0.000000000000000000e+00
    0.000000000000000000e+00 0.000000000000000000e+00 1.000000000000000000e+00
    0.000000000000000000e+00 0.000000000000000000e+00 2.000000000000000000e+00
    0.000000000000000000e+00 1.000000000000000000e+00 2.000000000000000000e+00
    0.000000000000000000e+00 2.000000000000000000e+00 3.000000000000000000e+00
    0.000000000000000000e+00 -1.000000000000000000e+00 0.000000000000000000e+00
    0.000000000000000000e+00 -2.000000000000000000e+00 0.000000000000000000e+00
    0.000000000000000000e+00 -3.000000000000000000e+00 0.000000000000000000e+00
    LINES 3 16
    5
    0
    1
    2
    3
    4
    5
    5
    6
    7
    8
    9
    3
    10
    11
    12
    '''

    out_filename = os.path.join(tmpdir,'polydata_lines_demo.vtk')

    points = np.array( [[0, 0, 0],
                        [1, 0, 0],
                        [2, 0, 0],
                        [2, 1, 0],
                        [3, 2, 0],
                        [0, 0, 0],
                        [0, 0, 1],
                        [0, 0, 2],
                        [0, 1, 2],
                        [0, 2, 3],
                        [0, -1, 0],
                        [0, -2, 0],
                        [0, -3, 0]])

    connectivity = {'LINES' : np.array([5,5,3])}

    vtk_export(out_filename, points,
                    dataset = 'POLYDATA',
                    connectivity = connectivity,
                    title='Title',
                    ftype=ftype)