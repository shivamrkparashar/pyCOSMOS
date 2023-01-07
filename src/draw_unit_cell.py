import numpy as np
import pandas as pd
import plotly.express as px

import plotly.io as pio
pio.renderers.default ="browser"

Lx, Ly, Lz = 1, 1, 1

def draw_unit_cell():

    A = np.array([0, 0, 0])
    B = np.array([Lx, 0, 0])
    C = np.array([0, Ly, 0])
    D = np.array([0, 0, Lz])

    x = np.array([A[0], B[0]])
    y = np.array([A[1], B[1]])
    z = np.array([A[2], B[2]])


    df = pd.DataFrame({"x": x, "y":y, "z":z})
    fig = px.line_3d(df, x="x", y="y", z="z")
    #px.line_3d(df, x="x", y="y", z="z")

    x = np.array([C[0], B[0]])
    y = np.array([C[1], B[1]])
    z = np.array([C[2], B[2]])
    df = pd.DataFrame({"x": x, "y":y, "z":z})
    #fig = px.line_3d(df, x="x", y="y", z="z")
    fig = px.line_3d(df, x="x", y="y", z="z")
    #fig = go.Figure(data=go.Scatter3d(x=x, y=y,z=z, mode='lines'))

    #xc, yc, zc = [0.5], [0.5], [0.5]
    #fig.add_scatter(df, x=xc, y=yc, z=zc)

    fig.show()

draw_unit_cell()