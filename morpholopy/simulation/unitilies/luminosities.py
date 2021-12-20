import numpy as np


#####
L_sun = 3.846e26
pc = 3.08567758e16
A_shell = 4 * np.pi * (10 * pc) ** 2
to_jy = 1e26 * L_sun / A_shell
zp_jy = 3631.0

magfac = to_jy / zp_jy
#####

################################################################################

# interp = 'pow'
# extrap = 'n'

# if interp == 'log':
#     if extrap == 'y':
#         interpfunc  	= BiLogInterp_extrap
#     else:
#         interpfunc	= BiLogInterp

# elif interp == 'pow':
#     if extrap == 'y':
#         interpfunc	= BiPowInterp_extrap
#     else:
#         interpfunc      = BiPowInterp

# else:
#     if extrap == 'y':
#         interpfunc  	= BiLinInterp_extrap
#     else:
#         interpfunc  	= BiLinInterp

################################################################################


# grids  = {}

# for pht in glob.glob(f'./photometry/{system}/*'):
#     grids[pht[-1]] = MakeGrid(pht)


def MakeGrid(bcfile):
    LumConv = np.genfromtxt(bcfile)
    cycles = LumConv[:, 1] == LumConv[0, 1]
    Z_points = LumConv[:, 0][cycles]
    idxs = np.arange(LumConv[:, 0].size)
    rangeedge = idxs[cycles][1]
    age_points = LumConv[:, 1][0:rangeedge]
    fluxs = LumConv[:, -1].copy()
    fluxs.shape = (Z_points.size, rangeedge)
    fluxs, Z_points, age_points = SortFluxGrid(fluxs, Z_points, age_points)
    return fluxs, Z_points, age_points


def SortFluxGrid(fs, Zs, ts):
    Zord = Zs.argsort()
    tord = ts.argsort()
    fs = fs.copy()[Zord, :]
    fs = fs.copy()[:, tord]
    return fs, Zs[Zord], ts[tord]


def BiPowInterp(x_nodes, y_nodes, gridvals, x, y):
    old_settings = np.seterr(all="ignore")

    x_nodes = np.clip(np.log10(x_nodes), -100, 100)
    y_nodes = np.clip(np.log10(y_nodes), -4, 100)

    x = np.clip(np.log10(x), -100, 100)
    y = np.clip(np.log10(y), -4, 100)

    gridvals = np.log10(gridvals)

    x_clip = np.clip(x, x_nodes.min(), (x_nodes.max()))
    y_clip = np.clip(y, y_nodes.min(), (y_nodes.max()))

    x_digi = np.clip(np.digitize(x_clip, x_nodes), -1, x_nodes.size - 1)
    y_digi = np.clip(np.digitize(y_clip, y_nodes), -1, y_nodes.size - 1)

    f_11 = gridvals[y_digi - 1, x_digi - 1]
    f_12 = gridvals[y_digi, x_digi - 1]
    f_21 = gridvals[y_digi - 1, x_digi]
    f_22 = gridvals[y_digi, x_digi]

    xdiff1 = x_nodes[x_digi] - x_clip
    xdiff2 = x_clip - x_nodes[x_digi - 1]
    xdiff3 = x_nodes[x_digi] - x_nodes[x_digi - 1]

    f_1 = (f_11 * xdiff1 + f_21 * xdiff2) / xdiff3
    f_2 = (f_12 * xdiff1 + f_22 * xdiff2) / xdiff3

    ydiff1 = y_nodes[y_digi] - y_clip
    ydiff2 = y_clip - y_nodes[y_digi - 1]
    ydiff3 = y_nodes[y_digi] - y_nodes[y_digi - 1]
    f_out = 10 ** ((ydiff1 * f_1 + ydiff2 * f_2) / ydiff3)

    return f_out
