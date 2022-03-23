import numpy as np
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
from swiftsimio.visualisation.projection_backends import backends, backends_parallel
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pl
import scipy.optimize as opt

def gauss_curve(x, A, a, b, c, x0, y0):
  return A * np.exp(
    -(a*(x[:,0]-x0)**2 + 2.*b*(x[:,0]-x0)*(x[:,1]-y0) + c * (x[:,1]-y0)**2)
  )

def calculate_HI_size(data, ang_momentum, galaxy_data, output_path, index):
  pos_parts = data[:,0:3].copy()

  face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum, axis="z")
  pos_face_on = np.matmul(face_on_rotation_matrix, pos_parts.T).T

  hsml = data[:, 7]
  mass = data[:, 8]

  x_min = pos_face_on[:,0].min()
  x_max = pos_face_on[:,0].max()
  y_min = pos_face_on[:,1].min()
  y_max = pos_face_on[:,1].max()

  x_range = x_max - x_min
  y_range = y_max - y_min

  x = (pos_face_on[:,0] - x_min) / x_range
  y = (pos_face_on[:,1] - y_min) / y_range

  image = backends["subsampled"](x=x, y=y, m=mass, h=hsml, res=128)
  imin = image.min()
  imax = image.max()

  xg, yg = np.meshgrid(np.linspace(x_min, x_max, 128), np.linspace(y_min, y_max, 128))
  xs = np.zeros((128*128, 2))
  xs[:,0] = xg.flatten()
  xs[:,1] = yg.flatten()
  A = imin
  sigX = 20 * x_range / 128
  sigY = 20 * y_range / 128
  a = 0.5 / sigX**2
  b = 0.
  c = 0.5 / sigY**2
  params, _ = opt.curve_fit(gauss_curve, xs, image.flatten(), p0 = (A, a, b, c, 0., 0.))
  print(params)

  fig, ax = pl.subplots(1, 2)
  ax[0].imshow(image, norm=matplotlib.colors.LogNorm(vmin=imin, vmax=imax))
  image = gauss_curve(xs, *params).reshape((128, 128))
  ax[1].imshow(image, norm=matplotlib.colors.LogNorm(vmin=imin, vmax=imax))
  pl.savefig(f"{output_path}/HI_size_{index}.png", dpi=300)
  pl.close()
