# JM: 10 Jun 2020
# checks to compare numbers from some checkpoint output with
# output after changes, for development of GEOMETRIC MPI code

import matplotlib
matplotlib.use("Agg")
matplotlib.rc("text", usetex = True)

from MITgcmutils import mnc
import matplotlib.pyplot as plt

# load the checkpoint calculation
nc = mnc.MNC("../../../mpi_testing/run/mnc_*/GEOMdiag2d.*")
for var in nc.variables:
  print(var)
ene = nc.variables["GEOMeE"][:]
Egen = nc.variables["GEOMEgen"][:]
Eadv = nc.variables["GEOMEadv"][:]
xT  = nc.variables["X"][:]
yT  = nc.variables["Y"][:]
nc.close()

# load the new calculation for comparison
nc = mnc.MNC("../mnc_*/GEOMdiag2d.*")
ene_new = nc.variables["GEOMeE"][:]
Egen_new = nc.variables["GEOMEgen"][:]
Eadv_new = nc.variables["GEOMEadv"][:]
xT  = nc.variables["X"][:]
yT  = nc.variables["Y"][:]
nc.close()

index = -1

fig = plt.figure(figsize = (10, 4))
plt.pcolormesh(xT / 1e3, yT / 1e3, (ene_new - ene)[index, 0, :, :], cmap = "RdBu_r")
plt.colorbar()
fig.savefig("./raw_errors_in_ene.png", dpi = 150, bbox_inches = "tight")
plt.close(fig)

fig = plt.figure(figsize = (10, 4))
plt.pcolormesh(xT / 1e3, yT / 1e3, (Egen_new - Egen)[index, 0, :, :], cmap = "RdBu_r")
plt.colorbar()
fig.savefig("./raw_errors_in_Egen.png", dpi = 150, bbox_inches = "tight")
plt.close(fig)

fig = plt.figure(figsize = (10, 4))
plt.pcolormesh(xT / 1e3, yT / 1e3, (Eadv_new - Eadv)[index, 0, :, :], cmap = "RdBu_r")
plt.colorbar()
fig.savefig("./raw_errors_in_Eadv.png", dpi = 150, bbox_inches = "tight")
plt.close(fig)
