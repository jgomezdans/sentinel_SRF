import glob

import h5py
import numpy as np
from netCDF4 import Dataset


def create_py6s_representation():
    for ii, bp in enumerate(modis_band_names):
        band_name = "ACCURATE_MODIS_%s_%d_%s" % (sensor, ii + 1, bp)
        splodge = "".join(["%0.5f, " % b for b in band_pass[bp]["rsr"]])
        code_block = "%s = ( %d, %0.5f, %0.5f, \n" % (
            band_name,
            ii + band_step,
            band_pass[bp]["min"],
            band_pass[bp]["max"],
        )
        code_block += "\t\t np.array([%s]))\n\n" % splodge
        print(code_block)


def olci_srf(step, fname="original_data/OLCI_SRF.nc"):
    f = Dataset(fname, "r")
    y = f["spectral_response_function"][:, :, :, :].mean(axis=(1, 2))
    x = f["spectral_response_function_wavelength"][:, :, :, :].mean(axis=(1, 2))
    n_bands = x.shape[0]  # should be 21...
    wv = np.arange(250, 4002.5, 2.5)  # 6s sampling
    for band in range(n_bands):
        y_interp = np.interp(wv, x[band, :], y[band, :], left=0.0, right=0)
        passer = y_interp >= 0.005
        code_block = "%s = ( %d, %0.5f, %0.5f, \n" % (
            "S3A_OLCI_%02d" % (band + 1),
            band + step,
            wv[passer][0] / 1000.0,
            wv[passer][-1] / 1000.0,
        )
        splodge = "".join(["%0.5f, " % r for r in y_interp[passer]])
        code_block += "\t\t np.array([%s]))\n\n" % splodge
        print(code_block)
    return band + step


def slstr_srf(step):

    files = glob.glob("data/SLSTR_FM02_S*")
    files.sort()
    wv = np.arange(250, 4002.5, 2.5)  # 6s sampling
    for band, fich in enumerate(files):

        f = Dataset(fich, "r")
        x = 1000.0 * f["wavelength"][:]
        y = f["response"][:]

        y_interp = np.interp(wv, x, y, left=0.0, right=0)
        passer = y_interp >= 0.005
        code_block = "%s = ( %d, %0.5f, %0.5f, \n" % (
            "S3A_SLSTR_%02d" % (band + 1),
            band + step,
            wv[passer][0] / 1000.0,
            wv[passer][-1] / 1000.0,
        )
        splodge = "".join(["%0.5f, " % r for r in y_interp[passer]])
        code_block += "\t\t np.array([%s]))\n\n" % splodge
        print(code_block)
    return band + step


def s2a_srf(step):

    d = np.loadtxt(
        "data/s2a_spectral_response_functions.csv", skiprows=1, delimiter=","
    )
    wv = np.arange(250, 4002.5, 2.5)  # 6s sampling
    x = d[:, 0]
    for band in range(1, 14):
        y = d[:, band]
        y_interp = np.interp(wv, x, y, left=0.0, right=0)
        # passer = y_interp >= 0.005
        passer = y_interp >= 0.001
        code_block = "%s = ( %d, %0.5f, %0.5f, \n" % (
            "S2A_MSI_%02d" % (band),
            band + step,
            wv[passer][0] / 1000.0,
            wv[passer][-1] / 1000.0,
        )
        splodge = "".join(["%0.5f, " % r for r in y_interp[passer]])
        code_block += "\t\t np.array([%s]))\n\n" % splodge
        print(code_block)
    return band + step


def himawari8_ahi(step):
    wv = np.arange(250, 4002.5, 2.5)  # 6s sampling

    with h5py.File("data/himawari_8_ahi_RSR.h5") as f:
        x = f["wavelength"][()]
        rsr = f["RSR"][()]
    for band in range(6):
        y = rsr[band, :]
        y_interp = np.interp(wv, x, y, left=0.0, right=0)
        # passer = y_interp >= 0.005
        passer = y_interp >= 0.001
        code_block = "%s = ( %d, %0.5f, %0.5f, \n" % (
            "H8_AHI_%02d" % (band + 1),
            band + step,
            wv[passer][0] / 1000.0,
            wv[passer][-1] / 1000.0,
        )
        splodge = "".join(["%0.5f, " % r for r in y_interp[passer]])
        code_block += "\t\t np.array([%s]))\n\n" % splodge
        print(code_block)
    return band + step


if __name__ == "__main__":
    step = s2a_srf(52)
    step = olci_srf(step + 1)
    step = slstr_srf(step + 1)
    step = 143
    step = himawari8_ahi(step + 1)
