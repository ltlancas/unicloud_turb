## generate.py
## Code to generate a Gaussian Random Field for velocities as a standard model for cloud
## scale star formation simulations. The field is guaranteed to be solenoidal.

import numpy as np
import h5py


def init_velocity_field_proj(sigma, kspec, kmin, kmax, N, seed=42):
    """
    Initializes a 3D, periodic gaussian random field with a given power spectrum

    Parameters
    ----------
    kspec : float
        The spectral index of the power spectrum
    sigma : float
        The standard deviation of the velocity field
    kmin : float
        The minimum wavenumber for the power spectrum
    kmax : float
        The maximum wavenumber for the power spectrum
    N : int
        The size of the grid (N x N x N)
    seed: int, optional
        Random seed (default 42)
    """
    kmin *= 2 * np.pi
    kmax *= 2 * np.pi
    k = np.fft.fftfreq(N, d=1.0 / N) * 2 * np.pi
    KS = np.meshgrid(k, k, k, indexing="ij")
    KS = np.array(KS)
    (KX, KY, KZ) = KS
    K = np.sqrt(KX**2 + KY**2 + KZ**2)
    mask = K > 0
    mask1 = K < kmin
    mask2 = K > kmax

    VK = []
    for i in range(3):
        np.random.seed(seed + i)
        vk = np.zeros((N, N, N), dtype=complex)
        phase = np.fft.fftn(np.random.normal(size=(N, N, N)))
        vk[mask] = (K[mask] ** (-1 * kspec)) * phase[mask]
        vk[mask1] = 0
        vk[mask2] *= np.exp(-1 * (K[mask2] / kmax) ** 2) * np.e
        VK.append(vk)
    VK = np.array(VK)

    VKP = np.zeros_like(VK, dtype=complex)
    for i in range(3):
        for j in range(3):
            if i == j:
                VKP[i] += VK[j]
            VKP[i][mask] -= (KS[i] * KS[j] * VK[j])[mask] / (K[mask] ** 2)

    (vx, vy, vz) = [np.fft.ifftn(vk).real for vk in VKP]
    sigma_res = np.sqrt(np.std(vx) ** 2 + np.std(vy) ** 2 + np.std(vz) ** 2)

    vx *= sigma / sigma_res
    vy *= sigma / sigma_res
    vz *= sigma / sigma_res
    return (vx, vy, vz)


def make_hermitian(arr, N):
    # guarantees that the fourier transform is real
    # first make sure maximum and zero frequency phases are zero
    a = np.copy(arr)
    a[0,0,0] = 0
    # set nyquist frequncies to zero
    a[N//2,:,:] = 0
    a[:,N//2,:] = 0
    a[:,:,N//2] = 0
    # deal with zeros in two places
    for i in range(1,N//2):
        a[i,0,0] = np.conjugate(a[-i,0,0])
        a[0,i,0] = np.conjugate(a[0,-i,0])
        a[0,0,i] = np.conjugate(a[0,0,-i])
    # deal with zero in one place
    for i in range(1,N//2):
        for j in range(1,N//2):
            a[ i, j, 0] = np.conjugate(a[-i,-j, 0])
            a[-i, j, 0] = np.conjugate(a[ i,-j, 0])
            a[ i, 0, j] = np.conjugate(a[-i, 0,-j])
            a[-i, 0, j] = np.conjugate(a[ i, 0,-j])
            a[ 0, i, j] = np.conjugate(a[ 0,-i,-j])
            a[ 0,-i, j] = np.conjugate(a[ 0, i,-j])
    # deal with the bulk of the grid
    for i in range(1,N//2):
        for j in range(1,N//2):
            for k in range(1,N//2):
                a[ i, j, k] = np.conjugate(a[-i,-j,-k])
                a[-i, j, k] = np.conjugate(a[ i,-j,-k])
                a[ i,-j, k] = np.conjugate(a[-i, j,-k])
                a[ i, j,-k] = np.conjugate(a[-i,-j, k])
    return a


def init_vec_potential(kspec, kmin, kmax, N, seed=42):
    """
    Initializes a 3D, periodic gaussian random field with a given power spectrum. Used to
    initialize a vector potential that can then be used in real or fourier space to
    create a divergence free velocity field

    Parameters
    ----------
    kspec : float
        The spectral index of the power spectrum
    kmin : float
        The minimum wavenumber for the power spectrum
    kmax : float
        The maximum wavenumber for the power spectrum
    N : int
        The size of the grid (N x N x N)
    seed: int, optional
        Random seed (default 42)
    """
    kmin *= 2 * np.pi
    kmax *= 2 * np.pi
    k = np.fft.fftfreq(N, d=1.0 / N) * 2 * np.pi
    KS = np.meshgrid(k, k, k, indexing="ij")
    (KX, KY, KZ) = KS
    KS = np.array(KS)
    K = np.sqrt(KX**2 + KY**2 + KZ**2)
    mask = K > 0
    mask1 = K < kmin
    mask2 = K > kmax

    AK = []
    for i in range(3):
        np.random.seed(seed + i)
        ak = np.zeros((N, N, N), dtype=complex)
        phase = np.exp(1.j*2*np.pi*np.random.normal(size=(N, N, N)))
        sigma_norm = (K**2 + 1e-300**2) ** (-0.5 * (kspec+1))
        norms = np.random.rayleigh(sigma_norm)
        ak[mask] = norms[mask]* phase[mask]
        ak = make_hermitian(ak, N)
        ak[mask1] = 0
        ak[mask2] *= np.exp(-1 * (K[mask2] / kmax) ** 2) * np.e
        AK.append(ak)
    AK = np.array(AK)
    (ax, ay, az) = [np.fft.ifftn(ak).real for ak in AK]

    return (ax, ay, az)


def init_velocity_field_kcurl(sigma, kspec, kmin, kmax, N, seed=42):
    """
    Initializes a 3D, periodic gaussian random field with a given power spectrum

    Parameters
    ----------
    kspec : float
        The spectral index of the power spectrum
    sigma : float
        The standard deviation of the velocity field
    kmin : float
        The minimum wavenumber for the power spectrum
    kmax : float
        The maximum wavenumber for the power spectrum
    N : int
        The size of the grid (N x N x N)
    seed: int, optional
        Random seed (default 42)
    """
    k = np.fft.fftfreq(N, d=1.0 / N) * 2 * np.pi
    KS = np.meshgrid(k, k, k, indexing="ij")
    (KX, KY, KZ) = KS
    KS = np.array(KS)
    K = np.sqrt(KX**2 + KY**2 + KZ**2)

    (ax, ay, az) = init_vec_potential(kspec, kmin, kmax, N, seed=seed)
    A = np.array([ax, ay, az])
    AK = np.array([np.fft.fftn(a) for a in A])
    VK = np.zeros_like(AK, dtype=complex)

    VK[0] = 1.j*(KS[1]*AK[2] - KS[2]*AK[1])
    VK[1] = 1.j*(KS[2]*AK[0] - KS[0]*AK[2])
    VK[2] = 1.j*(KS[0]*AK[1] - KS[1]*AK[0])

    (vx, vy, vz) = [np.fft.ifftn(vk).real for vk in VK]
    sigma_res = np.sqrt(np.std(vx)** 2 + np.std(vy)** 2 + np.std(vz)** 2)

    DIVVK = np.sum(1.j*KS * VK, axis=0)
    divv = (sigma/sigma_res)*np.fft.ifftn(DIVVK).real

    vx *= sigma / sigma_res
    vy *= sigma / sigma_res
    vz *= sigma / sigma_res
    return (vx, vy, vz)


def init_velocity_field_curl(sigma, kspec, kmin, kmax, N, seed=42):
    """
    Initializes a 3D, periodic gaussian random field with a given power spectrum

    Parameters
    ----------
    kspec : float
        The spectral index of the power spectrum
    sigma : float
        The standard deviation of the velocity field
    kmin : float
        The minimum wavenumber for the power spectrum
    kmax : float
        The maximum wavenumber for the power spectrum
    N : int
        The size of the grid (N x N x N)
    seed: int, optional
        Random seed (default 42)
    """
    k = np.fft.fftfreq(N, d=1.0 / N) * 2 * np.pi
    KS = np.meshgrid(k, k, k, indexing="ij")
    (KX, KY, KZ) = KS
    KS = np.array(KS)
    K = np.sqrt(KX**2 + KY**2 + KZ**2)

    (ax, ay, az) = init_vec_potential(kspec, kmin, kmax, N, seed=seed)

    # take the curl of the vector potential
    vx  = np.roll(az, -1, axis=1) - np.roll(az, 1, axis=1)
    vx -= np.roll(ay, -1, axis=2) - np.roll(ay, 1, axis=2)
    vy  = np.roll(ax, -1, axis=2) - np.roll(ax, 1, axis=2)
    vy -= np.roll(az, -1, axis=0) - np.roll(az, 1, axis=0)
    vz  = np.roll(ay, -1, axis=0) - np.roll(ay, 1, axis=0)
    vz -= np.roll(ax, -1, axis=1) - np.roll(ax, 1, axis=1)
    sigma_res = np.sqrt(np.std(vx)** 2 + np.std(vy)** 2 + np.std(vz)** 2)

    divv  = np.roll(vx, -1, axis = 0) -  np.roll(vx, 1, axis = 0)
    divv += np.roll(vy, -1, axis = 1) -  np.roll(vy, 1, axis = 1)
    divv += np.roll(vz, -1, axis = 2) -  np.roll(vz, 1, axis = 2)

    vx *= sigma / sigma_res
    vy *= sigma / sigma_res
    vz *= sigma / sigma_res
    return (vx, vy, vz)


def cell_center_coordinates(N: int = 128) -> tuple:
    """
    Returns a tuple of 3 arrays containing the Cartesian coordinate functions at the cell centers in the unit cube.

    Parameters
    ----------
    N: Number of cells per dimension

    Returns
    -------
    x, y, z:
    tuple of shape (N,N,N) arrays storing the coordinate functions, indexed such that
    increasing index corresponds to increasing coordinate value, and indices are ordered
    [x,y,z]
    """
    coords = (0.5 + np.arange(N)) / N
    return np.meshgrid(coords, coords, coords, indexing="ij")


def generate_velocity_cube(
    seed: int = 0,
    sigma: float = 1.0,
    kspec: float = 2.0,
    kmin: int = 0,
    kmax: int = 16,
    N: int = 128,
):
    """Generates an hdf5 file containing the Cartesian coordinate functions
    and velocity components.

    Parameters
    ----------
    seed: int, optional
        Random seed (default 0)
    kspec : float, optional
        The spectral index of the power spectrum
    sigma : float, optional
        The standard deviation of the velocity field
    kmin : float, optional
        The minimum wavenumber for the power spectrum (default: )
    kmax : float, optional
        The maximum wavenumber for the power spectrum
    N : int, optional
        The size of the grid (N x N x N)
    """
    # Lorentz's Birthday: 18th of July, 1853 - adding this to the random seed so seed 0 is really this.
    lorentz_bday = 18071853

    (vx, vy, vz) = init_velocity_field_proj(sigma, kspec, kmin, kmax, N, seed=seed + lorentz_bday)
    x, y, z = cell_center_coordinates(N)

    # Use H5py to create a HDF5 file that stores the velocity field information
    f = h5py.File(f"velocity_field_seed{seed}.h5", "w")
    f.create_dataset("vx", data=vx)
    f.create_dataset("vy", data=vy)
    f.create_dataset("vz", data=vz)
    f.create_dataset("x", data=x)
    f.create_dataset("y", data=y)
    f.create_dataset("z", data=z)
    f.close()


def generate_all_seeds():
    for i in range(100):
        generate_velocity_cube(seed=i)


if __name__ == "__main__":
    generate_all_seeds()
