import numpy as np



def read_physical_quantities(path):
    """
    Read a 'physical_quantities.txt'-style file and return
    density, mean_velocity, temperature as NumPy 1D arrays.

    Args:
        path (str): The file path to the 'physical_quantities.txt' file.

    Returns:
        tuple: A tuple containing three NumPy 1D arrays (density, mean_velocity, temperature).

    """
    # Determine how many header rows to skip (until the first numeric index line)
    skip = 0
    with open(path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                skip += 1
                continue
            parts = stripped.split()
            try:
                int(parts[0])  # first token should be the integer index
                break
            except (ValueError, IndexError):
                skip += 1

    data = np.loadtxt(path, usecols=(1, 2, 3), skiprows=skip)
    if data.ndim == 1:
        # Single data row case
        data = data[None, :]

    density = data[:, 0]
    mean_velocity = data[:, 1]
    temperature = data[:, 2]
    return density, mean_velocity, temperature

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def read_physical_quantity(path):
    """
    Read a 'physical_quantities.txt'-style file and return only density as a NumPy 1D array.
    """
    # Determine how many header rows to skip (until the first numeric index line)
    skip = 0
    with open(path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                skip += 1
                continue
            parts = stripped.split()
            try:
                int(parts[0])  # first token should be the integer index
                break
            except (ValueError, IndexError):
                skip += 1

    quantity = np.loadtxt(path, usecols=(1,), skiprows=skip)
    return quantity

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def read_space_mesh(path):
    """
    Read a 'space_mesh.txt'-style file ('Computational Points (index - x_comp):')
    and return only x_comp as a NumPy 1D array.
    """
    # Find start of numeric data (first line with int index and float x)
    skip = 0
    with open(path, 'r') as f:
        for line in f:
            s = line.strip()
            if not s:
                skip += 1
                continue
            parts = s.split()
            try:
                int(parts[0])
                float(parts[1])
                break
            except (ValueError, IndexError):
                skip += 1

    x_comp = np.loadtxt(path, usecols=(1,), skiprows=skip)
    return x_comp

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

