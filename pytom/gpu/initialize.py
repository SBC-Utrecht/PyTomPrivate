import os
global xp
global device
global map_coordinates


try:
    import cupy as xp
    import cupy.typing as xpt
    n_gpus = xp.cuda.runtime.getDeviceCount()
    all_gpus = list(range(n_gpus)) 
    ID = all_gpus #TODO: remove ID and move everything to all_gpus
    xp.cuda.Device().use()
    from cupyx.scipy.ndimage import map_coordinates
    device = f'gpu:{xp.cuda.runtime.getDevice()}'
except ImportError:
    print("Could not import cupy, using cpu instead")
    import numpy as xp
    import numpy.typing as xpt
    from scipy.ndimage import map_coordinates
    ID = []
    all_gpus = []
    device = 'cpu'
except Exception as e:
    print(f"found unexpected error:\n{e}\ndefaulting to cpu")
    import numpy as xp
    import numpy.typing as xpt
    from scipy.ndimage import map_coordinates
    ID = []
    all_gpus = []
    device = 'cpu'

def initialize_gpu(id):
    global xp, device
    import cupy as xp
    xp.cuda.Device(int(id)).use()
    from cupyx.scipy.ndimage import map_coordinates
    device = f'gpu:{id}'
