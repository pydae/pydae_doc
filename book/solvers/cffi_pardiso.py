import cffi
from ctypes.util import find_library

mkl_rt = find_library('mkl_rt.1')

ffi = cffi.FFI()
mkl = ffi.dlopen(mkl_rt)