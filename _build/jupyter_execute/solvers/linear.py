#!/usr/bin/env python
# coding: utf-8

# # Linear solvers

# In[1]:


import numpy as np
from scipy.sparse.linalg import spsolve,bicg,gmres
import scipy.sparse as sp
import pypardiso
import numba

N = 2000
A_sp = sp.rand(N, N, density=0.01, format='csr')
b_sp = np.random.rand(N)
A_de = A_sp.toarray()


# ## numpy

# In[2]:


get_ipython().run_line_magic('timeit', 'np.linalg.solve(A_de,b_sp)')
np.linalg.solve(A_de,b_sp)


# ## pypardiso

# In[3]:


x = pypardiso.spsolve(A_sp, b_sp)
x


# In[4]:


get_ipython().run_line_magic('pinfo', 'pypardiso.ps')


# In[ ]:





# In[5]:


import cffi
from ctypes.util import find_library
import ctypes

mkl_rt = find_library('mkl_rt.1')


ffi = cffi.FFI()
ffi.cdef("int pardiso (void *pt[64], int *maxfct, int *mnum, int *mtype, int *phase, int *n, double *a, int *ia, int *ja, int *perm, int *nrhs, int *iparm, double *dparm, double *b, double *x, int *error);")
mkl_pardiso = ffi.dlopen(mkl_rt)

a = sp.coo_matrix(A_de).data
ia = sp.coo_matrix(A_de).row
ja =sp.coo_matrix(A_de).col

pA = ffi.from_buffer(A_sp.data)
pia = ffi.from_buffer(A_sp.indices)
pja = ffi.from_buffer(A_sp.indptr)

pb = ffi.from_buffer(b_sp)
px = ffi.from_buffer(np.copy(b_sp))

pt_type = (ctypes.c_int64, np.int64)

pt = ffi.new("void *[64]")
pt = ctypes.POINTER(pt_type[0])
maxfct = ffi.new("int *", 1)
mnum = ffi.new("int *", 1)
mtype = ffi.new("int *", 11)
phase = ffi.new("int *", 13)
n = ffi.new("int *", A_sp.shape[0])
perm = ffi.new("int []", A_sp.shape[0])
nrhs = ffi.new("int *", 1)
iparm = ffi.new("int []", 64)
dparm = ffi.new("double []", 64)
#b = ffi.new("double []", b_sp)
#x = ffi.new("double []", np.copy(b_sp))
error = ffi.new("int *", 0)

mkl_pardiso.pardiso(pt, maxfct, mnum, mtype, phase, n, pA, pia, pja, perm, nrhs, iparm, dparm, pb, px, error)


# In[9]:


import numba

from ctypes import c_int, c_void_p, c_double

pardiso_sig = numba.cfunc("c_int, c_int, c_int, c_int, c_int, c_void_p, c_int, c_void_p, c_int, c_int, c_int, c_int, c_int, c_void_p, c_void_p, c_void_p")


# In[19]:


get_ipython().run_line_magic('timeit', 'pypardiso.spsolve(A_sp, b_sp)')
#pypardiso.spsolve(A_sp, b_sp)


# In[10]:


from numba import cfunc, types, carray

c_sig = types.void(types.CPointer(types.double),
                   types.CPointer(types.double),
                   types.intc, types.intc)

@cfunc(c_sig)
def my_callback(in_, out, m, n):
    in_array = carray(in_, (m, n))
    out_array = carray(out, (m, n))
    for i in range(m):
        for j in range(n):
            out_array[i, j] = 2 * in_array[i, j]


# In[ ]:


c_sig = types.void(types.CPointer(types.double),
                   types.CPointer(types.double),
                   types.intc, types.intc)
                   
# Initialize the MKL Pardiso solver
iparm = ffi.new("int[64]")
pt = ffi.new("int[16]")
mtype = ffi.new("int *", 11)  # Matrix type: real symmetric indefinite
#mkl_pardiso.pardisoinit(pt, mtype, iparm)

N = 5
A_sp = sp.rand(N, N, density=1, format='csr')
b_sp = np.random.rand(N)



# Set the parameters for the MKL Pardiso solver
maxfct = ffi.new("int *", 1)
mnum = ffi.new("int *", 1)
phase = ffi.new("int *", 13)  # Analysis and factorization
n = ffi.new("int *", A_sp.shape[0])
a = ffi.new(f"double[{len(A_sp.data)}]", list(A_sp.data))
ia = ffi.new(f"int[{len(A_sp.indptr)}]", list(A_sp.indptr+1))
ja = ffi.new(f"int[{len(A_sp.indices)}]", list(A_sp.indices+1))
perm = ffi.new(f"int[{len(b_sp)}]")
nrhs = ffi.new("int *", 1)
msglvl = ffi.new("int *", 0)
b = ffi.new(f"double[{len(b_sp)}]", list(b_sp))
x = ffi.new(f"double[{len(b_sp)}]")
error = ffi.new("int *")


# In[5]:


@numba.njit()
def pardiso(A_sp,b_sp):
    pypardiso.spsolve(A_sp, b_sp)

pardiso(A_sp,b_sp)


# In[ ]:


gmres(A_sp, b_sp)


# In[ ]:





# In[2]:


import mkl


# In[9]:


from ctypes.util import find_library
import ctypes


# ### Working with ctypes

# In[1]:


mkl_rt = find_library('mkl_rt.1')
libmkl = ctypes.CDLL(mkl_rt)

mkl_pardiso = libmkl.pardiso

pt_type = (ctypes.c_int64, np.int64)

mkl_pardiso.argtypes = [ctypes.POINTER(pt_type[0]),    # pt
                                      ctypes.POINTER(ctypes.c_int32),      # maxfct
                                      ctypes.POINTER(ctypes.c_int32),      # mnum
                                      ctypes.POINTER(ctypes.c_int32),      # mtype
                                      ctypes.POINTER(ctypes.c_int32),      # phase
                                      ctypes.POINTER(ctypes.c_int32),      # n
                                      ctypes.POINTER(None),                # a
                                      ctypes.POINTER(ctypes.c_int32),      # ia
                                      ctypes.POINTER(ctypes.c_int32),      # ja
                                      ctypes.POINTER(ctypes.c_int32),      # perm
                                      ctypes.POINTER(ctypes.c_int32),      # nrhs
                                      ctypes.POINTER(ctypes.c_int32),      # iparm
                                      ctypes.POINTER(ctypes.c_int32),      # msglvl
                                      ctypes.POINTER(None),                # b
                                      ctypes.POINTER(None),                # x
                                      ctypes.POINTER(ctypes.c_int32)]      # error

mtype=11
phase=13
size_limit_storage=5e7

mkl_pardiso.restype = None

pt = np.zeros(64, dtype=pt_type[1])
iparm = np.zeros(64, dtype=np.int32)
perm = np.zeros(0, dtype=np.int32)

mtype = mtype
phase = phase
msglvl = False

factorized_A = sp.csr_matrix((0, 0))
size_limit_storage = size_limit_storage
_solve_transposed = False

x = np.zeros_like(b_sp)
pardiso_error = ctypes.c_int32(0)
c_int32_p = ctypes.POINTER(ctypes.c_int32)
c_float64_p = ctypes.POINTER(ctypes.c_double)
A = A_sp
b = b_sp

ia = A.indptr + 1
ja = A.indices + 1

mkl_pardiso(pt.ctypes.data_as(ctypes.POINTER(pt_type[0])),  # pt
                          ctypes.byref(ctypes.c_int32(1)),  # maxfct
                          ctypes.byref(ctypes.c_int32(1)),  # mnum
                          ctypes.byref(ctypes.c_int32(mtype)),  # mtype -> 11 for real-nonsymetric
                          ctypes.byref(ctypes.c_int32(phase)),  # phase -> 13
                          ctypes.byref(ctypes.c_int32(A.shape[0])),  # N -> number of equations/size of matrix
                          A.data.ctypes.data_as(c_float64_p),  # A -> non-zero entries in matrix
                          ia.ctypes.data_as(c_int32_p),  # ia -> csr-indptr
                          ja.ctypes.data_as(c_int32_p),  # ja -> csr-indices
                          perm.ctypes.data_as(c_int32_p),  # perm -> empty
                          ctypes.byref(ctypes.c_int32(1 if b.ndim == 1 else b.shape[1])),  # nrhs
                          iparm.ctypes.data_as(c_int32_p),  # iparm-array
                          ctypes.byref(ctypes.c_int32(msglvl)),  # msg-level -> 1: statistical info is printed
                          b_sp.ctypes.data_as(c_float64_p),  # b -> right-hand side vector/matrix
                          x.ctypes.data_as(c_float64_p),  # x -> output
                          ctypes.byref(pardiso_error))  # pardiso error

x


# In[17]:


import mkl
mkl.set_num_stripes


# In[11]:


mkl_rt = find_library('mkl_rt.1')
libmkl = ctypes.CDLL(mkl_rt)

mkl_pardiso = libmkl.pardiso

pt_type = (ctypes.c_int64, np.int64)

maxfct = np.array([1],dtype=np.int32)
p_maxfct = ffi.from_buffer(maxfct)

mkl_pardiso.argtypes = [ctypes.POINTER(pt_type[0]),    # pt
                                      ctypes.POINTER(ctypes.c_int32),      # maxfct
                                      ctypes.POINTER(ctypes.c_int32),      # mnum
                                      ctypes.POINTER(ctypes.c_int32),      # mtype
                                      ctypes.POINTER(ctypes.c_int32),      # phase
                                      ctypes.POINTER(ctypes.c_int32),      # n
                                      ctypes.POINTER(None),                # a
                                      ctypes.POINTER(ctypes.c_int32),      # ia
                                      ctypes.POINTER(ctypes.c_int32),      # ja
                                      ctypes.POINTER(ctypes.c_int32),      # perm
                                      ctypes.POINTER(ctypes.c_int32),      # nrhs
                                      ctypes.POINTER(ctypes.c_int32),      # iparm
                                      ctypes.POINTER(ctypes.c_int32),      # msglvl
                                      ctypes.POINTER(None),                # b
                                      ctypes.POINTER(None),                # x
                                      ctypes.POINTER(ctypes.c_int32)]      # error

mtype=11
phase=13
size_limit_storage=5e7

mkl_pardiso.restype = None

pt = np.zeros(64, dtype=pt_type[1])
iparm = np.zeros(64, dtype=np.int32)
perm = np.zeros(0, dtype=np.int32)

mtype = mtype
phase = phase
msglvl = False

factorized_A = sp.csr_matrix((0, 0))
size_limit_storage = size_limit_storage
_solve_transposed = False

x = np.zeros_like(b_sp)
pardiso_error = ctypes.c_int32(0)
c_int32_p = ctypes.POINTER(ctypes.c_int32)
c_float64_p = ctypes.POINTER(ctypes.c_double)
A = A_sp
b = b_sp

ia = A.indptr + 1
ja = A.indices + 1

def eval_mkl(A_data,A_shape):
    #pt_type = (ctypes.c_int64, np.int64)
    pardiso_error = ctypes.c_void_p(1)
    # c_int32_p = ctypes.POINTER(ctypes.c_int32)
    # c_float64_p = ctypes.POINTER(ctypes.c_double)
    # mkl_pardiso(pt.ctypes.data_as(ctypes.POINTER(ctypes.c_int64)),  # pt
    #                         ctypes.byref(ctypes.c_int32(1)),  # maxfct
    #                         ctypes.byref(ctypes.c_int32(1)),  # mnum
    #                         ctypes.byref(ctypes.c_int32(11)),  # mtype -> 11 for real-nonsymetric
    #                         ctypes.byref(ctypes.c_int32(13)),  # phase -> 13
    #                         ctypes.byref(ctypes.c_int32(A_shape)),  # N -> number of equations/size of matrix
    #                         A_data.ctypes.data_as(c_float64_p),  # A -> non-zero entries in matrix
    #                         ia.ctypes.data_as(c_int32_p),  # ia -> csr-indptr
    #                         ja.ctypes.data_as(c_int32_p),  # ja -> csr-indices
    #                         perm.ctypes.data_as(c_int32_p),  # perm -> empty
    #                         ctypes.byref(ctypes.c_int32(1)),  # nrhs
    #                         iparm.ctypes.data_as(c_int32_p),  # iparm-array
    #                         ctypes.byref(ctypes.c_int32(1)),  # msg-level -> 1: statistical info is printed
    #                         b_sp.ctypes.data_as(c_float64_p),  # b -> right-hand side vector/matrix
    #                         x.ctypes.data_as(c_float64_p),  # x -> output
    #                         ctypes.byref(pardiso_error))  # pardiso error

    # return x

maxfct = ctypes.byref(ctypes.c_int32(1))
mnum = ctypes.byref(ctypes.c_int32(1)),  # mnum
mtype = ctypes.byref(ctypes.c_int32(mtype)),  # mtype -> 11 for real-nonsymetric
phase = ctypes.byref(ctypes.c_int32(phase)),  # phase -> 13
N = ctypes.byref(ctypes.c_int32(A.shape[0])),  # N -> number of equations/size of matrix

A_data = A_sp.data
A_shape = A_sp.shape[0]
eval_mkl(A_data,A_shape)


# In[58]:


numba.int32


# In[1]:


from cffi import FFI
import cffi
from ctypes.util import find_library
import ctypes

mkl_rt = find_library('mkl_rt.1.dll')

ffi = FFI()

# Define the MKL Pardiso functions and variables
ffi.cdef("""
    void pardisoinit (void   *pt, int    *mtype, int *iparm);
    void pardiso     (void   *pt, int    *maxfct, int *mnum, int *mtype, int *phase,
                      int    *n, double *a, int *ia, int *ja, int    *perm, int *nrhs,
                      int    *iparm, int *msglvl, double *b, double *x, int *error);
""")

# Load the MKL Pardiso library
mkl_pardiso = ffi.dlopen(mkl_rt)

# Initialize the MKL Pardiso solver
iparm = ffi.new("int[64]")
pt = ffi.new("int *",10)
mtype = ffi.new("int *", -2)  # Matrix type: real symmetric indefinite
#mkl_pardiso.pardisoinit(pt, mtype, iparm)

# Set the parameters for the MKL Pardiso solver
maxfct = ffi.new("int *", 1)
mnum = ffi.new("int *", 1)
phase = ffi.new("int *", 13)  # Analysis and factorization
n = ffi.new("int *", 3)
a = ffi.new("double[9]", [4.0, 1.0, 2.0, 1.0, 3.0, 1.0, 2.0, 1.0, 4.0])
ia = ffi.new("int[4]", [1, 3, 4, 4])
ja = ffi.new("int[9]", [1, 2, 3, 1, 2, 3, 1, 2, 3])
perm = ffi.new("int[3]")
nrhs = ffi.new("int *", 1)
msglvl = ffi.new("int *", 0)
b = ffi.new("double[3]", [5.0, 7.0, 5.0])
x = ffi.new("double[3]")
error = ffi.new("int *")

# Solve the linear system using MKL Pardiso
mkl_pardiso.pardiso(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, perm, nrhs, iparm, msglvl, b, x, error)
print(x)  # Output: (1.0, 2.0, 2.0)

# Print the solution
print('hola')  # Output: (1.0, 2.0, 2.0)


# In[16]:


@numba.njit()
def call_pardiso(A_data,A_indices,A_indptr,A_shape):

    pardiso_error = 0


    mkl_pardiso(0,  # pt
                           1,  # maxfct
                           1,  # mnum
                           11,  # mtype -> 11 for real-nonsymetric
                           13,  # phase -> 13
                           A_shape.ctypes.data_as(ctypes.POINTER(ctypes.c_uc_longlong)),  # N -> number of equations/size of matrix
                           A_data,  # A -> non-zero entries in matrix
                           A_indptr,  # ia -> csr-indptr
                           A_indices,  # ja -> csr-indices
                           perm,  # perm -> empty
                           1,  # nrhs
                           iparm,  # iparm-array
                           1,  # msg-level -> 1: statistical info is printed
                           b_sp,  # b -> right-hand side vector/matrix
                           x,  # x -> output
                           0)  # pardiso error


c_int32_p = ctypes.POINTER(ctypes.c_int32)
c_float64_p = ctypes.POINTER(ctypes.c_double)

call_pardiso(A_sp.data,A_sp.indices,A_sp.indptr,A_sp.shape[0])


# In[59]:


pardiso_error


# In[18]:


get_ipython().run_line_magic('pinfo', 'mkl_pardiso.pardiso')


# In[19]:


mkl_pardiso.pardiso(ffi.NULL, iparm, ffi.NULL, ffi.NULL, ffi.NULL, ffi.NULL, a, ia, ja, ffi.NULL, ffi.NULL, ffi.NULL, pt, b, ffi.NULL, ffi.NULL)


# In[9]:


n = 100
iparm = ffi.new("int[64]")
pt = ffi.new("int[64]")
a = ffi.new("double[100]")
ia = ffi.new("int[100]")
ja = ffi.new("int[100]")
b = ffi.new("double[100]")

# Set up the input parameters
iparm[0] = 1 # Use default values for all parameters
iparm[1] = 2 # Use the 2-level factorization algorithm

# Set up the matrix A
for i in range(n):
  a[i] = A[i]
  ia[i] = IA[i]
  ja[i] = JA[i]



# Set up the right-hand side vector b
for i in range(n):
  b[i] = B[i]


# In[13]:


from numba import njit
import numpy as np
from numba.extending import get_cython_function_address
import ctypes

# Datatype pointers to give to the cython LAPACK functions
_PTR = ctypes.POINTER

_dbl = ctypes.c_double
_int = ctypes.c_int64

_ptr_dbl = _PTR(_dbl)
_ptr_int = _PTR(_int)

# zgges is the complex QZ-decomposition
zgges_addr = get_cython_function_address('scipy.linalg.cython_lapack', 'zgges')
zgges_functype = ctypes.CFUNCTYPE(None,
                                  _ptr_int,  # JOBVSL
                                  _ptr_int,  # JOBVSR
                                  _ptr_int,  # SORT
                                  _ptr_int,  # SELCTG
                                  _ptr_int,  # N
                                  _ptr_dbl,  # A, complex
                                  _ptr_int,  # LDA
                                  _ptr_dbl,  # B, complex
                                  _ptr_int,  # LDB
                                  _ptr_int,  # SDIM
                                  _ptr_dbl,  # ALPHA, complex
                                  _ptr_dbl,  # BETA, complex
                                  _ptr_dbl,  # VSL, complex
                                  _ptr_int,  # LDVSL
                                  _ptr_dbl,  # VSR, complex
                                  _ptr_int,  # LDVSR
                                  _ptr_dbl,  # WORK, complex
                                  _ptr_int,  # LWORK
                                  _ptr_dbl,  # RWORK
                                  _ptr_int,  # BWORK
                                  _ptr_int)  # INFO
zgges_fn = zgges_functype(zgges_addr)

@njit
def numba_zgges(x, y):
    _M, _N = x.shape

    A = x
    B = y

    JOBVSL = np.array([ord('V')], dtype=np.int64)
    JOBVSR = np.array([ord('V')], dtype=np.int64)
    SORT = np.array([ord('N')], dtype=np.int64)
    SELCTG = np.empty(1, dtype=np.int64)

    N = np.array(_N, dtype=np.int64)
    LDA = np.array(_N, dtype=np.int64)
    LDB = np.array(_N, dtype=np.int64)
    SDIM = np.array(0, dtype=np.int64) # out

    ALPHA = np.empty(_N, dtype=np.complex128) # out
    BETA = np.empty(_N, dtype=np.complex128) # out
    LDVSL = np.array(_N, dtype=np.int64)
    VSL = np.empty((_N, _N), dtype=np.complex128) # out
    LDVSR = np.array(_N, dtype=np.int64)
    VSR = np.empty((_N, _N), dtype=np.complex128) # out

    WORK = np.empty((1,), dtype=np.complex128) #out
    LWORK = np.array(-1, dtype=np.int64)
    RWORK = np.empty(8*_N, dtype=np.float64)
    BWORK = np.empty(_N, dtype=np.int64)
    INFO = np.empty(1, dtype=np.int64)

    zgges_fn(JOBVSL.ctypes,
             JOBVSR.ctypes,
             SORT.ctypes,
             SELCTG.ctypes,
             N.ctypes,
             A.view(np.float64).ctypes,
             LDA.ctypes,
             B.view(np.float64).ctypes,
             LDB.ctypes,
             SDIM.ctypes,
             ALPHA.view(np.float64).ctypes,
             BETA.view(np.float64).ctypes,
             VSL.view(np.float64).ctypes,
             LDVSL.ctypes,
             VSR.view(np.float64).ctypes,
             LDVSR.ctypes,
             WORK.view(np.float64).ctypes,
             LWORK.ctypes,
             RWORK.ctypes,
             BWORK.ctypes,
             INFO.ctypes)

    print("Calculated workspace size as", WORK[0])
    WS_SIZE = np.int64(WORK[0].real)
    LWORK = np.array(WS_SIZE, np.int64)
    WORK = np.empty(WS_SIZE, dtype=np.complex128)
    zgges_fn(JOBVSL.ctypes,
             JOBVSR.ctypes,
             SORT.ctypes,
             SELCTG.ctypes,
             N.ctypes,
             A.view(np.float64).ctypes,
             LDA.ctypes,
             B.view(np.float64).ctypes,
             LDB.ctypes,
             SDIM.ctypes,
             ALPHA.view(np.float64).ctypes,
             BETA.view(np.float64).ctypes,
             VSL.view(np.float64).ctypes,
             LDVSL.ctypes,
             VSR.view(np.float64).ctypes,
             LDVSR.ctypes,
             WORK.view(np.float64).ctypes,
             LWORK.ctypes,
             RWORK.ctypes,
             BWORK.ctypes,
             INFO.ctypes)

    # The LAPACK function also returns SDIM, WORK, BWORK, but I don't need them here.
    return A, B, ALPHA, BETA, VSL.T, VSR.T, INFO

for _ in range(10):
    n = 10
    A = np.random.random([n, n])
    B = np.random.random([n, n])
    AA, BB, a, b, Q, Z, info = numba_zgges(np.asfortranarray(A, dtype='D'), 
                                           np.asfortranarray(B, dtype='D'))

    aa = Q @ AA @ Z.conj().T
    assert np.allclose(aa.real, A)
    assert np.allclose(aa.imag, 0)
    bb = Q @ BB @ Z.conj().T
    assert np.allclose(bb.real, B)
    assert np.allclose(bb.imag, 0)
    assert np.allclose(Q @ Q.conj().T, np.eye(n))
    assert np.allclose(Z @ Z.conj().T, np.eye(n))

    assert np.all(np.diag(BB) >= 0)


# In[50]:


n = 10
A = np.random.random([n, n])
B = np.random.random([n, n])
AA, BB, a, b, Q, Z, info = numba_zgges(np.asfortranarray(A, dtype='D'), 
                                        np.asfortranarray(B, dtype='D'))


# In[11]:


import numba as nb
import numpy as np
import ctypes
mkl_rt = find_library('splev')

lib = ctypes.cdll.LoadLibrary("splev.dll")

dble_p=ctypes.POINTER(ctypes.c_double)
int_p =ctypes.POINTER(ctypes.c_longlong)

SPLEV=lib.SPLEV
SPLEV.restype =  ctypes.c_void_p
SPLEV.argtypes = (dble_p,int_p,dble_p,int_p,dble_p,dble_p,int_p,int_p,int_p)

from numba import types
from numba.extending import intrinsic
from numba.core import cgutils

@intrinsic
def val_to_ptr(typingctx, data):
    def impl(context, builder, signature, args):
        ptr = cgutils.alloca_once_value(builder,args[0])
        return ptr
    sig = types.CPointer(nb.typeof(data).instance_type)(nb.typeof(data).instance_type)
    return sig, impl

@intrinsic
def ptr_to_val(typingctx, data):
    def impl(context, builder, signature, args):
        val = builder.load(args[0])
        return val
    sig = data.dtype(types.CPointer(data.dtype))
    return sig, impl

#with intrinsics, temporary arrays are allocated on stack
#faster but much more relevant for functions with very low runtime
@nb.njit()
def splev_wrapped(x, coeff,e):
    #There are just pointers passed to the fortran function.
    #The arrays have to be contiguous!
    t=np.ascontiguousarray(coeff[0])
    x=np.ascontiguousarray(x)
    
    c=coeff[1]
    k=coeff[2]
    
    y=np.empty(x.shape[0],dtype=np.float64)
    
    n_arr=val_to_ptr(nb.int64(t.shape[0]))
    k_arr=val_to_ptr(nb.int64(k))
    m_arr=val_to_ptr(nb.int64(x.shape[0]))
    e_arr=val_to_ptr(nb.int64(e))
    ier_arr=val_to_ptr(nb.int64(0))
    
    SPLEV(t.ctypes,n_arr,c.ctypes,k_arr,x.ctypes,
        y.ctypes,m_arr,e_arr,ier_arr)
    return y, ptr_to_val(ier_arr)

#without using intrinsics
@nb.njit()
def splev_wrapped_2(x, coeff,e):
    #There are just pointers passed to the fortran function.
    #The arrays have to be contiguous!
    t=np.ascontiguousarray(coeff[0])
    x=np.ascontiguousarray(x)
    
    c=coeff[1]
    k=coeff[2]
    y=np.empty(x.shape[0],dtype=np.float64)
    
    n_arr = np.empty(1,  dtype=np.int64)
    k_arr = np.empty(1,  dtype=np.int64)
    m_arr = np.empty(1,  dtype=np.int64)
    e_arr = np.empty(1,  dtype=np.int64)
    ier_arr = np.zeros(1,  dtype=np.int64)
    
    n_arr[0]=t.shape[0]
    k_arr[0]=k
    m_arr[0]=x.shape[0]
    e_arr[0]=e
    
    SPLEV(t.ctypes,n_arr.ctypes,c.ctypes,k_arr.ctypes,x.ctypes,
        y.ctypes,m_arr.ctypes,e_arr.ctypes,ier_arr.ctypes)
    return y, ier_arr[0]


# In[12]:


find_library('splev')


# In[ ]:




