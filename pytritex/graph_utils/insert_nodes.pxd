import numpy as np
cimport numpy as np
ctypedef np.int_t DTYPE_i
ctypedef np.float_t DTYPE_f


cpdef np.ndarray[DTYPE_i, ndim=1] insert_nodes(np.ndarray[DTYPE_f, ndim=2, mode="c"] edges,
                                               np.ndarray[DTYPE_i, ndim=1, mode="c"] path)