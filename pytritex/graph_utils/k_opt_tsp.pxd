import numpy as np
cimport numpy as np


cpdef float route_cost(np.ndarray graph, np.ndarray path)
cpdef np.ndarray tsp_2_opt(np.ndarray graph, np.ndarray route, long num_threads)
