# cython: boundscheck=False, wraparound=False

def strongest_path(long size, double [:] pref, double [:] spath):

    cdef Py_ssize_t i, j, k
    cdef double ji, path

    for i in range(size):
        for j in range(size):
            if i != j:
                if pref[i*size + j] > pref[j*size + i]:
                    spath[i*size + j] = pref[i*size + j]

    for i in range(size):
        for j in range(size):
            if i != j:
                # spath[j*size + i] does not change within the k loop (k != i)
                ji = spath[j*size + i]
                for k in range(size):
                    if (i != k) and (j != k):
                        # spath[j, k] = max(spath[j, k], min(spath[j, i], spath[i, k]))
                        path = spath[i*size + k]
                        if ji < path:
                            path = ji
                        if path > spath[j*size + k]:
                            spath[j*size + k] = path
