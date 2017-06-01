/harmonic/dynamical_matrix.py
_set_py_dynamical_matrix():
vec = vecs[k][i][l]改为vec = -vecs[k][j][l]
def set_dynamical_matrix(self, q, verbose=False):
        return self._set_py_dynamical_matrix(q, verbose=verbose)
        try:
            import phonopy._phonopy as phonoc
            self._set_c_dynamical_matrix(q)
        except ImportError:
            self._set_py_dynamical_matrix(q, verbose=verbose)
/harmonic/derivative_dynmat.py
vecs_multi = vecs[k, i, :multi]改为vecs_multi = -vecs[k, j, :multi]
def run(self, q, q_direction=None, lang='P'):
    if self._derivative_order is not None or lang != 'C':
        self._run_py(q, q_direction=q_direction)
    else:
        self._run_c(q, q_direction=q_direction)
/phonon/dos.py
i_x = np.arange(num_atom, dtype='int') * 3
i_y = np.arange(num_atom, dtype='int') * 3 + 1
i_z = np.arange(num_atom, dtype='int') * 3 + 2
filter=np.abs(self._frequencies)<1e-6
for i,x in enumerate(filter):
  self._eigenvectors[i,:,x]=0.0