from scipy.integrate import solve_ivp
from openpnm.algorithms._solution import TransientSolution
from openpnm.integrators import Integrator

class ScipyBDF(Integrator):
    """Custom integrator based on SciPy's BDF method (for stiff systems)"""

    def __init__(self, atol=1e-6, rtol=1e-6, verbose=False, linsolver=None):
        self.atol = atol
        self.rtol = rtol
        self.verbose = verbose
        self.linsolver = linsolver

    def solve(self, rhs, x0, tspan, saveat, **kwargs):
        options = {
            "atol": self.atol,
            "rtol": self.rtol,
            "t_eval": saveat,
        }
        sol = solve_ivp(rhs, tspan, x0, method="BDF", **options)
        if sol.success:
            return TransientSolution(sol.t, sol.y)
        raise Exception(sol.message)
