import numpy as np
from scipy.integrate import solve_ivp
from typing import Tuple, List, Callable, Any
from .physics import Spacetime
import logging

logging.basicConfig(level=logging.INFO)

class HamiltonianIntegrator:
    """
    Pembungkus untuk SciPy ODE Solver guna mengintegrasikan sistem Hamiltonian.
    Menggunakan metode RK45 yang sangat efisien dengan adaptive step-size.
    """
    def __init__(self, spacetime: Spacetime):
        self.spacetime = spacetime

    def hamiltonian_system(self, tau: float, state: np.ndarray) -> np.ndarray:
        """
        Fungsi ODE untuk solve_ivp. 
        Menerima array `state` 8D: [t, r, theta, phi, pt, pr, ptheta, pphi].
        Berdasarkan Persamaan Hamilton:
            dx/dtau =  dH/dp
            dp/dtau = -dH/dx
        """
        x = state[:4]
        p = state[4:]
        
        dx_dtau = self.spacetime.dHdp(x, p)
        dp_dtau = -self.spacetime.dHdx(x, p)
        
        return np.concatenate((dx_dtau, dp_dtau))

    def solve(self, 
              x0: np.ndarray, 
              p0: np.ndarray, 
              tau_span: Tuple[float, float], 
              events: List[Callable] = None,
              max_step: float = 0.5) -> Any:
        """
        Menjalankan integrasi ODE adaptif.
        """
        initial_state = np.concatenate((x0, p0))
        
        # Menggunakan metode 'RK45' sebagai standar terbaik untuk masalah ODE non-stiff
        solution = solve_ivp(
            fun=self.hamiltonian_system,
            t_span=tau_span,
            y0=initial_state,
            method='RK45',
            events=events,
            max_step=max_step,
            vectorized=False,
            dense_output=True, # Menyediakan interpolasi kurva yang halus
            rtol=1e-8,
            atol=1e-10
        )
        return solution

