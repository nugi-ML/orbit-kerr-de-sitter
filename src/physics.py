import numpy as np
from abc import ABC, abstractmethod
import logging

try:
    from .kerr_equations import eval_H, eval_dHdx, eval_dHdp, eval_gUU
except ImportError:
    raise ImportError("File kerr_equations.py tidak ditemukan. Harap jalankan 'python src/generate_equations.py' terlebih dahulu!")

class Spacetime(ABC):
    """Abstract Base Class untuk mendefinisikan geometri ruang-waktu."""
    @abstractmethod
    def H(self, x: np.ndarray, p: np.ndarray) -> float:
        pass

    @abstractmethod
    def dHdx(self, x: np.ndarray, p: np.ndarray) -> np.ndarray:
        pass

    @abstractmethod
    def dHdp(self, x: np.ndarray, p: np.ndarray) -> np.ndarray:
        pass

class KerrDeSitter(Spacetime):
    """
    Geometri Kerr-De Sitter berkinerja tinggi.
    Menggunakan rumus pre-computed statis tanpa overhead komputasi SymPy.
    """
    def __init__(self, M: float, a: float, lambda_cosmo: float):
        self.M = M
        self.a = a
        self.Lambda = lambda_cosmo
        logging.info(f"Menginisialisasi Kerr-De Sitter statis: M={M}, a={a}, Lambda={lambda_cosmo}")

    def H(self, x: np.ndarray, p: np.ndarray) -> float:
        return eval_H(x[0], x[1], x[2], x[3], p[0], p[1], p[2], p[3], self.M, self.a, self.Lambda)

    def dHdx(self, x: np.ndarray, p: np.ndarray) -> np.ndarray:
        return eval_dHdx(x[0], x[1], x[2], x[3], p[0], p[1], p[2], p[3], self.M, self.a, self.Lambda)

    def dHdp(self, x: np.ndarray, p: np.ndarray) -> np.ndarray:
        return eval_dHdp(x[0], x[1], x[2], x[3], p[0], p[1], p[2], p[3], self.M, self.a, self.Lambda)

    def get_initial_pr(self, x: np.ndarray, E: float, L_ang: float) -> float:
        """Mencari nilai momentum radial awal (pr) untuk geodesik null (H=0)."""
        gUU_tt, gUU_rr, gUU_tp, gUU_pp = eval_gUU(x[0], x[1], x[2], x[3], self.M, self.a, self.Lambda)
        
        if abs(gUU_rr) < 1e-10:
            raise ValueError(f"Singularitas koordinat (g^rr sangat kecil): {gUU_rr} di r={x[1]}")
            
        residue = gUU_tt * E**2 + gUU_pp * L_ang**2 + 2 * gUU_tp * (-E) * L_ang
        pr_sq = -residue / gUU_rr
        
        if pr_sq < 0:
            raise ValueError(f"Area terlarang foton: pr^2 = {pr_sq} < 0 (E={E}, L={L_ang})")
            
        return -np.sqrt(pr_sq)

