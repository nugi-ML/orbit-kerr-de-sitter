import os
import sys
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
from dataclasses import dataclass
from typing import List, Dict

# Memastikan direktori src dikenali oleh python path jika dijalankan dari luar
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.physics import KerrDeSitter
from src.integrator import HamiltonianIntegrator

logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')

@dataclass
class SimulationConfig:
    M: float = 1.0
    a: float = 0.99
    lambdas: List[float] = None
    r0: float = 500.0
    theta0: float = np.pi / 2
    phi0: float = 0.0
    E: float = 1.0
    b: float = 7.0
    tau_max: float = 2000.0  # Waktu affine maksimum untuk solver
    
    def __post_init__(self):
        if self.lambdas is None:
            self.lambdas = [0.0, 1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6]
        self.L_ang = self.b * self.E
        self.rh = self.M + np.sqrt(self.M**2 - self.a**2)

def make_event_horizon(rh: float):
    """Event untuk menghentikan simulasi jika foton menyentuh horizon"""
    def event(tau, state):
        r = state[1]
        return r - (rh * 1.1)
    event.terminal = True
    event.direction = -1
    return event

def make_event_escape(r0: float):
    """Event untuk menghentikan simulasi jika foton memantul kembali (scattering)"""
    def event(tau, state):
        r = state[1]
        # Pastikan ini hanya berhenti jika ia sudah sempat mendekat dan keluar lagi
        # Menggunakan sedikit toleransi r0 + 5.
        return r - (r0 + 5.0)
    event.terminal = True
    event.direction = 1
    return event

def run_simulation(cfg: SimulationConfig):
    # Persiapkan direktori output
    data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data')
    results_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results')
    
    print(f"{'a':<10} | {'Lambda':<10} | {'Defleksi (deg)':<20} | {'Waktu Eksekusi'}")
    print("-" * 65)

    summary_data = []
    trajectories_for_plot = []

    for lmd in cfg.lambdas:
        t_start = time.time()
        
        spacetime = KerrDeSitter(cfg.M, cfg.a, lmd)
        integrator = HamiltonianIntegrator(spacetime)
        
        x0 = np.array([0.0, cfg.r0, cfg.theta0, cfg.phi0])
        try:
            pr_init = spacetime.get_initial_pr(x0, cfg.E, cfg.L_ang)
        except ValueError as err:
            logging.warning(f"Lewati L={lmd}: {err}")
            continue
            
        p0 = np.array([-cfg.E, pr_init, 0.0, cfg.L_ang])
        
        # Stop conditions (events)
        ev_horizon = make_event_horizon(cfg.rh)
        ev_escape = make_event_escape(cfg.r0)
        
        # Eksekusi ODE Solver C-Backend!
        sol = integrator.solve(x0, p0, tau_span=(0, cfg.tau_max), events=[ev_horizon, ev_escape])
        
        t_end = time.time()
        
        # Ekstraksi hasil
        t_vals = sol.y[0]
        r_vals = sol.y[1]
        phi_vals = sol.y[3]
        
        X_vals = r_vals * np.cos(phi_vals)
        Y_vals = r_vals * np.sin(phi_vals)
        
        # Hitung Hamiltonian untuk validasi (idealnya H=0 sepanjang garis)
        H_vals = np.array([spacetime.H(sol.y[:4, i], sol.y[4:, i]) for i in range(len(sol.t))])
        
        # Simpan ke CSV per Lambda
        df = pd.DataFrame({
            'tau': sol.t,
            'r': r_vals,
            'phi': phi_vals,
            'X': X_vals,
            'Y': Y_vals,
            'H': H_vals
        })
        csv_filename = os.path.join(data_dir, f'trajectory_lambda_{lmd}.csv')
        df.to_csv(csv_filename, index=False)
        
        # Hitung defleksi
        status = "Scattered" if sol.status == 1 else "Terminated"
        if len(phi_vals) > 0:
            deflection = abs(abs(phi_vals[-1] - cfg.phi0) - np.pi)
            def_deg = np.degrees(deflection)
            
            summary_data.append({
                'M': cfg.M,
                'a': cfg.a,
                'Lambda': lmd,
                'deflection_deg': def_deg,
                'status': status,
                'execution_time_sec': t_end - t_start
            })
            
            print(f"{cfg.a:<10} | {lmd:<10} | {def_deg:<20.4f} | {t_end - t_start:.4f} dtk")
            trajectories_for_plot.append({'lambda': lmd, 'X': X_vals, 'Y': Y_vals})

    # Simpan summary CSV
    df_summary = pd.DataFrame(summary_data)
    df_summary.to_csv(os.path.join(data_dir, 'deflection_summary.csv'), index=False)
    
    # ---------------------------------------------------------
    # GENERATE PLOT OTOMATIS KE FOLDER results/
    # ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8,8))
    colors = plt.cm.viridis(np.linspace(0, 1, len(trajectories_for_plot)))
    
    horizon = plt.Circle((0,0), cfg.rh, color='gray', alpha=0.3, label="Horizon")
    ax.add_artist(horizon)
    ax.scatter(0, 0, s=40, c='black', marker='x')

    for i, res in enumerate(trajectories_for_plot):
        label_txt = rf"$\Lambda = {res['lambda']}$"
        ax.plot(res['X'], res['Y'], color=colors[i], label=label_txt, linewidth=1.5)

    ax.set_title(f"Lintasan Sinar Kerr-De Sitter (M={cfg.M}, a={cfg.a}, b={cfg.b})", fontsize=14)
    ax.set_xlabel(r'$X = r \cos(\phi)$')
    ax.set_ylabel(r'$Y = r \sin(\phi)$')
    ax.set_aspect('equal')
    ax.set_xlim(-10, 10)
    ax.set_ylim(-10, 10)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    plot_path = os.path.join(results_dir, f'trajectories_{cfg.a}.png')
    plt.savefig(plot_path, dpi=300)
    logging.info(f"Plot berhasil disimpan di {plot_path}")
    # plt.show() # Dinonaktifkan untuk run headless yang bersih

if __name__ == '__main__':
    logging.info("Memulai simulasi SciPy Integrator End-to-End...")
    cfg = SimulationConfig()
    run_simulation(cfg)

