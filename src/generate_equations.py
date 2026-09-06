import sympy as sp
from sympy.printing.lambdarepr import NumPyPrinter
import time
import os

def generate_equations():
    t0 = time.time()
    print("Mendefinisikan simbol...")
    t, r, th, ph = sp.symbols('t r theta phi')
    pt, pr, pth, pph = sp.symbols('pt pr ptheta pphi')
    M, a, L = sp.symbols('M a Lambda')

    Sigma = r**2 + (a * sp.cos(th))**2
    Delta_r = (r**2 + a**2) * (1 - (L * r**2)/3) - 2*M*r
    Delta_th = 1 + (L * a**2 * sp.cos(th)**2)/3
    Xi = 1 + (L * a**2)/3

    print("Membangun matriks metrik...")
    g_tt = - (Delta_r / Sigma) + (Delta_th * sp.sin(th)**2 / Sigma) * a**2
    g_rr = Sigma / Delta_r
    g_thth = Sigma / Delta_th
    g_pp = -(Delta_r/Sigma) * (a * sp.sin(th)**2 / Xi)**2 + \
           (Delta_th * sp.sin(th)**2 / Sigma) * ((r**2 + a**2)/Xi)**2
    g_tp = -(Delta_r/Sigma) * (-a * sp.sin(th)**2 / Xi) + \
           (Delta_th * sp.sin(th)**2 / Sigma) * (-a * (r**2 + a**2)/Xi)

    g_mat = sp.Matrix([
        [g_tt, 0, 0, g_tp],
        [0, g_rr, 0, 0],
        [0, 0, g_thth, 0],
        [g_tp, 0, 0, g_pp]
    ])
    
    print("Menghitung invers metrik (proses komputasi berat)...")
    g_inv = g_mat.inv()

    print("Merumuskan Hamiltonian...")
    p_vec = sp.Matrix([pt, pr, pth, pph])
    Ham_expr = sp.simplify(0.5 * (p_vec.T * g_inv * p_vec)[0])

    print("Menghitung turunan parsial terhadap x dan p...")
    vars_x = [t, r, th, ph]
    vars_p = [pt, pr, pth, pph]
    dHdx_expr = [sp.diff(Ham_expr, v) for v in vars_x]
    dHdp_expr = [sp.diff(Ham_expr, v) for v in vars_p]
    
    # Komponen inverse metrik yang dibutuhkan untuk set initial condition foton
    gUU_tt = sp.simplify(g_inv[0,0])
    gUU_rr = sp.simplify(g_inv[1,1])
    gUU_tp = sp.simplify(g_inv[0,3])
    gUU_pp = sp.simplify(g_inv[3,3])

    print("Mengekspor persamaan ke src/kerr_equations.py...")
    printer = NumPyPrinter()

    file_path = os.path.join(os.path.dirname(__file__), 'kerr_equations.py')
    with open(file_path, 'w', encoding='utf-8') as f:
        f.write("# FILE INI DI-GENERATE OTOMATIS OLEH src/generate_equations.py\n")
        f.write("# JANGAN MENGEDIT FILE INI SECARA MANUAL\n\n")
        f.write("import numpy as np\n\n")
        
        args = "t, r, theta, phi, pt, pr, ptheta, pphi, M, a, Lambda"
        
        # Penanganan khusus jika numpy tidak dipakai oleh pycode, diubah ke numpy explicit.
        f.write(f"def eval_H({args}):\n")
        f.write(f"    # Gunakan np dari modul numpy\n")
        # Replace global imports if printer output relies on internal math module sometimes
        f.write(f"    return {printer.doprint(Ham_expr).replace('numpy.', 'np.')}\n\n")
        
        f.write(f"def eval_dHdx({args}):\n")
        f.write(f"    return np.array([\n")
        for expr in dHdx_expr:
            f.write(f"        {printer.doprint(expr).replace('numpy.', 'np.')},\n")
        f.write(f"    ], dtype=np.float64)\n\n")
        
        f.write(f"def eval_dHdp({args}):\n")
        f.write(f"    return np.array([\n")
        for expr in dHdp_expr:
            f.write(f"        {printer.doprint(expr).replace('numpy.', 'np.')},\n")
        f.write(f"    ], dtype=np.float64)\n\n")
        
        f.write(f"def eval_gUU(t, r, theta, phi, M, a, Lambda):\n")
        f.write(f"    return (\n")
        f.write(f"        {printer.doprint(gUU_tt).replace('numpy.', 'np.')},\n")
        f.write(f"        {printer.doprint(gUU_rr).replace('numpy.', 'np.')},\n")
        f.write(f"        {printer.doprint(gUU_tp).replace('numpy.', 'np.')},\n")
        f.write(f"        {printer.doprint(gUU_pp).replace('numpy.', 'np.')}\n")
        f.write(f"    )\n")

    t1 = time.time()
    print(f"Selesai dalam {t1-t0:.2f} detik.")

if __name__ == '__main__':
    generate_equations()

