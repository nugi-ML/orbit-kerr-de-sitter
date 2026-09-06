# Simulasi Lintasan Foton pada Ruang-Waktu Kerr-de Sitter

Repositori ini berisi implementasi komputasi numerik untuk menyimulasikan lintasan foton (*null geodesics*) di sekitar lubang hitam berotasi dengan pengaruh konstanta kosmologis (Metrik Kerr-de Sitter) dalam kerangka Teori Relativitas Umum Einstein. 

Kode dirancang dengan arsitektur modular, efisiensi komputasi tinggi menggunakan solver C-backend, pemisahan data berbasis CSV, serta visualisasi otomatis siap publikasi.

---

## 🌟 Fitur Utama

- **Pre-computed Symbolic Formulation:**
  Penurunan metrik $g_{\mu\nu}$, invers $g^{\mu\nu}$, Hamiltonian $\mathcal{H}$, serta turunan parsial $\frac{\partial \mathcal{H}}{\partial x}$ dan $\frac{\partial \mathcal{H}}{\partial p}$ diturunkan secara analitik menggunakan **SymPy**, lalu diekspor menjadi fungsi numerik murni berbasis **NumPy**. Pendekatan ini menghilangkan *overhead* aljabar simbolik saat simulasi berlangsung.

- **High-Performance Adaptive ODE Integrator:**
  Menggunakan `scipy.integrate.solve_ivp` dengan metode **RK45 (Dormand-Prince)** dan *adaptive step-size*. Integrasi melambat secara presisi di dekat *periapsis* (titik terdekat dengan lubang hitam) dan berakselerasi di medan asimtotik jauh, memberikan percepatan hingga **~145× lebih cepat** dibandingkan metode RK4 fixed-step standar Python murni.

- **Verifikasi Konservasi Hamiltonian:**
  Dilengkapi modul validasi nilai $\mathcal{H} \approx 0$ sepanjang lintasan untuk mengontrol akumulasi galat numerik dan memastikan foton bergerak tepat pada *null geodesic*.

- **Pemisahan Data & Kode (Clean Data Architecture):**
  Hasil komputasi langsung disimpan secara terstruktur di direktori `data/` dalam format `.csv` berdasarkan variasi parameter spin $a$, sehingga memori RAM tetap hemat dan analisis data dapat dilakukan secara terpisah di Jupyter Notebook.

- **Otomatisasi Plot Siap Publikasi:**
  Script simulasi langsung menghasilkan grafik lintasan beresolusi tinggi (300 DPI) ke direktori `results/` dengan penamaan file otomatis sesuai nilai spin lubang hitam.

---

## 📁 Struktur Direktori

```text
orbit-foton-kerr/
├── data/                         # Direktori penyimpanan dataset CSV
│   ├── a_0.0/                    # Data simulasi untuk spin a = 0.0 (Schwarzschild-de Sitter)
│   │   ├── deflection_summary.csv
│   │   └── trajectory_lambda_*.csv
│   ├── a_0.2/                    # Data simulasi untuk spin a = 0.2
│   ├── a_0.4/                    # Data simulasi untuk spin a = 0.4
│   ├── a_0.6/                    # Data simulasi untuk spin a = 0.6
│   ├── a_0.8/                    # Data simulasi untuk spin a = 0.8
│   ├── a_0.99/                   # Data simulasi untuk spin a = 0.99 (Extreme Kerr-de Sitter)
├── notebooks/                    # Eksplorasi dan analisis interaktif
│   ├── Analisis_Data.ipynb       # Notebook visualisasi dan analisis data CSV
│   ├── Penelitian_S1.ipynb       # Notebook eksperimen awal (arsip)
│   └── Penelitian_S1_Refactored.ipynb # Notebook terstruktur versi awal
├── results/                      # Output visualisasi grafik lintasan (PNG 300 DPI)
│   ├── trajectories_0.0.png
│   ├── trajectories_0.2.png
│   ├── trajectories_0.4.png
│   ├── trajectories_0.6.png
│   ├── trajectories_0.8.png
│   └── trajectories_0.99.png
├── src/                          # Modul utama program (Source Code)
│   ├── generate_equations.py     # Script penurun persamaan matematis SymPy
│   ├── kerr_equations.py         # File fungsi persamaan statis NumPy (auto-generated)
│   ├── physics.py                # Kelas representasi geometri ruang-waktu Kerr-de Sitter
│   ├── integrator.py             # Wrapper SciPy RK45 solver untuk sistem Hamiltonian
│   └── simulation.py             # Script eksekusi simulasi, ekspor data, & plotting
├── requirements.txt              # Daftar dependensi pustaka Python
└── README.md                     # Dokumentasi teknis proyek
```

---

## ⚙️ Persyaratan Sistem & Instalasi

### 1. Prasyarat
- Python 3.12 (direkomendasikan)
- Package manager: [uv](https://github.com/astral-sh/uv) (sangat direkomendasikan karena kecepatan instalasinya) atau `pip`.

### 2. Instalasi dengan `uv` (Rekomendasi)

```powershell
# Clone repositori / buka direktori proyek
cd c:\Users\nugip\Documents\orbit-foton-kerr

# Buat virtual environment berbasis Python 3.12
uv venv --python 3.12

# Aktifkan environment (PowerShell)
.\.venv\Scripts\Activate.ps1

# Instal seluruh dependensi
uv pip install -r requirements.txt
```

### 3. Instalasi dengan `pip` Standar

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

---

## 🚀 Panduan Penggunaan

### 1. Menjalankan Simulasi Utama
Untuk menjalankan simulasi pelacakan sinar, menghitung sudut defleksi foton untuk berbagai nilai konstanta kosmologis $\Lambda$, mengekspor data ke folder `data/`, dan menyimpan grafik ke `results/`:

```powershell
python src/simulation.py
```

> **Catatan Pengaturan Parameter:**  
> Anda dapat mengubah parameter massa $M$, parameter spin $a$, impact parameter $b$, atau daftar nilai $\Lambda$ langsung di dalam kelas `SimulationConfig` pada file [src/simulation.py](file:///c:/Users/nugip/Documents/orbit-foton-kerr/src/simulation.py):
> ```python
> @dataclass
> class SimulationConfig:
>     M: float = 1.0
>     a: float = 0.6            # Nilai spin lubang hitam
>     lambdas: List[float] = None  # Default: [0.0, -1e-8, -5e-8, ..., -5e-6]
>     b: float = 7.0            # Impact parameter
>     r0: float = 500.0         # Jarak awal sumber foton
> ```

### 2. Membaca & Menganalisis Data di Jupyter Notebook
1. Buka [notebooks/Analisis_Data.ipynb](file:///c:/Users/nugip/Documents/orbit-foton-kerr/notebooks/Analisis_Data.ipynb).
2. Pastikan memilih kernel Python dari environment `.venv`.
3. Tentukan variabel spin `a_val` di sel pertama:
   ```python
   a_val = 0.0  # Sesuaikan dengan folder data/a_{a_val} yang ingin dibaca
   ```
4. Jalankan sel untuk memuat tabel ringkasan defleksi dan membuat grafik konservasi Hamiltonian $\mathcal{H}$ vs radius $r$.

### 3. Memperbarui Persamaan Simbolik (Opsional)
Jika Anda memodifikasi bentuk analitik metrik dasar atau ingin memperbarui penurunan turunan:

```powershell
python src/generate_equations.py
```
*Proses ini memakan waktu $\sim$10-11 menit (sekali eksekusi saja) untuk menyelesaikan inversi simbolik $4 \times 4$ dan derivasi kalkulus rantai, lalu memperbarui file `src/kerr_equations.py`.*

---

## 📐 Landasan Teori & Formulasi Matematis

### 1. Metrik Kerr-de Sitter (Koordinat Boyer-Lindquist)
Elemen garis $ds^2 = g_{\mu\nu} dx^\mu dx^\nu$ di bidang ekuator/ruang-waktu Kerr-de Sitter didefinisikan dengan fungsi-fungsi metrik:

$$ds^2 = -\frac{\Delta_r}{\Sigma}\left(dt - \frac{a \sin^2\theta}{\Xi} d\phi\right)^2 + \frac{\Sigma}{\Delta_r} dr^2 + \frac{\Sigma}{\Delta_\theta} d\theta^2 + \frac{\Delta_\theta \sin^2\theta}{\Sigma}\left(a dt - \frac{r^2 + a^2}{\Xi} d\phi\right)^2$$

dengan fungsi-fungsi metrik:

$$\Sigma = r^2 + a^2 \cos^2\theta$$

$$\Delta_r = (r^2 + a^2)\left(1 - \frac{\Lambda r^2}{3}\right) - 2Mr$$

$$\Delta_\theta = 1 + \frac{\Lambda a^2 \cos^2\theta}{3}$$

$$\Xi = 1 + \frac{\Lambda a^2}{3}$$

### 2. Dinamika Geodesik Foton (Hamiltonian Formalism)
Lintasan foton (*null geodesic*) memenuhi kondisi super-Hamiltonian:

$$\mathcal{H}(x^\mu, p_\mu) = \frac{1}{2} g^{\mu\nu}(x) p_\mu p_\nu = 0$$

Persamaan gerak kanonik Hamilton:

$$\frac{dx^\mu}{d\tau} = \frac{\partial \mathcal{H}}{\partial p_\mu}$$

$$\frac{dp_\mu}{d\tau} = -\frac{\partial \mathcal{H}}{\partial x^\mu}$$

dengan $\tau$ merupakan parameter afin lintasan foton.

### 3. Kondisi Awal & Sudut Defleksi
Untuk foton yang datang dari jarak jauh $r_0$ dengan energi $E$ dan parameter impak $b$ ($L = bE$), momentum radial awal $p_r$ ditentukan dari constraint $\mathcal{H} = 0$:

$$p_r = -\sqrt{-\frac{g^{tt}E^2 + g^{\phi\phi}L^2 - 2g^{t\phi}EL}{g^{rr}}}$$

Sudut defleksi $\hat{\alpha}$ dihitung dari perubahan azimuth total foton yang lolos (*scattering*):

$$\hat{\alpha} = |\Delta\phi - \pi| = |\phi_{\text{out}} - \phi_{\text{in}} - \pi|$$

---

## 📊 Hasil Komputasi

Contoh hasil defleksi sudut untuk spin $a = 0.8$ dan berbagai variasi konstanta kosmologis $\Lambda$:

| $a$ | $\Lambda$ | Defleksi ($^\circ$) | Waktu Komputasi |
| :---: | :---: | :---: | :---: |
| 0.8 | $0$ | 46.4746 | 7.92 dtk |
| 0.8 | $1 \times 10^{-8}$ | 46.4748 | 7.88 dtk |
| 0.8 | $5 \times 10^{-8}$ | 46.4755 | 7.96 dtk |
| 0.8 | $1 \times 10^{-7}$ | 46.4763 | 7.73 dtk |
| 0.8 | $5 \times 10^{-7}$ | 46.4831 | 8.17 dtk |
| 0.8 | $1 \times 10^{-6}$ | 46.4918 | 8.28 dtk |
| 0.8 | $5 \times 10^{-6}$ | 46.5730 | 7.52 dtk |

Grafik lintasan yang dihasilkan dapat dilihat pada direktori [results/](file:///c:/Users/nugip/Documents/orbit-foton-kerr/results/).

---

## 🛠️ Stack Teknologi

- **Bahasa:** Python 3.12
- **Matematika Simbolik:** [SymPy](https://www.sympy.org/)
- **Komputasi Numerik & Matriks:** [NumPy](https://numpy.org/)
- **Penyelesaian Persamaan Diferensial (ODE):** [SciPy (`scipy.integrate.solve_ivp`)](https://scipy.org/)
- **Struktur Data & Analisis:** [Pandas](https://pandas.pydata.org/)
- **Visualisasi & Plotting:** [Matplotlib](https://matplotlib.org/)
- **Lingkungan Interaktif:** [Jupyter / IPython Kernel](https://jupyter.org/)

