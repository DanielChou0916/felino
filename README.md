# Felino

Felino is an open-source phase-field fracture framework built on the [MOOSE framework](https://mooseframework.inl.gov/).  
This repository provides an extension for geomaterial fracture modeling, supporting asymmetric tensile/compressive behavior.

---

## 🚀 Installation Guide

The following instructions assume you are using **Ubuntu (or WSL with Ubuntu)**.  
Felino is built on top of the MOOSE framework. The instruction in this page includes installation of MOOSE and Felino.
And it is presented with the Assumption that the number of cpu on user's device is 4.
### 1️⃣ Install Miniforge (Lightweight Conda)

```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge
```

### 2️⃣ Setup Conda and Restart Shell

```bash
export PATH=$HOME/miniforge/bin:$PATH
conda init --all
exit
```

⚠️ Close and reopen your terminal after running `exit`.

---

### 3️⃣ Update Conda and Add MOOSE Channel

```bash
conda update --all --yes
conda config --add channels https://conda.software.inl.gov/public
```

---

### 4️⃣ Create and Activate MOOSE Environment

```bash
conda create -n moose moose-dev=2025.07.22=mpich
conda activate moose
```

---

### 5️⃣ Clone and Checkout MOOSE

```bash
mkdir -p ~/projects
cd ~/projects
git clone https://github.com/idaholab/moose.git
cd moose
git checkout master
```

---

### 6️⃣ Build and Test MOOSE

```bash
cd ~/projects/moose/test
make -j 4
./run_tests -j 4
```

✅ If you see `All tests passed`, MOOSE installation is successful.

---

## 🐾 Felino Installation

After MOOSE is installed, build Felino:

```bash
cd ~/projects
git clone https://github.com/DanielChou0916/felino.git
cd felino
make -j 4
```

---

## ▶️ Running a Test

```bash
cd ~/projects/felino
chmod +x felino_test.sh
./felino_test.sh
```

When you see `[Felino] (='v'=) Installation test PASSED.`, the test has run successfully.


---
## ▶️ Running a Tutorial for v1.0.2 Geomaterial Analysis
📌 This example demonstrates the new geomaterial fracture models added in v1.0.2.
```bash
cd ~/projects/felino
mpiexec -n 4 ./felino-opt -i tutorials/tur4_uniaxial_compression/Ce/Ce_disp.i
```

---

## 🔄 Updating Felino

```bash
cd ~/projects/felino
git pull
make clean
make -j 4
```

---

## 📖 Documentation & Tutorials

📌 Official website: [Felino Website](https://danielchou0916.github.io/felino.github.io/)  
📌 Tutorials and benchmark examples are provided under the `tutorials/` folder.

## Citation

If you use **Felino** in your work, please cite:

> Chou, D. T. (2026). *Felino: Extension of an open-source phase-field framework to geomaterial fracture*. **Journal of Open Source Software, 11**(121), 9406. https://doi.org/10.21105/joss.09406

```bibtex
@article{chou2026felino,
  title={Felino: Extension of an open-source phase-field framework to geomaterial fracture},
  author={Chou, Daniel T.},
  journal={Journal of Open Source Software},
  volume={11},
  number={121},
  pages={9406},
  year={2026},
  doi={10.21105/joss.09406}
}
