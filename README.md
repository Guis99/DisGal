# DisGal
I (finally) take the time to learn Discontinuous Galerkin

# Notes

## Prologue: Testing MathJax

Ladyzhenskaya–Babuška–Brezzi condition:
$$ \underset{q\in Q}{\text{inf}} \,\, \underset{v\in V}{\text{sup}} \frac{b(v,q)}{\| v_{V} \| \| q_{Q} \|} \geq \beta $$

# Running locally

Requirements:
- A supported C++ compiler
- CMake (non-windows)
- Make (non-windows)
- A working Python installation
    - Python's pip package manager

At the moment, only the code in the "Diffusion" subdirectory is available to run. First, clone this repository via a method of your choice. Once done, run 

```
cd DisGal
git submodule update --init --recursive
```

from the terminal to load dependencies. 

Then, please follow the below platform-specific instructions:

# Platform-specific instructions
## Linux (x86) and Mac (ARM)
- Run the following commands in your terminal:

```
cd Diffusion
mkdir build
```

At this point, if you are running Linux you can choose to enable/disable multithreading. To do so, open CMakeLists.txt and either remove or keep the '-fopenmp' and '-DMULTITHREAD' flags on line 18.

```
cd build
cmake ..
make
cd ..
```

- Verify that there is a binary called "diffDG" in the build folder
- Skip to the next section below

## Windows
- Run the following commands in your terminal:

```
cd Diffusion
mkdir build
cd ..
```

At this point, you can choose to enable/disable multithreading. To do so, open CMakeLists.txt and either remove or keep the '-fopenmp' and '-DMULTITHREAD' in the batch file called 'buildPoisson_dG.bat'

```
buildPoisson_dG.bat
cd Diffusion
```

- Verify that there is a binary called "diffDG.exe" in the build folder

# Instructions (continued)
- Make sure your current directory is DisGal/Diffusion
- Set up Python environment
    - Make sure pip is installed
```
pip -r install ../requirements.txt
```   
- Open plot1.py
- Run plot1.py
- Try changing the parameters. Set different boundary conditions, refine the resolution, pass more threads to the solver, ... Have some fun with it!
