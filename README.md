# \# MeshSim (C++) + buffer\_certify tool

# 

# This repository contains two programs:

# 

# 1\) \*\*MeshSim (C++)\*\*  

# &#x20;  A fast simulation + analytical model implementation with a standalone CLI:

# &#x20;  - `meshsim\_cli`

# 

# 2\) \*\*buffer\_certify\*\* (\*\*separate tool\*\*)  

# &#x20;  Included for the specific “buffer certification” workflow.  

# &#x20;  \*\*It is not part of MeshSim / MeshSorter\*\*; it simply lives in this repository for convenience.

# 

# \---

# 

# \## Files

# 

# Core MeshSim:

# \- `MeshSim.hpp`

# \- `MeshSim.cpp`

# 

# MeshSim CLI:

# \- `MeshSim\_cli.cpp`  → builds `meshsim\_cli`

# 

# Buffer certify tool:

# \- `buffer\_hypothesis\_certify\_3phases\_iut.cpp` → builds `buffer\_certify`

# 

# Build:

# \- `CMakeLists.txt`

# 

# \---

# 

# \## Build

# 

# \### Linux (Ubuntu/Debian)

# 

# Install prerequisites:

# 

# ```bash

# sudo apt update

# sudo apt install -y cmake build-essential

