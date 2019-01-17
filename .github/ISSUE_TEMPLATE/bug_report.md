---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. logs or outputs
```
terminal outputs
```
2. namelist
```
noahmp.namelist
```
3. the domain file headers
```
ncdump -h domain.nc
```
4. the initialization or restart file headers
```
ncdump -h init.nc or ncdump -h restart.nc
```
5. the input file headers
```
ncdump -h input.20000101T000000
```
6. parameter tables
default or enter your update here

**Expected behavior**
A clear and concise description of what you expected to happen.

**Environment (please complete the following information):**
 - OS: [e.g. Linux, Windows]
 - Parallel environment: [e.g. MPI, OpenMP, CUDA]

**Additional context**
Add any other context about the problem here.
