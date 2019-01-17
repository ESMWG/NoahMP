# The Noah with Multi-Parameterization options (NoahMP) land surface model

This readme file describes the building and running processes of the NoahMP land surface model.

## Building
To build the NoahMP land surface model, users need to set environment variables, and then configure and compile the model as follows.

### Environmental variables

Environmental variables are used to specify the NetCDF installation path.

You can only set the `NETCDF` environment variable:
```bash
export NETCDF=$your_netcdf_installation_directory
```
In this case, the NetCDF develoment header files and runtime libraries will be assumed under "`$NETCDF/include`" and "`$NETCDF/lib`", respectively.

If these settings do not match your configurations, you have to specify  the NetCDF include and library pathes separately:
```bash
export NETCDF_INC=$your_netcdf_include_directory
export NETCDF_LIB=$your_netcdf_library_directory
```

NoahMP uses two NetCDF libraries for Fortran and C respectively: 
`libnetcdff` and `libnetcdf`. Make sure that these two files are in your NetCDF library path. If the user's NetCDF library combined them together (only has one), the user will need to manually change this part in order to successfully compile NoahMP. 
See the section below on porting about how to change this.

Notes: If you are going to create model output file that is more than 2Gb,
      you should consider using NetCDF large file support function. To activate
      this, one must set the environment variable `WRFIO_NCD_LARGE_FILE_SUPPORT`.
      In a bash environment, do
      
```bash
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
```

### Configuring

Run the configuration script to configure the compiling:
```bash
./configure
```
Several options will be prompted. Select the one that match your operating system (Linux, Mac OS X Darwin), compiler (GCC/Gfortran, Intel, PGI), and intended parallel environment (seq for sequential and dm for MPI).

After the configuration, a file named `makefile.in` will be generated in the source code directory, specifying the detailed compiling settings. Modify if necessary.

### Compiling

Run `make` to start the compilation:
```bash
make
``` 

If compiled successfully, an executable file named `noahmp.exe` will be created under "run" directory.

## Running

The NoahMP model use a namelist called `noahmp.namelist` as well as some additional parameter files (.TBL files), which are located under the "run" directory. Users need to copy those files to the directory where the model is going to be executed.

For a NoahMP cold start run (i.e. not from a restart file), the user needs to turn off the flag `from_restart=.false.` in `noahmp.namelist` and provide an initialization file to the `INIT_FILE` option.

For running NoahMP from restart file, the user needs to switch on the flag `from_restart=.true.` in `noahmp.namelist` and provide an existing restart file to the option `RESTART_FILE`. Running from a restart condition is common when the land surface has been 
'spun-up' by running NoahMP in an offline or 'uncoupled' capacity.

## Debugging

There are two debugging methods shown as follow. You can combine the two methods or use separately. Both methods require recompiling.

1) Print diagnostic information at run time.

Set the environment variable `HYDRO_D=1` and rebuild the model.

```bash
export HYDRO_D=1
```
A "1" for `HYDRO_D` results in NoahMP producing some run-time diagnostic information. 
When `HYDRO_D` is set to "0 "or not defined, the diagnostic information will not be produced 
during run-time. 

2) Compile a debugable executable for use by debuggers.

Edit "`makefile.in`", and append the compiler debug options "`-g`" to the fortran compiler "`F90`".

## Example

An example presenting the namelist, domain, input, restart files can be found in [https://github.com/esmwg/NoahMP-Training].

## Porting

The NoahMP model does not presently support OpenMP. The default support platform is Linux 
with the Intel compiler and Linux with the GFortran compiler. However, NoahMP is fairly easy to port to other systems.  
The basic steps to do so are as follows:

1) Provide your own `makefile.in` under the "arch" directory.
2) Edit `configure`, add your own compling options in the "Compiling Option" section.
3) Edit `configure`, add the path of your own `makefile.in` in the "Compiling Configurations" section.
4) ./configure
5) make

If there is no error, then users can compile  NoahMP on the new platform.

## Contribution

Fork this repository, code, and submit a pull request.

Read the GitHub help page [https://help.github.com/articles/fork-a-repo] for more information.
