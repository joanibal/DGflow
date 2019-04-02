default:
# Check if the config.mk file is in the config dir.
	f2py -c --f90flags='-g -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan' -m dg_solver DGSolver.f90

production:
	f2py -c --f90flags='-O3 -ffast-math -funroll-loops' -m dg_solver DGSolver.f90
