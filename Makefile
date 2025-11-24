make: clean rust cython

rust: cleanrust
	maturin develop --release

cython: cleancython
	cd bioloptimizer && python3 -B setup.py build_ext --inplace

clean: cleanrust cleancython
	rm -rf lib

cleanrust:
	rm -rf target
	pip uninstall biorust -y

cleancython:
	rm -rf bioloptimizer/build
	rm -rf bioloptimizer/*.so
	rm -rf bioloptimizer/cythonlib/biocython.c
