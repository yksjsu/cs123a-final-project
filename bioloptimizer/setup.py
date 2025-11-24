from setuptools import setup, Extension
from Cython.Build import cythonize


extensions = [
    Extension(
        "biocython",
        sources=["cythonlib/biocython.pyx", "cythonlib/biocython_def.c"],
        include_dirs=["."],
    )
]

setup(
    ext_modules=cythonize(extensions, compiler_directives={"language_level": "3"})
)