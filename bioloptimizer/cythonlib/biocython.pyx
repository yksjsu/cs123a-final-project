cdef extern from "biocython.h":
    int add_numbers(int a, int b)
    int check_import()


def add_numbers(a: int, b: int) -> int:
    return add_numbers(a, b)

def check_import() -> int:
    return check_import()