#include "biocython.h"
#include "Python.h"


int check_import() {
    printf("Cython module imported!\n");
    return 0;
}


int add_numbers(int a, int b) {
    return a + b;
}

