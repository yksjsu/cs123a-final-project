# Bioloptimizer

## Overview

The goal of this project is to achieve runtime comparable to the Biopython Python package while also providing a visualization feature.

## Setup

Please ensure that you have a recent version of Python installed (i.e. Python 3.12.x, Python 3.13.x, Python 3.14.x, etc).

Create a new virtual environment called ".env" and activate the environment.

Install the pip packages in the requirements.txt file.

## Running Instructions

Use ```make clean``` to clean build files, and ```make rust``` to build the necessary packages for rust.

Use ```python3 -B -m bioloptimizer --help``` as a way to run the CLI interface's help message.

If you wish to try a sample command for running global, linear alignment, try this command:

```bash
python3 -B -m bioloptimizer linalign ACTGTTGT ACTC --no-is-local
```

If you wish to try a sample command for running global, affine alignment, try this command:

```bash
python3 -B -m bioloptimizer affalign ACTGTTGT ACTC --no-is-local
```

These commands use default match, mismatch, and linear gap / gap open penalty / gap extension penalty scores.