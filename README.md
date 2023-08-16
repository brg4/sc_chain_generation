# sc\_chain\_generation

## Description

**A fortran program to generate circular chain of spherocylinder segments within a spherical volume and in the presence of spherical crowders.**

Gilbert et al., *Frontiers in Cell and Developmental Biology*, 2023, DOI:[10.3389/fcell.2023.1214962](https://doi.org/10.3389/fcell.2023.1214962)

## Repository directory structure

 - `src` - *source files and Makefile for program*
 - `example` - *example to demonstrate program*

## Installation

Program installation requires a fortran compiler. The program was tested with gfortran-v12.1.0.

Run `make gen_sc_chain` in the `src` directory. The executable (**gen_sc_chain**) will be in `src`.

## Usage

Prepare an input file.

Run with: `./gen_sc_chain --i_f=<input file> --o_d=<output directory> --o_l=<output label> --s=<RNG seed> --l=<log file> --n_t=<num threads> --bin --dat --xyz`

The last three flags are optional and dictate what subsets of output files should be written.

- `bin` - *binary files storing array of doubles*
- `dat` - *plain-text files storing array as comma-separated values*
- `xyz` - *xyz files of atomic coordinates for simple visualization*

See `example` directory for a demonstration.

## Support
brg4@illinois.edu

## Authors and acknowledgment
Benjamin R. Gilbert - brg4@illinois.edu

## Project status
This project is under development.
