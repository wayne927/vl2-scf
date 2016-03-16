#!/usr/bin/python

# reads an SCF file and prints the nmax, lmax, and 
# the first N coefficients in the file
# Usage: read_scf.py filename N

import struct as st
import numpy as np
import sys
import os

filename = sys.argv[1]
num_coeffs = int(sys.argv[2])

file = open(filename, 'rb')

dat = file.read(8)

dat = st.unpack('2i', dat[0:8])
nmax = dat[0];
lmax = dat[1];

filestat = os.stat(filename)
filesize = filestat.st_size

print nmax, lmax

nbytes = -1

# Is this the SCF file for the main halo?
# 2*(nmax+1)*(lmax+1)^2 numbers in double precision
if (2*(nmax+1)*(lmax+1)**2*8+8 == filesize) :
    nbytes = 8
    suffix = 'd'

# Or is it one of the subhalos?
# (nmax+1) numbers in single precision
elif ((nmax+1)*4+8 == filesize) :
    nbytes = 4
    suffix = 'f'

else :
    print "Inconsistent SCF file! Can't print this."
    sys.exit(1)

dat = file.read(num_coeffs*nbytes)

dat = st.unpack(str(num_coeffs)+suffix, dat[:])


for num in dat:
    print "%e" % num

