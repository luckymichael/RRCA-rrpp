gfortran -static -fdefault-real-8 -O2 rrpp.f krige.f utl.f -o rrpp
gfortran -static -fdefault-real-8 -O2 acct.f utl.f -o acct
gfortran -static -fdefault-real-8 -O2 acct32.f90 utl.f -o acct32