#!/bin/awk -f

BEGIN {
  print "V 2000.1.2"
  print "TBEGIN ? ? ?"

  srand(1);
  for(i=0; i<1; i++){
    x=(2*rand()-1)*500 # meters
    y=(2*rand()-1)*500 # meters
    z=(2*rand()-1)*500 # meters
    zenith=rand()*180  # degrees
    azimuth=rand()*360 # degrees
    l=500              # length, m
    energy=5        # GeV
    t=0                # ns
    N=1    # number of photons

    print "EM 1 1 1970 0 0 0"
    print "TR 1 0 e-    ", x, y, z, zenith, azimuth, 0, energy, t
  ##  print "TR 1 0 amu  ", x, y, z, zenith, azimuth, l, energy, t
 ##   print "TR 1 0 hadr ", x, y, z, zenith, azimuth, 0, energy, t
 ##   print "TR 1 0 ph   ", x, y, z, zenith, azimuth, 0, N, t
  ##  print "EE"
  }
  print "TEND ? ? ?"
  print "END"
}