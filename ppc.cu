#include <map>
#include <deque>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sys/time.h>

#ifndef __CUDACC__
#define XCPU
#endif

#ifdef XCPU
#include <cmath>
#include <cstring>
#endif



using namespace std;

namespace xppc{

#include "ini.cxx"
#include "pro.cu"  // You might need to adapt this later

  void initialize(float enh = 1.f){ m.set(); q.eff*=enh; }

  unsigned int pmax, pmxo, pn, pk, hquo;

  void setq(){
    char * HQUO=getenv("HQUO");
    hquo=HQUO==NULL?1:atoi(HQUO);
    cerr<<"HQUO(photons/max number of hits)="<<hquo<<endl;
  }

  // Keep only this XCPU section, remove everything in the #else
  dats *e;
  int nblk, nthr, ntot;

  void ini(){
    setq();
    rs_ini();
    pn=0; pk=0;

    ntot=nblk*nthr;
    pmax=ntot*NPHO;
    pmxo=pmax/OVER;
    pmax=pmxo*OVER;
    d.hnum=pmax/hquo;

    d.gdev=0; d.gnum=1;
    d.gini=0; d.gspc=pmax; d.gtot=pmax; d.gdiv=1;

    {
      d.hits = q.hits = new hit[d.hnum];
      d.pz = q.pz = new photon[pmxo];
      d.bf = new pbuf[pmax];
    }

    {
      d.z=&z; d.oms=q.oms; e=&d;
    }

    {
      unsigned int size=d.rsize, need=seed+1;
      if(size<need) cerr<<"Error: not enough multipliers: asked for "<<seed<<"-th out of "<<size<<"!"<<endl;
    }
  }

  void fin(){
    delete d.pz;
    delete d.hits;
    delete d.bf;
  }


  void start(){}
  void stop(){}
  void choose(int device){
    sv+=device;
    seed=device;
    nblk=NBLK, nthr=NTHR;
  }
  void listDevices(){}

#include "f2k.cxx"
}

#ifndef XLIB
using namespace xppc;

float zdh;

float zshift(float4 r){
  zdh=d.dh;
  return zshift(d, r, zdh);
}

// void initialize(){
//     m.set();
// }

float photon_yield(string loss_type, int energy_gev, float track_length) {
    float rho = 0.9216f; // density of ice in icecube
    float logE = logf(energy_gev);
    float em_cascade_value=5.321*0.910f/rho;  // important for em cascade photons
    float eff_tracl_length = 0;
    float num_photons = 0.0f;
    float ppm = 2000.0f;
    if (loss_type == "amu-") {
        // calcualte photons for muon track
        float additional_track = 1+ max(0.0f, 0.1880f+0.0206f*logE)*0.910f/rho;
        num_photons = track_length>0?track_length*additional_track:0;
    } 
    if (loss_type == "em") {
        num_photons=energy_gev*em_cascade_value;
    }
    return num_photons*ppm;
}

int main(int argc, char* argv[]) {
  start(); // doesnt do anythin i think
    // Arguements: loss type str, energy gev, track length m
    cout << "initalizing: " << endl;
    initialize();
    string loss_type= argv[1];
    int energy_gev = stoi(argv[2]);
    float track_length = stof(argv[3]);

    cout << "loss type: " << loss_type << ", has energy: " << energy_gev << ", and track length: " << track_length << endl;
    float photons = photon_yield(loss_type, energy_gev, track_length);
    cout << "photons: " << photons << endl;

    return 0;
}
#endif

