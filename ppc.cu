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

#ifdef USE_I3_LOGGING
#include "icetray/I3Logging.h"
#else
#define log_info_stream(msg) \
  do { std::cerr << msg << std::endl; } while (0)
#endif

using namespace std;

namespace xppc{
  #include "ini.cxx"
  #include "pro.cu"

  void initialize(float enh = 1.f){ m.set(); q.eff*=enh; }

  unsigned int pmax, pmxo, pn, pk, hquo;

  void setq(){
    char * HQUO=getenv("HQUO");
    hquo=HQUO==NULL?1:atoi(HQUO);
    cerr<<"HQUO(photons/max number of hits)="<<hquo<<endl;
  }

  dats *e;  // pointer to a copy of "d" on device
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




  void print();

  void kernel(unsigned int num){
    unsigned int & old = num;
    if(old>0){
      d.hidx=0;
      for(d.blockIdx=0, d.gridDim=nblk, blockDim.x=nthr; d.blockIdx<d.gridDim; d.blockIdx++)
	    for(threadIdx.x=0; threadIdx.x<blockDim.x; threadIdx.x++) propagate(e, num);  
      if(d.hidx>d.hnum){ cerr<<"Error: data buffer overflow occurred: "<<d.hidx<<">"<<d.hnum<<"!"<<endl; d.hidx=d.hnum; }
      log_info_stream("photons: "<<old<<"  hits: "<<d.hidx);
    }



  if(old>0) print();

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


using namespace xppc;

float zdh;

float zshift(float4 r){
  zdh=d.dh;
  return zshift(d, r, zdh);
}

int main(int arg_c, char *arg_a[]){
  start();
  if(arg_c<=1){
    listDevices();
    fprintf(stderr, "Use: %s [device] (f2k muons)\n"
	    "     %s [str] [om] [num] [device] (flasher)\n", arg_a[0], arg_a[0]);
  }
  else if(0==strcmp(arg_a[1], "-")){
    initialize();
    ices & w = z.w[WNUM/2];
    cerr<<"For wavelength="<<q.wvs[w.wvl].w<<" [nm]  np="<<(1/w.coschr)<<"  cm="<<1/w.ocm<<" [m/ns]"<<endl;
    float4 r;
    r.w=0;
    if(arg_c==4){
      r.x=atof(arg_a[2]);
      r.y=atof(arg_a[3]);
    }
    else r.x=0, r.y=0;
    for(int i=0; i<d.size; i++){
      float z=d.hmin+d.dh*i;
      r.z=z; for(int j=0; j<10; j++) r.z=z+zshift(r); z=r.z;
      cout<<z<<" "<<w.z[i].abs<<" "<<w.z[i].sca*(1-d.g)<<" "<<d.az[i].ra*d.sum<<endl;
    }
  }
  else if(0==strcmp(arg_a[1], "=")){
    initialize();
    ices & w = z.w[WNUM/2];
    cerr<<"For wavelength="<<q.wvs[w.wvl].w<<" [nm]  np="<<(1/w.coschr)<<"  cm="<<1/w.ocm<<" [m/ns]"<<endl;
    float4 r;
    r.w=0;
    string in;
    while(getline(cin, in)){
      if(3==sscanf(in.c_str(), "%f %f %f", &r.x, &r.y, &r.z)){
	float dz=zshift(r);
	cout<<in<<" "<<dz<<" "<<zdh<<endl;
      }
    }
  }
  else if(0==strcmp(arg_a[1], "_")){
    initialize();
    float4 r;
    r.w=0;
    for(r.x=-750.f; r.x<751.f; r.x+=3.f) for(r.y=-750.f; r.y<751.f; r.y+=3.f) for(float z=-750.f; z<751.f; z+=6.f){
	  r.z=z; for(int j=0; j<10; j++) r.z=z+zshift(r);
	  cout<<z<<" "<<r.x<<" "<<r.y<<" "<<(r.z-z)<<endl;
	}
  }
  else if(arg_c<=2){
    int device=0;
    if(arg_c>1) device=atoi(arg_a[1]);
    initialize();
    choose(device);
    fprintf(stderr, "Processing f2k muons from stdin on device %d\n", device);
    f2k();
  }
  else{
    int str=0, dom=0, device=0, itr=0;
    unsigned long long num=1000000ULL;

    if(arg_c>1) str=atoi(arg_a[1]);
    if(arg_c>2) dom=atoi(arg_a[2]);
    if(arg_c>3){
      num=(unsigned long long) atof(arg_a[3]);
      char * sub = strchr(arg_a[3], '*');
      if(sub!=NULL) itr=(int) atof(++sub);
    }
    if(arg_c>4) device=atoi(arg_a[4]);
    initialize();
    choose(device);
    fprintf(stderr, "Running flasher simulation on device %d\n", device);
    flasher(str, dom, num, itr);
  }

  stop();
}