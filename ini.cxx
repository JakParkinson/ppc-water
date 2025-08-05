#define LMAX 80      // number of dust loggers
#define LYRS 172     // number of depth points

#define CTX  11
#define CTY  10
#define DIR1 9.3f
#define DIR2 129.3f

#define CX   21
#define CY   19
#define NSTR 240

#ifdef XCPU
#define OVER 1
#define NBLK 1
#define NTHR 512
#else
#define OVER 128     // minimum number of photons per particle (optimized)
#endif

#define NPHO   1024  // maximum number of photons propagated by one thread

#define WNUM   512   // number of wavelength slices
#define MAXLYS 172   // maximum number of ice layers
#define MAXGEO 16384 // maximum number of OMs
#define MAXRND 131072 // max. number of random number multipliers

// #define DTMN
#define XXX 1.e-5f
#define FPI 3.141592653589793f
#define OMR 0.16510f  // DOM radius [m]

static const float doma=FPI*OMR*OMR; // DOM cross-sectional area [m^2]
static const float omav=0.335822;    // average DOM angular (lab) sensitivity
static const float dppm=2450.08;     // photons per meter for nominal IceCube DOM
static const float fcv=FPI/180.f;

static unsigned int ovr=1;

float xrnd();

struct DOM{
  float R, F;
  float r[3];
};

struct ikey{
  int str, dom;

  bool isinice() const{
    return str>86 || (str>0 && dom>=1 && dom<=60);
  }

  bool operator< (const ikey & rhs) const {
    return str == rhs.str ? dom < rhs.dom : str < rhs.str;
  }

  bool operator!= (const ikey & rhs) const {
    return str != rhs.str || dom != rhs.dom;
  }
};

struct OM:DOM,ikey{};
vector<OM> i3oms;
map<int, pair<float, float> > strs;

template <int n> class V{
  float x[n];

public:
  V(){
    for(int i=0; i<n; i++) x[i]=i<n-1?0:1;
  }

  float & operator[](int i){
    return x[i];
  }

  float dot(V<n> q){
    float sum=0;
    for(int i=0; i<n; i++) sum+=x[i]*q[i];
    return sum;
  }
};

struct name:ikey{
  int omt;  // module type
  int type; // wavelength behavior
  float rde, hv;
  float azi; // cable
  V<3> tilt;

  name(){}
  name(ikey k, int m, int t, float r, float h, V<3>& axis, float ph):ikey(k){
    omt=m; type=t; rde=r; hv=h;
    tilt=axis; azi=ph;
  }
};

struct mesh{
  vector< V<3> > dirs;

  int ini(const string & file){
    ifstream inFile(file.c_str(), ifstream::in);
    if(!inFile.fail()){
      string in;
      while(getline(inFile, in)){
	V<3> dir;
	if(3==sscanf(in.c_str(), "%*d %f %f %f", &dir[0], &dir[1], &dir[2])){
	  dirs.push_back(dir);
	}
      }
      inFile.close();
      cerr<<"Configured "<<dirs.size()<<" uniform directions"<<endl;
    }
    return dirs.size();
  }
} ico;

bool nextgen=false;

struct itype{
  float area, beta, rde, fx, Rr, Rz, cable;
  vector< V<3> > dirs;

  bool def;
  float ave;  // average angular sensitivity
  float mas;  // maximum angular sensitivity
  vector<float> s; // ang. sens. coefficients

  itype(): def(false){ }

  void add(string file){
    def=true; mas=1, ave=0;
    area=1, beta=0.33f, Rr=OMR, Rz=OMR, cable=0;

    ifstream inFile(file.c_str(), ifstream::in);
    if(!inFile.fail()){
      float aux;
      if(inFile >> aux) mas=aux;
      while(inFile >> aux) s.push_back(aux);
      if(s.empty()) s.push_back(1.f);

      if(mas>0){
	ave=s[0];
	for(unsigned int i=2; i<s.size(); i+=2) ave+=s[i]/i;
      }
      else{
	ave=(1-s[0])/2;
      }
      if(ave>0) cerr<<"Loaded "<<(s.size()+1)<<" angsens coefficients"<<endl;
      else{ cerr<<"File "<<file<<" did not contain valid data"<<endl; exit(1); }
      inFile.close();
    }
    else{ cerr<<"Could not open file "<<file<<endl; exit(1); }

    ave=omav; // for backwards compatibility
    rde=(mas>0?mas:1.f)/(doma*ave);
  }

  void add(float th, float ph){
    float ct=cos(th*fcv), st=sin(th*fcv);
    float cp=cos(ph*fcv), sp=sin(ph*fcv);
    V<3> dir;
    dir[0]=st*cp, dir[1]=st*sp, dir[2]=ct;
    dirs.push_back(dir);
  }

  float aS(float x){
    float al=acos(x);
    return al-sin(2*al)/2;
  }

  float f(float x){ // angular sensitivity curve (peaks at 1)
    if(beta<-1) return sqrt(1-x*x);
    else{
      float sum=x>0?x:0;
      if(beta<1){
	float y=sqrt(1-x*x)/beta;
	if(y>1){
	  y=1/y;
	  float c=1-beta*beta;
	  sum+=(aS(y)/c-aS(y*fabs(x)/sqrt(c))*fabs(x))/FPI;
	}
      }
      return sum;
    }
  }

  float xarea(float dot){ // OM cross-sectional area for direction with cos(zenith)=dot
    return Rz>0?FPI*Rr*sqrt(Rz*Rz-dot*dot*(Rz*Rz-Rr*Rr)):-4*Rr*Rz*sqrt(1-dot*dot);
  }

  void fraq(){ // max (over all directions) of PMT to sensor cross-sectional area
    if(def) return;

    float fr_flat=0, sum_ave=0;
    for(vector< V<3> >::iterator i=ico.dirs.begin(); i!=ico.dirs.end(); ++i){
      float tot=0;
      for(vector< V<3> >::const_iterator j=dirs.begin(); j!=dirs.end(); ++j) tot+=f(i->dot(*j));
      float fra=xarea((*i)[2]); fra=fra>0?tot/fra:0;
      if(fr_flat<fra) fr_flat=fra; sum_ave+=tot;
    }
    sum_ave/=ico.dirs.size();

    fx=fr_flat;
    rde=area*fr_flat/sum_ave;
  }

  int getPMT(V<3> dir, V<3> pos, V<3> tilt, float rnd, float ph = -1.f){
    if(def){
      bool flag;
      if(mas>0){
	float sum;
	{
	  float x = dir.dot(tilt);
	  float y=1;
	  sum=s[0];
	  for(unsigned int i=1; i<s.size(); i++){ y*=x; sum+=s[i]*y; }
	}

	flag=mas*rnd<sum;
      }
      else{
	flag=pos.dot(tilt)>s[0];
      }
      return flag?0:-1;
    }
    else{
      if(ph>=0){ // rotating photon by -ph instead of PMTs
	ph-=cable;
	float cp=cos(fcv*ph);
	float sp=sin(fcv*ph);
	float nx=dir[0]*cp+dir[1]*sp;
	float ny=dir[1]*cp-dir[0]*sp;
	dir[0]=nx, dir[1]=ny;
      }

      float crx=fx*xarea(dir[2]);

      unsigned int k=0;
      float tot=0;
      for(vector< V<3> >::const_iterator j=dirs.begin(); j!=dirs.end(); ++j, k++){
	tot+=f(-dir.dot(*j))/crx;
	if(tot>rnd) break;
      }

      return k<dirs.size()?(int) k:-1;
    }
  }
};

struct spec{
  bool swv; // single wavelength distribution
  int num;  // number of kept wavelength bins
  float bin, wmin, wmax;

  spec(): swv(false), num(450), wmin(250.f), wmax(700.f){
    bin=(wmax-wmin)/num;
  }

  void ssw(float wva){ // set single wavelength
    swv=true; num=10; wmin=wva-5.f; wmax=wva+5.f;
    bin=(wmax-wmin)/num;
  }

  float wav(float i){
    return wmin+bin*i;
  }

  float np(float wv, float * ng = NULL){ // phase and group refrative indices
    float np=1.55749-wv*(1.57988-wv*(3.99993-wv*(4.68271-2.09354*wv)));
    if(ng!=NULL) *ng=np*(1+0.227106-wv*(0.954648-wv*(1.42568-0.711832*wv)));
    return np;
  }

  float cherenkov(float wva){
    if(swv){ // over 10 nm
      return 1/(wmax-wmin);
    }
    else{ // um^-2 (here) * nm*cm^2 (later) = 0.1 m
      const float c=FPI*0.2/137.03599911;
      float wv=wva*1.e-3, n=np(wv);
      return c*(1-1/(n*n))/(wv*wv);
    }
  }
} qwv;

struct irde{
  float rde;
  vector<float> dat, sum, rat;

  int oms;  // number of OMs of this type
  float rmax; // max RDE for this type

  irde(): rde(0.f), oms(0), rmax(0.f){}

  void addr(float r){ // enter new OM
    if(r>rmax) rmax=r;
    oms++;
  }

  float binf(float p){
    int i, j=0; while(sum[++j]<=0); j--; while(sum[++j]<p); i=j-1;
    return i + (sum[j]>sum[i] ? (p-sum[i])/(sum[j]-sum[i]) : 0);
  }

  void read(const vector<float> & qw, const vector<float> & qf, float blo = -1.f, float bhi = -1.f){
    dat.clear();

    int k=0, n=qw.size();
    float wlo=qw[0]-(blo<0?qw[1]-qw[0]:blo);
    float whi=qw[n-1]+(bhi<0?qw[n-1]-qw[n-2]:bhi);
    float w1=wlo, w2=qw[0], f1=0, f2=qf[0];

    for(float w=qwv.wmin+qwv.bin/2; w<qwv.wmax; w+=qwv.bin){
      while(qw[k]<w && k<n){
	w1=w2, f1=f2; k++;
	if(k<n) w2=qw[k], f2=qf[k];
	else w2=whi, f2=0;
      }
      dat.push_back(wlo<w && w<whi ? (f1*(w2-w)+f2*(w-w1))/(w2-w1) : 0);
    }
  }
};

map<ikey, int> omts;
map<int, itype> types;
map<pair<int,int>, irde> irdes;
irde env;

map<ikey, float> hvs;
map<ikey, pair<float, int> > rdes;
map<ikey, V<3> > cx;
map<ikey, float > dx;

struct hit{
  unsigned int i, n, z;
  float t;
  float pth, pph, dth, dph;
};

#ifdef XCPU
struct int4{
  int x, y, z, w;
};

struct float2{
  float x, y;
};

struct float3:float2{
  float z;
};

struct float4:float3{
  float w;
};
#endif

struct pbuf{
  float4 r;        // location, time
  float4 n;        // direction
  unsigned int q;  // wavelength bin
  unsigned int i;  // segment number
  int fla, ofla;
};

struct photon{
  float4 r;    // location, time
  float4 n;    // direction, track length
  unsigned int q; // track segment
  unsigned int num; // number of photons in this bunch
  int type;    // source type
  float f;     // fraction of light from muon alone (without cascades)
  union{
    struct{
      float a, b;  // longitudinal development parametrization coefficients
      float beta;  // velocity of incident particle
      float tau;   // luminescence decay time
    };
    struct{
      float ka, up; // 2d-gaussian rms and zenith of cone
      float fldr;   // horizontal direction of the flasher led #1
      short fla, ofla;
    };
    int4 c;  // for quick copy
  };
};

struct ices{
  int wvl;               // wavelength number of this block
  float ocm;             // 1 / speed of light in medium
  float coschr, sinchr;  // cos and sin of the cherenkov angle
  struct{
    float abs;           // absorption
    float sca;           // scattering
  } z [MAXLYS];
};

struct aniz{
  float k1;
  float k2;
  float ra;
  float rb;
};

struct line{
  short n, max;
  float x, y, r;
  float h, d;
  float dl, dh;
};

struct datz{
  ices w[WNUM];
  float2 lp[LMAX][LYRS];
  unsigned int rm[MAXRND];
  unsigned long long rs[MAXRND];
} z;

struct dats{
  unsigned int hidx;

#ifndef XCPU
  unsigned int tn, tx;  // kernel time clocks
  unsigned int ab;      // if TOT was abnormal
  unsigned int mp;      // kernel block counter
  short bmp[4];         // list of 4 faulty MPs
#endif
  short blockIdx, gridDim;  // bad/current MP; number of MPs

  unsigned int gdev;  // number of this GPU
  unsigned int gnum;  // number of all GPUs
  unsigned int gini, gspc, gtot, gdiv;

  float rx;
  float hifl;

  unsigned int hnum;    // size of hits buffer
  int size;   // size of kurt table
  int rsize;  // count of multipliers
  int gsize;  // count of initialized OMs

  short tmod; // tilt model: 1: 1d, 2: 2d
  short vthk; // 0: uniform layers, 1: variable thickness
  float dh, rdh, hmin; // step, 1/step, and min depth

  float ocv;  // 1 / speed of light in vacuum
  float sf;   // scattering function: 0=HG; 1=SAM
  float g, gr; // g=<cos(scattering angle)> and gr=(1-g)/(1+g)

  float xR;   // DOM oversize scaling factor
  float SF, G, GR; // hole ice sf, g, gr
  float hr, hr2, hs, ha; // hole ice radius, radius^2, effective scattering and absorption coefficients

  float azx, azy;  // ice anisotropy direction

  int cn[2];
  float cl[2], crst[2];

  float cb[2][2];

  int lnum, lpts, l0;
  float lmin, lrdz, r0;
  float lnx, lny;
  float lr[LMAX];

  float mmin[2], mstp[2];
  int mnum[2], mcut[2];

  float k1, k2, kz, fr; // ice absorption anisotropy parameters
  aniz az[MAXLYS];

  unsigned short ls[NSTR];
  unsigned short is[CX][CY];
  char mcol[CTX][CTY];

  float sum, bfr[12];

  line sc[NSTR];

  datz * z;
  hit * hits;
  photon * pz;
  pbuf * bf;
  DOM * oms;
} d;

struct doms{
  DOM oms[MAXGEO];
  name names[MAXGEO];
  struct{
    float w, i, f;

    float x(){
      return i+(f-i)*xrnd();
    }
  } wvs [WNUM];

  hit * hits;
  photon * pz;
  pbuf * bf;

  float eff;  // OM efficiency correction
} q;

unsigned short sname(int n){
  static map<int, unsigned short> overflow;
  name & s = q.names[n];
  if(s.str>86){
    static unsigned short next_num=86+10+1;
    map<int, unsigned short>::iterator it=overflow.find(s.str);
    if(it==overflow.end()){
      if(next_num>=0x8000){ cerr<<"Number of string labels exceeds capacity of "<<(unsigned int) 0x8000<<endl; exit(1); }
      overflow.insert(make_pair(s.str, next_num++));
    }
    return overflow[s.str];
  }
  else return (unsigned short)( (s.str>78 && s.dom>10) ? s.str+10 : s.str );
}

static const float zoff=1948.07;
unsigned int sv=0;

void rs_ini(){
  union{
    unsigned long long da;
    struct{
      unsigned int ax;
      unsigned int dx;
    };
  } s;

  s.ax=362436069;
  s.dx=1234567;

  s.ax+=sv;

  for(int i=0; i<d.rsize; i++) z.rs[i]=s.da;
}


struct ini {
  void set() {
  //  float d = 0.0f;
  

    string ppcdir("");
    {
      char * env = getenv("PPCTABLESDIR"); // gets PPCTABLESDIR location
      if(env!=NULL) {
        ppcdir=string(env)+"/"; // sets PPCTABLESDIR str
        cerr<<"Configuring ppc in \""<<ppcdir<<"\""<<endl;
      }
      else {
        cerr << "NO PPCTABLESDIR SET!" << endl;
      }
    }

    string omdir("");
    {
      char * env = getenv("NEXTGENDIR");
      if (env!=NULL) {
        omdir=string(env) + "/";
        cerr << "using user set NEXTGENDIR: " << omdir << endl;
      } else {
        omdir = ppcdir;
        cerr << "using PPCTABLESDIR as omdir: " << omdir << endl;
      }
    }

    string icedir("");
    {
      char * env = getenv("ICEMODELDIR");
      if(env!=NULL) icedir=string(env)+"/";
      else icedir=ppcdir;
      cerr<<"Configuring icemodel in \""<<icedir<<"\""<<endl;
    }

    string tiltdir("");
    {
      char * env = getenv("TILTMODELDIR");
      if(env!=NULL) tiltdir=string(env)+"/";
      else tiltdir=icedir;
      cerr<<"Configuring tiltmodel in \""<<tiltdir<<"\""<<endl;
    }

    string holeice("");
    {
      char * env = getenv("PPCHOLEICE");
      if(env!=NULL) holeice=string(env);
      else holeice=ppcdir+"as.dat";
      cerr<<"Configuring holeice from \""<<holeice<<"\""<<endl;
    }

    float dk1, dk2, dkz; // ice anisotropy parameters

    {
      ifstream inFile((icedir+"cfg.txt").c_str(), ifstream::in);
      if(!inFile.fail()){
        string in;
        float aux;
        vector<float> v;

        while (getline(inFile, in)) {
          if(sscanf(in.c_str(), "%f", &aux)==1) v.push_back(aux);
        }
      

        {
          char * OVSZ=getenv("OVSZ");

          if(OVSZ!=NULL){
            float ovsz=atof(OVSZ); if(v.size()>=1) v[0]=ovsz;
            cerr<<"Using oversize as set with OVSZ="<<ovsz<<endl;
          }
        }

        if(v.size()>=4) { // at least 4 lines in cfg.txt
          int xR = lroundf(v[0]); // first line is oversize radius
          d.xR = xR;
          ovr*=xR*xR; // squared because 2D pancake

          q.eff = v[1]; // second line is DOM effeciency correction
          d.sf  = v[2]; // third line is scattering function: 0=HG; 1=SAM
          d.g   = v[3]; // fourth line is g=<cos(scattering angle)>
          d.gr = (1-d.g)/(1+d.g);

          cerr << "Configured xR (oversize) = " << xR << " eff= " << q.eff << " sf scattering= " << d.sf << " g=<cos(scattering angle)>: " << d.g << endl;
          
          if (v.size()<12) {
            d.SF = d.sf; // hole ice properties
            d.G = d.g;
            d.GR = d.gr;
          } else {
            d.SF=v[10], d.G=v[11], d.GR=(1-d.G)/(1+d.G);
          }

          if (v.size()>=10) {
            float xH=v[7]; // hole ice radius in units of [DOM radius]
            float hS = v[8]; // hole ice effective scattering length [m]
            float hA = v[9]; //hole ice absorption length [m]
            d.hr = OMR*xH; // hole ice radius in meters;
            d.hr2 = d.hr*d.hr; 
            d.hs = 1/(hS*(1-d.G));
            d.ha = 1/hA;
            if (xH > 0) cerr << "Simulating hole ice with DOM radius: " << xH << " , and sca length: " << hS << "m " <<" ("<<d.SF<<","<<d.G<<") abs="<<hA<<endl;
          } else {
            d.hr=0, d.hr2=0, d.hs=0, d.ha=0;
            cerr << "NOT simulating hole ice" << endl;
          }

          if(v.size()>=7){
            const float thx=v[4];
            d.azx=cos(fcv*thx), d.azy=sin(fcv*thx);
            dk1=exp(v[5]); dk2=exp(v[6]); dkz=1/(dk1*dk2);
            cerr<<"Ice anisotropy is k("<<thx<<")="<<dk1<<","<<dk2<<","<<dkz<<endl;
          }
          else {
            dk1=1, dk2=1, dkz=1, d.azx=1, d.azy=0;
            cerr << "NOT simulating anisotropy" << endl;
          }

          if(v.size()>=15){
            // new absorption anisotropy
            d.k1=exp(v[12]); d.k2=exp(v[13]); d.kz=exp(v[14]);
            cerr<<"New Ice anisotropy is "<<d.k1<<","<<d.k2<<","<<d.kz<<endl;
            if(d.k1>=d.k2 && d.k2==d.kz){
              float r=d.k1/d.k2;
              float s=sqrt(r*r-1);
              r=(s>XXX?log(r+s)/s:1+s/2)/d.k2;
              d.k1*=r, d.k2*=r, d.kz*=r;
              cerr<<"Renorm. NI anisotropy "<<d.k1<<","<<d.k2<<","<<d.kz<<endl;
            }
            else{
              cerr<<"Warning: taking absorption anisotropy as given (not renormalizing)!"<<endl;
            }

          }
          else d.k1=1, d.k2=1, d.kz=1;

          if(v.size()>=16){
            // scaling for absorption anisotropy (old implementation)
            d.fr=v[15];
            cerr<<"Ice absorption anisotropy scaling is "<<d.fr<<endl;
          }
          else d.fr=1;


          if(v.size()>=28){
            for(int i=0; i<12; i++) d.bfr[i]=v[16+i];

            {
              char * BFRA=getenv("BFRA");
              float bfra=BFRA==NULL?1.0:atof(BFRA);

              char * BFRB=getenv("BFRB");
              float bfrb=BFRB==NULL?1.0:atof(BFRB);

              if(BFRA!=NULL || BFRB!=NULL) cerr<<"Setting BFRA="<<bfra<<" BFRB="<<bfrb<<endl;
              for(int i=0; i<12; i+=4) d.bfr[i]*=i<8?sqrt(bfra):bfra*bfrb;
            }

            {
              float step=0.01, sum=0;
              for(float x=step/2; x<1; x+=step){
                float y=sqrt(1-x*x);
                float sx=max(0.f, d.bfr[0]*exp(-d.bfr[1]*pow(atan(d.bfr[3]*y), d.bfr[2])));
                float sy=max(0.f, d.bfr[4]*exp(-d.bfr[5]*pow(atan(d.bfr[7]*y), d.bfr[6])));
                float mx=max(0.f, d.bfr[8]*atan(d.bfr[11]*y*x)*exp(-d.bfr[9]*y+d.bfr[10]*x));
                sum+=sx*sx+sy*sy+mx*mx;
              }
              sum*=step/2; d.sum=sum;
            }
            cerr<<"Initialized BFR diffusion patterns; s_eff="<<d.sum<<" m^-1"<<endl;
          }
          else{
            d.sum=0;
            for(int i=0; i<12; i++) d.bfr[i]=i%4<1?0:1;
          }
	      }
        else { cerr<<"File cfg.txt did not contain valid data"<<endl; exit(1); }
        inFile.close();
      }
      else { cerr<<"Could not open file cfg.txt"<<endl; exit(1); }
    }

    {
      ifstream inFile((omdir+"om.conf").c_str(), ifstream::in);
      if(!inFile.fail()){
        string in;
          while(getline(inFile, in)){
            int m; // module number type
            unsigned int n; // number of pmts
            itype t; // area thing (shoudl be 1 i think...)
            float th, ph; // theta and phi of pmt 
            float other; // cable degrees ?
            int read=sscanf(in.c_str(), "%*s %d %f %f %f %f %d %f %f %f", &m, &t.area, &t.beta, &t.Rr, &t.Rz, &n, &th, &ph, &other);
            t.cable=read>=9?other:0;
            if(read>=8){
              t.add(th, ph); 
              for(unsigned int i=1; i<n && getline(inFile, in); i++)
                if(2==sscanf(in.c_str(), "%f %f %f", &th, &ph, &other)) t.add(th, ph);
              if(t.dirs.size()==n) types.insert(make_pair(m, t));
            }
          }
        inFile.close();
      }

      if (!types.empty()) { //if no modules have bee ndeclraed
        if(ico.ini(omdir+"om.dirs")<1){ cerr<<"Error: could not initialize an array of directions"<<endl; exit(1); } // if length of dirs of om.dirs  is less than 1 return error
        nextgen=true; // using nextgendir
        for(map<int, itype>::iterator j=types.begin(); j!=types.end(); ++j) {
          j->second.fraq(); // calculate rde
          cerr<<" OM Type "<<j->first<<" with "<<j->second.dirs.size()<<" PMTs added ("<<j->second.rde<<")"<<endl;
	      }
      }
    }

    if(types.find(-1)==types.end()){ // read angular sensitivity parameters
    /// i think above -1 type is default so if default then use holeice
      types[-1].add(holeice);
    }

    { //initialize random numbers:
      int size;
      vector<unsigned int> rx; // rx is vector of random ints

      ifstream inFile((ppcdir+"rnd.txt").c_str(), ifstream::in);

      if(!inFile.fail()){
        string in;
        while(getline(inFile, in)){
          stringstream str(in);
          unsigned int a;
          if(str>>a) rx.push_back(a);
        }
        if(rx.size()<1){ cerr<<"File rnd.txt did not contain valid data"<<endl; exit(1); }
        inFile.close();
      }
      
      else{ cerr<<"Could not open file rnd.txt"<<endl; exit(1); }

      size=rx.size();
      if(size>MAXRND){
        cerr<<"Error: too many random multipliers ("<<size<<"), truncating to "<<MAXRND<<endl;
        size=MAXRND;
      }

      cerr<<"Loaded "<<size<<" random multipliers"<<endl;

      #ifndef DTMN
            timeval tv; gettimeofday(&tv, NULL);
            sv=1000000*(unsigned long long)tv.tv_sec+tv.tv_usec;
      #endif

      d.rsize=size;
      for(int i=0; i<size; i++) z.rm[i]=rx[i];
    }

    {
      ifstream inFile((ppcdir+"geo-f2k").c_str(), ifstream::in);
      if(!inFile.fail()){
        if(!i3oms.empty()){
          i3oms.clear();
          cerr<<"Warning: overwriting existing geometry!"<<endl;
        }
        OM om;
        string mbid;
        unsigned long long omid;
        while(inFile>>mbid>>hex>>omid>>dec>>om.r[0]>>om.r[1]>>om.r[2]>>om.str>>om.dom){
          // first two columns of geo-f2k are throwaway
          om.R=OMR, om.F=1;
          om.r[2]+=zoff;
          i3oms.push_back(om);
        }
        inFile.close();
      }
    }


    {
      map<ikey, pair<float, float> > geco;
      for(vector<OM>::iterator i=i3oms.begin(); i!=i3oms.end(); ++i){
        map<ikey, pair<float, float> >::iterator j=geco.find(*i);
        if(j!=geco.end()){
          pair<float, float>& gc = j->second;
          i->r[0]+=gc.first, i->r[1]+=gc.second;
        }
      }
    }

    // {
    //   ifstream inFile((ppcdir+"str-f2k").c_str(), ifstream::in);
    //   if(!inFile.fail()){
    //     cerr<<"Using fixed positions of hole ice columns from str-f2k!"<<endl;

    //     int str;
    //     float x, y;
    //     while(inFile>>str>>x>>y) strs[str]=make_pair(x, y);
    //     inFile.close();
    //   }
    
    // }

    {
      ifstream inFile((ppcdir+"eff-f2k").c_str(), ifstream::in);
      if(!inFile.fail()){
        if(!rdes.empty()){
          rdes.clear();
          cerr<<"Warning: overwriting existing RDE table!"<<endl;
        }
        int typ;
        ikey om;
        float eff;

        string in;
        while(getline(inFile, in)){
          int num=sscanf(in.c_str(), "%d %d %f %d", &om.str, &om.dom, &eff, &typ);
          if(num<4) typ=0;
          if(num>=3) if(om.isinice()) rdes.insert(make_pair(om, make_pair(eff, typ)));
        }
        inFile.close();
      }
    }

    {
      ifstream inFile((omdir+"om.map").c_str(), ifstream::in);
      if(!inFile.fail()){
        if(!omts.empty()){
          omts.clear();
          cerr<<"Warning: overwriting existing OM type table!"<<endl;
        }
        ikey om;
        int m;
        while(inFile>>om.str>>om.dom>>m) if(om.isinice()) omts.insert(make_pair(om, m));
        inFile.close();
      }
    }

    //omitting  hvs-f2k, cx.dat, dx.dat because i dont want to 

    float Rr=0, Rz=0;

    { // initialize geo
      vector<DOM> oms;
      vector<name> names;
      int nhqe = 0;

      sort(i3oms.begin(), i3oms.end());
      for(vector<OM>::iterator i=i3oms.begin(); i!=i3oms.end(); ++i) if(i->isinice()){
        ikey om(*i);

        int m, t;
        float r, h;
        V<3> tilt;
        float azi = -1.f;

        if(omts.empty()) m=-1;
        else{
          map<ikey, int>::iterator j=omts.find(om);
          m=j==omts.end()?-1:j->second;
          map<int, itype>::const_iterator it=types.find(m);
          if(it==types.end()) m=-1;

          it=types.find(m);
          if(it!=types.end()){
            const itype & t=it->second;
            i->R=t.Rr, i->F=t.Rz/t.Rr;
          }
        }

        Rr=max(Rr, i->R), Rz=max(Rz, i->R*fabs(i->F));


        oms.push_back(*i);

        {
          map<ikey, pair<float, int> >::iterator j=rdes.find(om);
          if(j!=rdes.end()){
            nhqe++;
            r=j->second.first;
            t=j->second.second;
          }
          else r=1, t=0;

          irdes[make_pair(m,t)].addr(r);
        }

        if(hvs.empty()) h=1200;
        else{
          map<ikey, float>::iterator j=hvs.find(om);
          h=j==hvs.end()?0:j->second;
        }


        {
          map<ikey, V<3> >::iterator ci=cx.find(om);
          if(ci!=cx.end()) tilt=ci->second;

          map<ikey, float >::iterator di=dx.find(om);
          if(di!=dx.end()) azi = di->second;
        }

	      names.push_back(name(om, m, t, r, h, tilt, azi));
      }

      if(nhqe>0) cerr<<"Loaded "<<nhqe<<" RDE coefficients"<<endl;
      for(map<pair<int,int>, irde>::const_iterator i=irdes.begin(); i!=irdes.end(); ++i)
	    cerr<<" Found "<<i->second.oms<<" OMs of type "<<i->first.first<<","<<i->first.second<<endl;

      int gsize = oms.size();
      if(gsize>MAXGEO){
        cerr<<"Error: too many OMs ("<<gsize<<"), truncating to "<<MAXGEO<<endl;
        gsize=MAXGEO;
      }

      for(int n=0; n<gsize; n++){ q.oms[n]=oms[n]; q.names[n]=names[n]; }

      d.gsize=gsize;

    }

    Rr*=d.xR, Rz*=d.xR;



  }
} m;