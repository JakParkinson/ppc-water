#define XXX 1.e-5f
#define FPI 3.141592653589793f
#define OMR 0.16510f  // DOM radius [m]

static unsigned int ovr=1; // oversize radius, default 1


static const float doma=FPI*OMR*OMR; // DOM cross-sectional area [m^2]
static const float omav=0.335822;    // average DOM angular (lab) sensitivity
static const float dppm=2450.08;     // photons per meter for nominal IceCube DOM
static const float fcv=FPI/180.f;


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


struct dats {
  float ocv;  // 1 / speed of light in vacuum
  float sf;   // scattering function: 0=HG; 1=SAM
  float g, gr; // g=<cos(scattering angle)> and gr=(1-g)/(1+g)

  float xR;   // DOM oversize scaling factor
  float SF, G, GR; // hole ice sf, g, gr
  float hr, hr2, hs, ha; // hole ice radius, radius^2, effective scattering and absorption coefficients

  float azx, azy;  // ice anisotropy direction
  float k1, k2, kz, fr; // ice absorption anisotropy parameters

  float sum, bfr[12];

} d;

bool nextgen=false;

struct itype {
  float area, beta, rde, fx, Rr, Rz, cable;
  vector <V<3>> dirs;

  bool def;
  float ave; // average angular sensitivty
  float mas; //maximum angular sensitiyivty
  vector <float> s; // angular sensitivty coefficients

  itype(): def(false) { }

  void add(string file) {
    def=true;
    mas = 1;
    ave = 0;
    area = 1;
    beta = 0.33f; // 'standard DOM' defaults:
    Rr = OMR;
    Rz = OMR;
    cable  = 0;

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

  void add(float th, float ph){ // adds a pmt dir i think
    float ct=cos(th*fcv), st=sin(th*fcv);
    float cp=cos(ph*fcv), sp=sin(ph*fcv);
    V<3> dir;
    dir[0]=st*cp, dir[1]=st*sp, dir[2]=ct;
    dirs.push_back(dir);
  }

  float aS(float x) {
    float al = acos(x);
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

  float xarea(float dot) { // OM cross-sectional area for direction with cos(zenith)=dot
    if (Rz>0) { // ie 
      return FPI*Rr*sqrt(Rz*Rz-dot*dot*(Rz*Rz-Rr*Rr));
    } else {
      return -4*Rr*Rz*sqrt(1-dot*dot);
    }

  }


}

struct doms {
  float eff; // OM efficiency correction (set by cfg.txt)
} q;

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
    }


  } 
} m;