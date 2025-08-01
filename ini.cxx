static unsigned int ovr=1; // oversize radius, default 1



struct dats {
  float SF, G, GR; // hole ice sf, g, gr
  float xR; // DOM oversize scaling factor
} d;

struct ini {
  void set() {
    float d = 0.0f;
  

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
      } s

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

        if (v.size()<12) {
          d.SF = d.sf; // hole ice properties
          d.G = d.g;
          d.GR = d.gr;
        }
      }


    }


  } 
}m;