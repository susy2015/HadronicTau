#ifndef EFFICIENCY_H
#define EFFICIENCY_H

class Efficiency{
 public:
// The acceptance is from a measurement after baseline
// Stat. uncertainty of the number need to be provided as well
  static double reco(int ptbin, int actbin);
  static double iso(int ptbin, int actbin);
  static double SBaccMix(int searchbin);
  static double taumucorMix(int njetbin, int metbin);
  static double mtwMix(int njetbin, int metbin);
  static double isotrkeffMix_NjetNbjet(int njetbin, int nbjetbin);
  static double SBtaumucorMix(int searchbin);
  static double SBmtwMix(int searchbin);
  static double SBisotrkeffMix(int searchbin);
  static double mistag(int ptbin);  
  static double HTMHT_trgEff(double HT, int metbin);
  static double MuonTrkSF(int etabin);
  static double MuonIDSF(int ptbin, int etabin);

  static int Ptbin(double pt); 
  static int Ptbin1(double pt);
  static int PtbinIDSF(double pt);
  static int Actbin(double act);
  static int Njetbin(int njet);
  static int NJetbin(int nJet);
  static int NBjetbin(int nBjet);
  static int MetFinerbin(double met);
  static int metbin(double met);
  static int Trgmetbin(double met);
  static int etabin(double eta);
  static int etabinIDSF(double eta);

  static double accErr(int searchbin);
  static double mtwErr(int searchbin);
  static double mtwMetjecUp(int searchbin);
  static double mtwMetjecLow(int searchbin);
  static double mtwMetjerUp(int searchbin);
  static double mtwMetjerLow(int searchbin);
  static double taumucorErr(int searchbin);
  static double isotrkeffErrUp(int searchbin);
  static double isotrkeffErrLow(int searchbin);
  static double accpdfCentral(int searchbin);
  static double accpdfUp(int searchbin);
  static double accpdfDown(int searchbin);
  static double accscaleUp(int searchbin);
  static double accscaleDown(int searchbin);
  static double HTMHT_trgEffUp(double HT, int metbin);
  static double HTMHT_trgEffDown(double HT, int metbin);
};

int Efficiency::Ptbin(double pt){
  int bin =0;
  //if(pt>=10. && pt<20.) bin =1;
  if(pt>=20. && pt<30.) bin =1;  
  if(pt>=30. && pt<40.) bin =2;
  if(pt>=40. && pt<50.) bin =3;
  if(pt>=50. && pt<70.) bin =4;
  if(pt>=70. && pt<100.) bin =5;
  if(pt>=100.) bin =6;
  return bin;
 }

int Efficiency::Ptbin1(double pt){
  int bin =0;
  if(pt>=30. && pt<60.) bin =1;
  if(pt>=60. && pt<90.) bin =2;
  if(pt>=90. && pt<120.) bin =3;
  if(pt>=120. && pt<150.) bin =4;
  if(pt>=150. && pt<180.) bin =5;
  if(pt>=180. && pt<210.) bin =6;
  if(pt>=210. && pt<240.) bin =7;
  if(pt>=240. && pt<270.) bin =8;
  if(pt>=270. && pt<300.) bin =9;
  if(pt>=300. && pt<330.) bin =10;
  if(pt>=330. && pt<360.) bin =11;
  if(pt>=360. && pt<390.) bin =12;
  if(pt>=390. && pt<420.) bin =13;
  if(pt>=420. && pt<450.) bin =14;
  if(pt>=450. && pt<500.) bin =15;
  if(pt>=500. && pt<1000.) bin =16;
  if(pt>=1000.) bin =17;
  return bin;
}

int Efficiency::PtbinIDSF(double pt){
  int bin =0;
  if(pt>=20. && pt<25.) bin =1;
  if(pt>=25. && pt<30.) bin =2;
  if(pt>=30. && pt<40.) bin =3;
  if(pt>=40. && pt<50.) bin =4;
  if(pt>=50. && pt<60.) bin =5;
  if(pt>=60.) bin =6;
  return bin;
}

int Efficiency::Actbin(double act){
  int bin =0;
  if(act>=0.02 && act<0.05) bin =1;
  if(act>=0.05 && act<0.15) bin =2;
  if(act>=0.15 && act<1) bin =3;
  if(act>=1) bin =4;
  return bin;
}

int Efficiency::Njetbin(int njet){
  int bin = 0;
  if(njet==5) bin =1;
  if(njet==6) bin =2;
  if(njet==7) bin =3;
  if(njet==8) bin =4;
  if(njet>=9) bin =5;
  return bin;
}

int Efficiency::NJetbin(int nJet){
  int bin = 0;
  if(nJet==5) bin =1;
  if(nJet==6) bin =2;
  if(nJet==7) bin =3;
  if(nJet==8) bin =4;
  if(nJet==9) bin =5;
  if(nJet>=10) bin =6;
  return bin;
}

int Efficiency::NBjetbin(int nBjet){
  int bin = 0;
  if(nBjet==2) bin =1;
  if(nBjet>=3) bin =2;
  return bin;
}

int Efficiency::metbin(double met){
  int bin = 0;
  if(met>=350 && met<500)bin =1;
  if(met>=500 && met<650)bin =2;
  if(met>=650)bin =3;
  return bin;
}

int Efficiency::MetFinerbin(double met){
  int bin = 0;
  if(met>=250 && met<275)bin =1;
  if(met>=275 && met<300)bin =2;
  if(met>=300 && met<350)bin =3;
  if(met>=350 && met<400)bin =4;
  if(met>=400 && met<450)bin =5;
  if(met>=450)bin =6;
}

int Efficiency::Trgmetbin(double met){
  int bin = 0;
  if(met>=275 && met<400)bin =1;
  if(met>=400)bin =2;
  return bin;
}
int Efficiency::etabin(double eta){
  int bin =0;
  if(eta>=-2.1 && eta<-1.6) bin =1;
  if(eta>=-1.6 && eta<-1.1) bin =2;
  if(eta>=-1.1 && eta<-0.6) bin =3;
  if(eta>=-0.6 && eta<0.0) bin =4;
  if(eta>=0.0 && eta<0.6) bin =5;
  if(eta>=0.6 && eta<1.1) bin =6;
  if(eta>=1.1 && eta<1.6) bin =7;
  if(eta>=1.6 && eta<2.1) bin =8;
  if(eta>=2.1 && eta<2.4) bin =9;
  return bin;
}

int Efficiency::etabinIDSF(double eta){
  int bin =0;
  if(eta>=0.9 && eta<1.2) bin =1;
  if(eta>=0.9 && eta<2.1) bin =2;
  if(eta>=2.1)            bin =3;
  return bin;
}

double Efficiency::reco(int ptbin, int actbin){
  double mu_recoeff[7][5] = {{0.97661,0.972552,0.963392,0.954908,0.960335},{0.982219,0.990507,0.9824,0.978234,0.980237},{0.989838,0.989538,0.990029,0.983527,0.98285},{0.984342,0.98666,0.989442,0.992909,0.964613},{0.992758,0.990257,0.991292,0.986346,0.96769},{0.988949,0.987315,0.985676,0.963952,0.958604},{0.961157,0.974281,0.971273,0.941513,0.981896}};

 double reco = mu_recoeff[ptbin][actbin];
 return reco;
}

double Efficiency::iso(int ptbin, int actbin){
  double mu_isoeff[7][5] ={{0.976009,0.96345,0.938758,0.830007,0.480129},{0.976622,0.974385,0.955499,0.837204,0.551807},{0.991752,0.988272,0.968003,0.840972,0.607771},{0.989423,0.990211,0.97113,0.832533,0.684722},{0.992414,0.988981,0.983772,0.892006,0.774208},{0.997727,0.996613,0.991423,0.953812,0.887645},{0.99716,0.997529,0.999378,0.981586,0.942588}};

  double iso = mu_isoeff[ptbin][actbin];
  return iso;
 }

double Efficiency::SBaccMix(int searchbin){
  double mu_acc[59]={0.784, 0.822, 0.854, 0.978, 0.820, 0.786, 0.796, 0.848, 0.944, 0.824, 0.791, 0.920, 0.783, 0.797, 0.813, 0.872, 0.804, 0.812, 0.733, 0.956, 0.796, 0.718, 0.985, 0.854, 0.897, 0.902, 0.837, 0.772, 0.871, 0.947, 0.964, 0.967, 0.797, 0.859, 0.961, 0.859, 0.887, 0.816, 0.916, 0.818, 0.889, 0.942, 0.923, 0.815, 0.865, 0.924, 0.875, 0.827, 0.910, 0.900, 0.959, 0.840, 0.870, 0.882, 0.963, 0.913, 0.928, 0.903, 0.950};
  //double mu_acc[59]={0.789, 0.816, 0.835, 0.914, 0.820, 0.794, 0.812, 0.859, 0.925, 0.834, 0.818, 0.874, 0.784, 0.794, 0.804, 0.864, 0.820, 0.802, 0.759, 0.957, 0.815, 0.916, 0.976, 0.809, 0.854, 0.845, 0.795, 0.774, 0.876, 0.946, 0.969, 0.939, 0.807, 0.868, 0.954, 0.759, 0.883, 0.834, 0.911, 0.868, 0.889, 0.941, 0.926, 0.822, 0.856, 0.920, 0.867, 0.798, 0.903, 0.883, 0.940, 0.794, 0.831, 0.888, 0.963, 0.916, 0.933, 0.894, 0.942};
   double acc = mu_acc[searchbin];
  return acc;
}
double Efficiency::mistag(int ptbin){
  double mistag[18] = {0,  0.126438,  0.153027,  0.164752,  0.169983,  0.155869,  0.138723,  0.125972,  0.120945,  0.113134,  0.103171,  0.0996899,  0.0976882,  0.100576,  0.0939329,  0.0790066,  0.073263,  0.073263};//CSV:0.8

  double rate = mistag[ptbin];
  return rate;
}

double Efficiency::mtwMix(int njetbin, int metbin){
  double mtwcorrMix[6][4] = {{0.892074, 0.905599, 0.919644, 0.896667},{0.907548, 0.944439, 0.95243, 0.92549},{0.910834, 0.943973, 0.949883, 0.776983},{0.907159, 0.926849, 0.88718, 0.91926},{0.9006, 0.933275, 0.876988, 0.912926},{0.886225, 0.922286, 0.925417, 0.875297}};
double mtwMix = mtwcorrMix[njetbin][metbin];
  return mtwMix;
}
double Efficiency::SBmtwMix(int searchbin){
  double mtwcorrMix[59] = {0.895, 0.922, 0.938, 0.897, 0.918, 0.934, 0.945, 0.944, 0.883, 0.926, 0.888, 0.892, 0.908, 0.937, 0.842, 0.888, 0.921, 0.915, 0.932, 0.924, 0.931, 0.917, 0.814, 0.907, 0.907, 0.902, 0.922, 0.939, 0.899, 0.922, 0.956, 0.916, 0.913, 0.923, 0.912, 0.972, 0.840, 0.926, 0.956, 0.936, 0.906, 0.930, 0.946, 0.900, 0.935, 0.933, 0.825, 0.933, 0.923, 0.893, 0.921, 0.893, 0.945, 0.887, 0.938, 0.877, 0.831, 0.860, 0.943};
  //double mtwcorrMix[59] = {0.895, 0.922, 0.948, 0.900, 0.919, 0.935, 0.941, 0.946, 0.870, 0.925, 0.893, 0.893, 0.908, 0.935, 0.827, 0.889, 0.921, 0.917, 0.912, 0.946, 0.931, 0.923, 0.806, 0.911, 0.909, 0.928, 0.924, 0.945, 0.900, 0.924, 0.961, 0.903, 0.916, 0.923, 0.907, 0.972, 0.843, 0.933, 0.959, 0.945, 0.908, 0.933, 0.944, 0.898, 0.937, 0.932, 0.822, 0.938, 0.934, 0.895, 0.918, 0.903, 0.954, 0.886, 0.939, 0.873, 0.832, 0.854, 0.942};
double mtwMix = mtwcorrMix[searchbin];
  return mtwMix;
}

double Efficiency::isotrkeffMix_NjetNbjet(int njetbin, int nbjetbin){
  double isotrk[6][3] = {{0.300, 0.151, 0.150},{0.311, 0.285, 0.143},{0.375, 0.322, 0.277},{0.392, 0.332, 0.304},{0.362, 0.354, 0.323},{0.358, 0.355, 0.359}};//mix (ttbar singlelep, Wjets, single top)
  double eff = isotrk[njetbin][nbjetbin];
  return eff;
}
double Efficiency::SBisotrkeffMix(int searchbin){
  double isotrk[59] = {0.340, 0.348, 0.267, 0.314, 0.412, 0.456, 0.476, 0.222, 0.485, 0.422, 0.528, 0.287, 0.305, 0.301, 0.311, 0.168, 0.362, 0.408, 0.323, 0.269, 0.417, 0.512, 0.760, 0.259, 0.391, 0.345, 0.318, 0.159, 0.394, 0.414, 0.332, 0.541, 0.383, 0.397, 0.321, 0.415, 0.534, 0.356, 0.304, 0.282, 0.377, 0.429, 0.289, 0.331, 0.347, 0.373, 0.457, 0.321, 0.619, 0.311, 0.385, 0.341, 0.263, 0.405, 0.253, 0.294, 0.331, 0.392, 0.677};
  //double isotrk[59] = {0.343, 0.349, 0.260, 0.316, 0.415, 0.460, 0.468, 0.219, 0.460, 0.423, 0.529, 0.287, 0.305, 0.301, 0.301, 0.165, 0.358, 0.397, 0.307, 0.254, 0.412, 0.531, 0.760, 0.261, 0.392, 0.346, 0.308, 0.164, 0.392, 0.418, 0.354, 0.451, 0.389, 0.395, 0.325, 0.397, 0.558, 0.343, 0.307, 0.280, 0.373, 0.436, 0.293, 0.337, 0.344, 0.367, 0.477, 0.297, 0.623, 0.306, 0.388, 0.345, 0.263, 0.407, 0.271, 0.296, 0.342, 0.404, 0.718};
  double eff = isotrk[searchbin];
  return eff;
}

double Efficiency::taumucorMix(int njetbin, int metbin){
  double SBtaumucorr[6][4]= {{0.177389, 0.219184, 0.213738, 0.228602},{0.152886, 0.210099, 0.293804, 0.249091},{0.148117, 0.185513, 0.227056, 0.231356},{0.14619, 0.189357, 0.215143, 0.283387},{0.147645, 0.17074, 0.145709, 0.397917},{0.139241, 0.178362, 0.161042, 0.153092}};
  double SBtaumucor = SBtaumucorr[njetbin][metbin];
  return SBtaumucor;
}
double Efficiency::SBtaumucorMix(int searchbin){
  double SBtaumucorr[59]= {0.129, 0.161, 0.160, 0.200, 0.126, 0.157, 0.157, 0.186, 0.134, 0.135, 0.224, 0.172, 0.131, 0.167, 0.168, 0.135, 0.136, 0.173, 0.178, 0.250, 0.169, 0.242, 0.196, 0.129, 0.152, 0.203, 0.145, 0.215, 0.130, 0.161, 0.189, 0.157, 0.132, 0.148, 0.147, 0.152, 0.135, 0.168, 0.272, 0.206, 0.120, 0.145, 0.156, 0.145, 0.162, 0.151, 0.108, 0.140, 0.173, 0.134, 0.121, 0.137, 0.153, 0.119, 0.133, 0.130, 0.136, 0.096, 0.148};
  //double SBtaumucorr[59]= {0.128, 0.162, 0.157, 0.237, 0.124, 0.157, 0.164, 0.139, 0.127, 0.132, 0.219, 0.172, 0.130, 0.173, 0.170, 0.123, 0.128, 0.170, 0.192, 0.544, 0.151, 0.290, 0.185, 0.142, 0.171, 0.188, 0.162, 0.236, 0.130, 0.159, 0.192, 0.179, 0.133, 0.153, 0.134, 0.150, 0.130, 0.170, 0.273, 0.205, 0.119, 0.144, 0.163, 0.150, 0.175, 0.153, 0.094, 0.134, 0.162, 0.141, 0.134, 0.149, 0.164, 0.125, 0.126, 0.124, 0.137, 0.092, 0.134};
  double SBtaumucor = SBtaumucorr[searchbin];
  return SBtaumucor;
}

double Efficiency::HTMHT_trgEff(double HT, int metbin){
  double trig_low[3] = {0.974, 0.993, 1.000};
  double trig_up[3] = {0.910, 0.989, 1.000};
  double trgEff = HT>1000?trig_up[metbin]:trig_low[metbin];
  return trgEff;
}

double Efficiency::accErr(int searchbin){
  double mu_accMix[59] ={0.004, 0.010, 0.028, 0.014, 0.010, 0.019, 0.051, 0.087, 0.032, 0.022, 0.036, 0.036, 0.004, 0.010, 0.040, 0.062, 0.018, 0.020, 0.053, 0.042, 0.034, 0.155, 0.010, 0.006, 0.016, 0.045, 0.023, 0.049, 0.004, 0.009, 0.024, 0.022, 0.018, 0.022, 0.033, 0.111, 0.044, 0.037, 0.043, 0.075, 0.004, 0.010, 0.034, 0.020, 0.024, 0.043, 0.057, 0.037, 0.044, 0.008, 0.018, 0.025, 0.033, 0.021, 0.030, 0.019, 0.040, 0.033, 0.049};
 double err = mu_accMix[searchbin];
 return err;
}
double Efficiency::accpdfUp(int searchbin){
  double mu_accMix[59] ={0.783, 0.821, 0.851, 0.977, 0.818, 0.784, 0.796, 0.855, 0.942, 0.823, 0.785, 0.920, 0.782, 0.796, 0.815, 0.877, 0.804, 0.811, 0.737, 0.954, 0.813, 0.721, 0.984, 0.853, 0.899, 0.902, 0.836, 0.773, 0.870, 0.946, 0.966, 1.000, 0.794, 0.857, 0.963, 0.864, 0.882, 0.824, 0.914, 0.815, 0.889, 0.941, 0.916, 0.812, 0.862, 0.925, 0.875, 0.829, 0.912, 0.900, 0.969, 0.837, 0.866, 0.894, 0.962, 0.911, 0.929, 0.903, 0.950};
  double err = mu_accMix[searchbin];
  return err;
}
double Efficiency::accpdfDown(int searchbin){
  double mu_accMix[59] ={0.786, 0.823, 0.858, 0.979, 0.823, 0.788, 0.800, 0.846, 0.944, 0.825, 0.796, 0.921, 0.785, 0.798, 0.808, 0.882, 0.806, 0.816, 0.727, 0.959, 0.812, 0.714, 0.987, 0.856, 0.902, 0.902, 0.836, 0.770, 0.872, 0.947, 0.967, 1.000, 0.799, 0.862, 0.962, 0.838, 0.892, 0.836, 0.918, 0.821, 0.890, 0.943, 0.931, 0.817, 0.867, 0.921, 0.874, 0.826, 0.905, 0.900, 0.967, 0.844, 0.874, 0.892, 0.964, 0.916, 0.927, 0.902, 0.949};
  double err = mu_accMix[searchbin];
  return err;
}
double Efficiency::accscaleUp(int searchbin){
  double mu_accMix[59] ={0.781, 0.811, 0.828, 0.979, 0.812, 0.789, 0.769, 0.831, 0.937, 0.810, 0.781, 0.914, 0.781, 0.793, 0.835, 0.857, 0.799, 0.782, 0.721, 0.944, 0.790, 0.720, 0.963, 0.854, 0.894, 0.881, 0.832, 0.803, 0.870, 0.945, 0.963, 1.000, 0.796, 0.851, 0.957, 0.861, 0.884, 0.809, 0.904, 0.809, 0.886, 0.939, 0.924, 0.820, 0.858, 0.925, 0.874, 0.825, 0.891, 0.898, 0.959, 0.836, 0.860, 0.882, 0.965, 0.911, 0.917, 0.900, 0.952};
   double err = mu_accMix[searchbin];
  return err;
}
double Efficiency::accscaleDown(int searchbin){
  double mu_accMix[59] ={0.784, 0.813, 0.830, 0.977, 0.814, 0.795, 0.759, 0.823, 0.934, 0.806, 0.785, 0.905, 0.782, 0.792, 0.832, 0.857, 0.798, 0.785, 0.724, 0.947, 0.787, 0.716, 0.960, 0.856, 0.895, 0.891, 0.830, 0.800, 0.871, 0.947, 0.963, 1.000, 0.798, 0.854, 0.958, 0.849, 0.890, 0.816, 0.901, 0.814, 0.889, 0.942, 0.922, 0.824, 0.863, 0.919, 0.875, 0.827, 0.887, 0.901, 0.955, 0.841, 0.868, 0.881, 0.961, 0.915, 0.937, 0.904, 0.948};
  double err = mu_accMix[searchbin];
  return err;
}
double Efficiency::mtwErr(int searchbin){
  double mtwcorrMix[59] = {0.003, 0.006, 0.009, 0.025, 0.006, 0.007, 0.009, 0.009, 0.025, 0.020, 0.042, 0.022, 0.003, 0.004, 0.058, 0.043, 0.010, 0.020, 0.025, 0.030, 0.017, 0.027, 0.128, 0.003, 0.010, 0.035, 0.013, 0.012, 0.003, 0.012, 0.010, 0.043, 0.008, 0.011, 0.024, 0.011, 0.033, 0.018, 0.009, 0.018, 0.003, 0.007, 0.014, 0.009, 0.009, 0.020, 0.042, 0.018, 0.028, 0.006, 0.017, 0.016, 0.012, 0.016, 0.017, 0.016, 0.087, 0.034, 0.032};
 
double err = mtwcorrMix[searchbin];
  return err;
}
double Efficiency::mtwMetjecUp(int searchbin){
  double mtwcorrMix[59] = {0.891, 0.920, 0.936, 0.895, 0.913, 0.933, 0.945, 0.943, 0.860, 0.922, 0.886, 0.892, 0.904, 0.933, 0.841, 0.887, 0.920, 0.915, 0.932, 0.922, 0.930, 0.917, 0.813, 0.903, 0.905, 0.902, 0.921, 0.939, 0.895, 0.921, 0.955, 0.916, 0.910, 0.923, 0.911, 0.972, 0.839, 0.925, 0.954, 0.936, 0.902, 0.928, 0.943, 0.893, 0.935, 0.933, 0.825, 0.933, 0.923, 0.888, 0.921, 0.887, 0.945, 0.885, 0.938, 0.870, 0.830, 0.844, 0.929};
double err = mtwcorrMix[searchbin];
  return err;
}
double Efficiency::mtwMetjecLow(int searchbin){
  double mtwcorrMix[59] = {0.900, 0.925, 0.939, 0.898, 0.920, 0.936, 0.946, 0.945, 0.884, 0.928, 0.895, 0.895, 0.911, 0.940, 0.843, 0.889, 0.923, 0.918, 0.935, 0.949, 0.931, 0.918, 0.814, 0.911, 0.909, 0.903, 0.923, 0.939, 0.904, 0.925, 0.961, 0.918, 0.915, 0.925, 0.917, 0.973, 0.844, 0.927, 0.958, 0.937, 0.910, 0.933, 0.950, 0.905, 0.935, 0.934, 0.835, 0.934, 0.925, 0.896, 0.922, 0.905, 0.946, 0.888, 0.939, 0.880, 0.831, 0.861, 0.946};
  double err = mtwcorrMix[searchbin];
  return err;
}
double Efficiency::mtwMetjerUp(int searchbin){
  double mtwcorrMix[59] = {0.895, 0.922, 0.938, 0.897, 0.918, 0.934, 0.945, 0.944, 0.883, 0.926, 0.888, 0.892, 0.908, 0.937, 0.842, 0.888, 0.921, 0.915, 0.932, 0.924, 0.931, 0.917, 0.814, 0.907, 0.907, 0.902, 0.922, 0.939, 0.899, 0.922, 0.956, 0.916, 0.913, 0.923, 0.912, 0.972, 0.840, 0.926, 0.956, 0.936, 0.906, 0.930, 0.946, 0.900, 0.935, 0.933, 0.825, 0.933, 0.923, 0.893, 0.921, 0.893, 0.945, 0.887, 0.938, 0.877, 0.831, 0.860, 0.943};
  double err = mtwcorrMix[searchbin];
  return err;
}
double Efficiency::mtwMetjerLow(int searchbin){
  double mtwcorrMix[59] = {0.895, 0.922, 0.938, 0.897, 0.918, 0.934, 0.945, 0.944, 0.883, 0.926, 0.888, 0.892, 0.908, 0.937, 0.842, 0.888, 0.921, 0.915, 0.932, 0.924, 0.931, 0.917, 0.814, 0.907, 0.907, 0.902, 0.922, 0.939, 0.899, 0.922, 0.956, 0.916, 0.913, 0.923, 0.912, 0.972, 0.840, 0.926, 0.956, 0.936, 0.906, 0.930, 0.946, 0.900, 0.935, 0.933, 0.825, 0.933, 0.923, 0.893, 0.921, 0.893, 0.945, 0.887, 0.938, 0.877, 0.831, 0.860, 0.943};
  double err = mtwcorrMix[searchbin];
  return err;
}
double Efficiency::taumucorErr(int searchbin){
  double SBtaumucorr[59]= {0.003, 0.007, 0.015, 0.024, 0.008, 0.014, 0.025, 0.026, 0.031, 0.015, 0.052, 0.020, 0.003, 0.010, 0.020, 0.032, 0.018, 0.025, 0.034, 0.139, 0.033, 0.086, 0.090, 0.005, 0.011, 0.037, 0.015, 0.038, 0.003, 0.013, 0.041, 0.045, 0.013, 0.015, 0.028, 0.049, 0.035, 0.024, 0.098, 0.045, 0.003, 0.010, 0.028, 0.011, 0.018, 0.039, 0.038, 0.027, 0.041, 0.007, 0.019, 0.017, 0.025, 0.015, 0.032, 0.016, 0.036, 0.021, 0.056};

  double err = SBtaumucorr[searchbin];
  return err;
}
double Efficiency::isotrkeffErrUp(int searchbin){
  double isotrk[59] = {0.006, 0.019, 0.049, 0.060, 0.017, 0.033, 0.102, 0.155, 0.084, 0.047, 0.058, 0.081, 0.006, 0.019, 0.046, 0.091, 0.034, 0.049, 0.075, 0.112, 0.087, 0.122, 0.128, 0.011, 0.042, 0.090, 0.049, 0.055, 0.008, 0.026, 0.108, 0.142, 0.025, 0.038, 0.089, 0.138, 0.183, 0.061, 0.072, 0.169, 0.009, 0.032, 0.082, 0.023, 0.047, 0.105, 0.092, 0.058, 0.124, 0.019, 0.059, 0.065, 0.100, 0.041, 0.079, 0.052, 0.081, 0.072, 0.145};
  double err = isotrk[searchbin];
  return err;
}
double Efficiency::isotrkeffErrLow(int searchbin){
  double isotrk[59] = {0.006, 0.019, 0.046, 0.056, 0.017, 0.033, 0.100, 0.115, 0.083, 0.046, 0.058, 0.073, 0.006, 0.019, 0.043, 0.070, 0.033, 0.048, 0.069, 0.095, 0.084, 0.123, 0.172, 0.011, 0.041, 0.083, 0.046, 0.047, 0.008, 0.025, 0.097, 0.147, 0.025, 0.037, 0.081, 0.130, 0.189, 0.058, 0.066, 0.135, 0.009, 0.031, 0.074, 0.023, 0.045, 0.098, 0.090, 0.055, 0.135, 0.019, 0.057, 0.061, 0.086, 0.041, 0.070, 0.048, 0.075, 0.069, 0.173};
  double err = isotrk[searchbin];
  return err;
}
double Efficiency::HTMHT_trgEffUp(double HT, int metbin){
  double trig_low[3] = {0.002, 0.002, 0.000};
  double trig_up[3] = {0.010, 0.004, 0.000};
  double trgEff = HT>1000?trig_up[metbin]:trig_low[metbin];
  return trgEff;
}
double Efficiency::HTMHT_trgEffDown(double HT, int metbin){
  double trig_low[3] = {0.003, 0.002, 0.004};
  double trig_up[3] = {0.011, 0.006, 0.005};
  double trgEff = HT>1000?trig_up[metbin]:trig_low[metbin];
  return trgEff;
}
double Efficiency::MuonTrkSF(int etabin){
  double trk_SF[10] = {98.79, 99.39, 99.70, 99.54, 99.37, 99.59, 99.76, 99.61, 99.30, 98.19};
  double trkSF = trk_SF[etabin]/100;
  return trkSF;
}
double Efficiency::MuonIDSF(int ptbin, int etabin){
  double mu_IDSF[7][4] = {{0.958256, 0.940449, 0.93575, 0.984418},{0.972896, 0.992422, 0.987885, 0.979471},{0.987123, 0.986222, 0.984449, 0.970577},{0.991571, 0.992328, 0.989676, 0.971951},{0.991701, 0.991498, 0.989048, 0.970364},{0.985412, 0.986677, 0.984837, 0.96496},{0.9868, 0.990375, 0.98724, 0.958166}};
 double IDSF = mu_IDSF[ptbin][etabin];
 return IDSF;
}

#endif
