#ifndef EFFICIENCY_H
#define EFFICIENCY_H

class Efficiency{
 public:
// The acceptance is from a measurement after baseline
// Stat. uncertainty of the number need to be provided as well
  static double acc(int njetbin);
  static double reco(int ptbin, int actbin);
  static double iso(int ptbin, int actbin);
  static double taumucorTT(int searchbin);
  static double taumucorWJet(int searchbin);
  static double accTT(int njetbin, int mt2bin);
  static double accWJet(int njetbin, int mt2bin);
  static int Ptbin(double pt); 
  static int Ptbin1(double pt);
  static int Actbin(double act);
  static int Njetbin(int njet);
  static int NJetbin(int nJet);
  static int MT2bin(double mt2);
  static double mistag(int ptbin);  
  static double SBmtwTT(int searchbin);
  static double SBmtwWJet(int searchbin);
  static double SBmtwMix(int searchbin);
  static double SBisotrkeffTT(int searchbin);
  static double SBisotrkeffWJet(int searchbin);
  static double isotrkeffTT(int njetbin);
  static double MT2Bin_mtwcorr() {return 0.921;}//TTbarSingleLep
  // static double MT2Bin_mtwcorr() {return 0.876;}//WJet
   static double MT2Bin_isotrkcorr(){return  0.354;} //(ttbarSingleLep) 
   //static double MT2Bin_isotrkcorr(){return 0.281;} //(WJet)
   static double MT2Bin_taumucorr(){return 0.156;}//TTbarSingleLep
   //static double MT2Bin_taumucorr(){return 0.144;}//(WJet)

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

int Efficiency::Actbin(double act){
  int bin =0;
  if(act>=5. && act<10.) bin =1;
  if(act>=10. && act<20.) bin =2;
  if(act>=20. && act<40.) bin =3;
  if(act>=40. && act<60.) bin =4;
  if(act>=60. && act<80.) bin =5;
  if(act>=80. && act<100.) bin =6;
  if(act>=100.) bin =7;
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

int Efficiency::MT2bin(double mt2){
  int bin = 0;
  if(mt2>=150 && mt2<200)bin =1;
  if(mt2>=200 && mt2<300)bin =2;
  if(mt2>=300 && mt2<400)bin =3;
  if(mt2>=400)bin =4;
  return bin;
}

double Efficiency::reco(int ptbin, int actbin){
  double mu_recoeff[7][8] = {{0.961671,0.932999,0.94894,0.974745,0.961365,0.93704,0.978094,0.959886},{0.982242,0.991145,0.985873,0.988798,0.975173,0.979503,0.972802,0.979673},{0.989233,0.982045,0.993872,0.982023,0.990218,0.992962,0.97708,0.976282},{0.995599,0.998264,0.993941,0.993821,0.983841,0.979065,0.992041,0.983542},{0.99238,0.986337,0.983247,0.990857,0.98217,0.99713,0.986692,0.98711},{0.987953,0.986693,0.97548,0.986735,0.989243,0.989195,0.983622,0.987648},{0.980322,0.9926,0.961815,0.959966,0.978283,0.969026,0.991232,0.981346}};

 double reco = mu_recoeff[ptbin][actbin];
 return reco;
}


double Efficiency::iso(int ptbin, int actbin){
  double mu_isoeff[7][8] = {{0.957087,0.732571,0.748459,0.805492,0.785299,0.759865,0.745552,0.709808},{0.979091,0.916709,0.814723,0.880013,0.860471,0.800555,0.839904,0.779696},{0.97982,0.943115,0.835159,0.873117,0.892093,0.901176,0.871294,0.842513},{0.984215,0.968748,0.926695,0.938785,0.916061,0.910372,0.892436,0.867568},{0.993318,0.986006,0.980538,0.968903,0.975371,0.953154,0.948252,0.913556},{0.997035,0.99858,1,0.989644,0.987459,0.986717,0.986926,0.971902},{0.996723,0.997938,1,0.999195,0.994069,0.998615,0.995497,0.992187}};

  double iso = mu_isoeff[ptbin][actbin];
  return iso;
 }


double Efficiency::acc(int njetbin){
  double mu_acc[6] = {0.479, 0.714, 0.794, 0.837, 0.860, 0.884};//ttbarSingleLep
  double acc = mu_acc[njetbin];
  return acc;
}

double Efficiency::accTT(int njetbin, int mt2bin){
  //double mu_acc[6][4] = {{0.473, 0.488, 0.521, 0.586},{0.779, 0.652, 0.631, 0.677},{0.850, 0.760, 0.755, 0.797},{0.879, 0.825, 0.798, 0.821},{0.889, 0.857, 0.841, 0.831},{0.897, 0.891, 0.865, 0.853}};
  //double mu_accTT[6][5] = {{0.622579, 0.593932, 0.628354, 0.654941, 0.648248},{0.697259, 0.803373, 0.732478, 0.680996, 0.708742},{0.792691, 0.861641, 0.795518, 0.769224, 0.799289},{0.851117, 0.885878, 0.835696, 0.804419, 0.821783},{0.855674, 0.894845, 0.860378, 0.843171, 0.832868},{0.836107, 0.900845, 0.891478, 0.866341, 0.85302}};//with HT cut
    double mu_accTT[6][5] = {{0.640498, 0.605869, 0.639647, 0.661332, 0.655662},{0.702658, 0.805425, 0.734796, 0.687402, 0.716927},{0.795322, 0.862231, 0.796648, 0.77136, 0.803028},{0.852213, 0.886316, 0.836255, 0.805418, 0.824536},{0.856476, 0.894955, 0.860801, 0.844134, 0.834918},{0.835806, 0.901035, 0.891526, 0.866662, 0.853679}};//with HT cut & new temp.

  double acc = mu_accTT[njetbin][mt2bin];
  return acc;
}

double Efficiency::accWJet(int njetbin, int mt2bin){
  double mu_accWJet[6][5] = {{0.876704, 0.635278, 0.708017, 0.646517, 0.635638},{0.889719, 0.862303, 0.68848, 0.621839, 0.61824},{0.865448, 0.774457, 0.746186, 0.730277, 0.747207},{0.965311, 0.826688, 0.768086, 0.779976, 0.82561},{0.982316, 0.744423, 0.767649, 0.759841, 0.826374},{0.997319, 0.678293, 0.726884, 0.796241, 0.8659}};

 double acc = mu_accWJet[njetbin][mt2bin];
  return acc;
}

double Efficiency::mistag(int ptbin){
  //  double mistag[18] = {0., 0.0929165, 0.0977224, 0.097206, 0.0984691, 0.0988062, 0.0769717, 0.0807692, 0.075177, 0.0796671, 0.0665965, 0.0717385, 0.0602786, 0.0747461, 0.0429148, 0.0590428, 0.047948, 0.047948};
  double mistag[18] = {0, 0.0563114, 0.0654405, 0.0842772, 0.0854902, 0.0999826, 0.0955245, 0.0847876, 0.0753087, 0.0704465, 0.0702774, 0.056887, 0.0747066, 0.0765881, 0.0600503, 0.0456406, 0.0472811, 0.0472811};
  double rate = mistag[ptbin];
  return rate;
}


double Efficiency::SBmtwTT(int searchbin){
  //double SBmtwcorr[45] = {0.924, 0.939, 0.940, 0.940, 0.929, 0.946, 0.957, 0.918, 0.925, 0.946, 0.970, 0.927, 0.936, 0.936, 0.956, 0.924, 0.945, 0.950, 0.957, 0.947, 0.962, 0.903, 0.932, 0.913, 0.908, 0.920, 0.940, 0.906, 0.926, 0.935, 0.892, 0.950, 0.909, 0.927, 0.917, 0.914, 0.922, 0.913, 0.876, 0.935, 0.886, 0.920, 0.896, 0.893, 0.848};//TTbarSingleLep
  //double SBmtwcorr[45] = {0.899, 0.929, 0.936, 0.943, 0.912, 0.931, 0.949, 0.917, 0.921, 0.942, 0.969, 0.908, 0.928, 0.932, 0.955, 0.903, 0.936, 0.943, 0.956, 0.944, 0.962, 0.890, 0.928, 0.913, 0.904, 0.918, 0.940, 0.903, 0.921, 0.935, 0.889, 0.949, 0.904, 0.925, 0.916, 0.911, 0.920, 0.915, 0.875, 0.935, 0.884, 0.920, 0.895, 0.893, 0.848};//TTbarSingleLep with HT cut
  double SBmtwcorr[45] = {0.900, 0.930, 0.937, 0.943, 0.913, 0.931, 0.950, 0.916, 0.922, 0.943, 0.968, 0.909, 0.929, 0.933, 0.955, 0.905, 0.936, 0.943, 0.955, 0.943, 0.962, 0.891, 0.929, 0.914, 0.904, 0.918, 0.939, 0.904, 0.921, 0.933, 0.891, 0.948, 0.904, 0.926, 0.915, 0.912, 0.920, 0.915, 0.875, 0.934, 0.885, 0.920, 0.896, 0.895, 0.847};//TTbarSingleLep with HT cut & new temp

  double SBmtw = SBmtwcorr[searchbin];
  return SBmtw;
}

double Efficiency::SBmtwWJet(int searchbin){
  double SBmtwcorrWJet[45] = {0.871, 0.902, 0.903, 0.893, 0.868, 0.875, 0.884, 0.929, 0.863, 0.879, 0.888, 0.908, 0.893, 0.884, 0.916, 0.907, 0.898, 0.928, 0.879, 0.913, 0.890, 0.861, 0.840, 0.914, 0.870, 0.844, 0.944, 0.855, 0.864, 0.904, 0.829, 0.855, 0.863, 0.828, 0.796, 0.827, 0.916, 0.956, 0.885, 0.918, 0.930, 0.849, 0.879, 0.969, 0.972};
  double SBmtwWJet = SBmtwcorrWJet[searchbin];
  return SBmtwWJet;
}

double Efficiency::SBmtwMix(int searchbin){
  double SBmtwcorrMix[64] = {0.863, 0.881, 0.885, 0.878, 0.934, 0.852, 0.897, 0.910, 0.950, 0.923, 0.897, 0.842, 0.954, 0.914, 0.912, 0.872, 0.891, 0.878, 0.917, 0.873, 0.883, 0.924, 0.838, 0.947, 0.884, 0.788, 0.921, 0.938, 0.955, 0.881, 0.854, 0.866, 0.941, 0.885, 0.918, 0.956, 0.884, 0.867, 0.900, 0.900, 0.895, 0.854, 0.857, 0.972, 0.876, 0.863, 0.851, 0.853, 0.874, 0.863, 0.881, 0.936, 0.835, 0.846, 0.974, 0.986, 0.868, 0.787, 0.781, 0.931, 0.940, 0.864, 0.995, 0.883};
  double SBmtwMix = SBmtwcorrMix[searchbin];
  return SBmtwMix;
}

double Efficiency::SBisotrkeffTT(int searchbin){
  //  double SBisotrkeffTT[45] = {0.274, 0.281, 0.285, 0.290, 0.312, 0.301, 0.273, 0.334, 0.346, 0.307, 0.333, 0.255, 0.278, 0.276, 0.340, 0.322, 0.285, 0.299, 0.303, 0.368, 0.337, 0.275, 0.248, 0.299, 0.386, 0.436, 0.391, 0.324, 0.358, 0.367, 0.421, 0.366, 0.395, 0.385, 0.473, 0.315, 0.354, 0.444, 0.381, 0.317, 0.321, 0.371, 0.392, 0.400, 0.200};//ttbarSingleLep

  double SBisotrkeffTT[45] = {0.308, 0.319, 0.307, 0.290, 0.331, 0.341, 0.298, 0.337, 0.363, 0.320, 0.343, 0.290, 0.303, 0.286, 0.340, 0.340, 0.322, 0.310, 0.303, 0.368, 0.321, 0.296, 0.252, 0.295, 0.390, 0.436, 0.388, 0.330, 0.367, 0.371, 0.418, 0.367, 0.396, 0.381, 0.473, 0.333, 0.346, 0.438, 0.383, 0.317, 0.318, 0.374, 0.392, 0.400, 0.200};//ttbarSingleLep with HT cut

 double SBisotrkeff = SBisotrkeffTT[searchbin];
 return SBisotrkeff;
}

double Efficiency::SBisotrkeffWJet(int searchbin){

  double SBisotrkeffWJet[45] = {0.360, 0.346, 0.321, 0.261, 0.441, 0.257, 0.232, 0.301, 0.355, 0.471, 0.459, 0.404, 0.394, 0.403, 0.336, 0.155, 0.284, 0.158, 0.642, 0.596, 0.651, 0.397, 0.068, 0.384, 0.328, 0.640, 0.589, 0.372, 0.159, 0.267, 0.766, 0.292, 0.219, 0.219, 0.219, 0.818, 0.679, 0.934, 0.312, 0.312, 0.520, 0.520, 0.574, 0.574, 0.574};//WJets with HT cut                                                      

double SBisotrkeff = SBisotrkeffWJet[searchbin];
return SBisotrkeff;
}


double Efficiency::isotrkeffTT(int njetbin){

  //  double isotrk[7] = {0.175, 0.267, 0.311, 0.365, 0.375, 0.310, 0.457};//ttbarInc
    double isotrk[7] = {0.163, 0.288, 0.335, 0.359, 0.376, 0.383, 0.381};//ttbarSingleLep
  double eff = isotrk[njetbin];
  return eff;
}
double Efficiency::taumucorTT(int searchbin){
  //  double SBtaumucorr[45]= {0.154, 0.181, 0.211, 0.267, 0.165, 0.186, 0.207, 0.259, 0.172, 0.187, 0.248, 0.153, 0.191, 0.209, 0.263, 0.159, 0.214, 0.250, 0.305, 0.170, 0.260, 0.156, 0.201, 0.216, 0.134, 0.162, 0.159, 0.136, 0.177, 0.185, 0.180, 0.239, 0.139, 0.149, 0.197, 0.140, 0.149, 0.175, 0.219, 0.214, 0.137, 0.184, 0.139, 0.134, 0.099};
  // double SBtaumucorr[45]= {0.142, 0.170, 0.205, 0.270, 0.149, 0.174, 0.205, 0.260, 0.162, 0.181, 0.249, 0.141, 0.179, 0.202, 0.263, 0.150, 0.203, 0.247, 0.307, 0.172, 0.261, 0.145, 0.199, 0.215, 0.134, 0.162, 0.160, 0.136, 0.174, 0.186, 0.171, 0.236, 0.137, 0.149, 0.198, 0.137, 0.149, 0.174, 0.219, 0.216, 0.135, 0.185, 0.140, 0.134, 0.099};//with HT cut
  double SBtaumucorr[45]= {0.143, 0.172, 0.207, 0.266, 0.151, 0.176, 0.207, 0.260, 0.163, 0.182, 0.251, 0.142, 0.181, 0.203, 0.265, 0.151, 0.205, 0.248, 0.306, 0.172, 0.259, 0.146, 0.199, 0.214, 0.133, 0.163, 0.159, 0.138, 0.176, 0.185, 0.170, 0.236, 0.136, 0.148, 0.194, 0.137, 0.149, 0.175, 0.227, 0.217, 0.134, 0.187, 0.141, 0.135, 0.102};//with HT cut & new temp.
  double SBtaumucor = SBtaumucorr[searchbin];
  return SBtaumucor;
}
 double Efficiency::taumucorWJet(int searchbin){
   double SBtaumucorr[45]= {0.144, 0.190, 0.203, 0.296, 0.159, 0.161, 0.233, 0.262, 0.152, 0.163, 0.227, 0.144, 0.138, 0.170, 0.186, 0.101, 0.193, 0.171, 0.108, 0.060, 0.155, 0.151, 0.126, 0.142, 0.153, 0.146, 0.243, 0.163, 0.267, 0.120, 0.092, 0.188, 0.166, 0.160, 0.348, 0.328, 0.132, 0.284, 0.081, 0.210, 0.115, 0.172, 0.175, 0.041, 0.073};
   double SBtaumucor = SBtaumucorr[searchbin];
   return SBtaumucor;
 }
#endif
