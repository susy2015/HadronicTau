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
  static double SBtaumucorMix(int searchbin);
  static double taumucorMix(int njetbin, int metbin);
  static double taumucorMix3D(int nhtbin, int njetbin, int metbin);
  static double accTT(int njetbin, int mt2bin);
  static double accWJet(int njetbin, int mt2bin);
  static double accMix_NjetMT2(int njetbin, int mt2bin);
  static double accMix_NjetMet(int njetbin, int metbin);
  static double accMix_NjetMT2Met(int metbin, int njetbin, int mt2bin);
  static double SBaccMix(int searchbin);
  static int Ptbin(double pt); 
  static int Ptbin1(double pt);
  static int Actbin(double act);
  static int Njetbin(int njet);
  static int NJetbin(int nJet);
  static int NBjetbin(int nBjet);
  static int Htbin(double ht);
  static int MT2bin(double mt2);
  static int MetFinerbin(double met);
  static int metbin(double met);
  static int topbin(int top);
  static double mistag(int ptbin);  
  static double SBmtwTT(int searchbin);
  static double SBmtwWJet(int searchbin);
  static double SBmtwMix(int searchbin);
  static double mtwMix(int njetbin, int metbin);
  static double SBisotrkeffTT(int searchbin);
  static double SBisotrkeffWJet(int searchbin);
  static double SBisotrkeffMix(int searchbin);
  static double isotrkeffMix(int njetbin);
  static double isotrkeffMix_NjetNbjet(int njetbin, int nbjetbin);

  static double accErr(int searchbin);
  static double mtwErr(int njetbin, int metbin);
  static double mtwMetjecUp(int njetbin, int metbin);
  static double mtwMetjecLow(int njetbin, int metbin);
  static double mtwMetjerUp(int njetbin, int metbin);
  static double mtwMetjerLow(int njetbin, int metbin);
  static double taumucorErr(int njetbin, int metbin);
  static double isotrkeffErr(int njetbin, int nbjetbin);
  static double accpdfCentral(int searchbin);
  static double accpdfUp(int searchbin);
  static double accpdfDown(int searchbin);
  static double accscaleUp(int searchbin);
  static double accscaleDown(int searchbin);

  //static double MT2Bin_mtwcorr() {return 0.533;}//TTbarDiLep
  static double MT2Bin_mtwcorr() {return 0.921;}//Mix v3
  //static double MT2Bin_isotrkcorr(){return 0.294;} //(ttbarDiLep)
  static double MT2Bin_isotrkcorr(){return 0.352;} //Mix v3
   //static double MT2Bin_taumucorr(){return 0.104;}//TTbarDiLep
  static double MT2Bin_taumucorr(){return 0.154;}//(Mix v3)

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
  if(nBjet>=2) bin =1;
  return bin;
}

int Efficiency::MT2bin(double mt2){
  int bin = 0;
  if(mt2>=250 && mt2<300)bin =1;
  if(mt2>=300 && mt2<350)bin =2;
  if(mt2>=350 && mt2<400)bin =3;
  if(mt2>=400)bin =4;
  return bin;
}
int Efficiency::metbin(double met){
  int bin = 0;
  if(met>=275 && met<350)bin =1;
  if(met>=350 && met<450)bin =2;
  if(met>=450)bin =3;
  return bin;
}
int Efficiency::topbin(int top){
  int bin = 0;
  if(top>=2)bin =1;
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
int Efficiency::Htbin(double ht){
  int bin = 0;
  if(ht>=550 && ht<650) bin =1;
  if(ht>=650 && ht<800) bin =2;
  if(ht>=800) bin =3;
  return bin;
}
double Efficiency::reco(int ptbin, int actbin){
  double mu_recoeff[7][5] = {{0.960853,0.976249,0.967035,0.958783,0.947575},{0.982002,0.98174,0.986945,0.982036,0.964605},{0.986154,0.983751,0.985766,0.979747,0.974022},{0.993597,0.993206,0.989084,0.98649,0.968653},{0.993947,0.994536,0.985616,0.987728,0.971316},{0.993105,0.985795,0.982129,0.968903,0.97366},{0.967846,0.968589,0.965088,0.956531,0.971364}};

 double reco = mu_recoeff[ptbin][actbin];
 return reco;
}

double Efficiency::iso(int ptbin, int actbin){
  double mu_isoeff[7][5] = {{0.959676,0.954265,0.939252,0.84329,0.473471},{0.977308,0.972333,0.960718,0.823429,0.492998},{0.986146,0.977891,0.963588,0.811057,0.644872},{0.988595,0.98678,0.969607,0.830557,0.647136},{0.993762,0.995003,0.984036,0.89737,0.791012},{0.995213,0.995036,0.99284,0.961394,0.849533},{0.999452,0.996991,0.996952,0.982389,0.962313}};

  double iso = mu_isoeff[ptbin][actbin];
  return iso;
 }


double Efficiency::acc(int njetbin){
  double mu_acc[6] = {0.479, 0.714, 0.794, 0.837, 0.860, 0.884};//ttbarSingleLep
  double acc = mu_acc[njetbin];
  return acc;
}

double Efficiency::accTT(int njetbin, int mt2bin){

  double mu_accTT[6][5] = {{0.638848, 0.620549, 0.647076, 0.671148, 0.680209},{0.716494, 0.806024, 0.738966, 0.695279, 0.727193},{0.801679, 0.862727, 0.797992, 0.779789, 0.811038},{0.855841, 0.886209, 0.837599, 0.808209, 0.829777},{0.868522, 0.895633, 0.862863, 0.845216, 0.840523},{0.824104, 0.903998, 0.891471, 0.869352, 0.858733}};//SingleLep v3 TTtemp

  double acc = mu_accTT[njetbin][mt2bin];
  return acc;
}

double Efficiency::accWJet(int njetbin, int mt2bin){

  double mu_accWJet[6][5] = {{0.976083, 0.688995, 0.713259, 0.710482, 0.752297},{0.889607, 0.831978, 0.705857, 0.629649, 0.627766},{0.843985, 0.79728, 0.73966, 0.754925, 0.740138},{0.973495, 0.839425, 0.773697, 0.786713, 0.830726},{0.975866, 0.734035, 0.770906, 0.769405, 0.819778},{0.990525, 0.917348, 0.757306, 0.812438, 0.91809}};//WJet v3 WJettemp
 double acc = mu_accWJet[njetbin][mt2bin];
  return acc;
}

double Efficiency::accMix_NjetMT2(int njetbin, int mt2bin){
  //  double mu_accMix[6][3] = {{0.651664, 0.608061, 0.697368} , {0.743137, 0.699504, 0.739363} , {0.797933, 0.780114, 0.819086} , {0.834814, 0.810208, 0.850696} , {0.861942, 0.839966, 0.841126} , {0.891635, 0.870958, 0.860843}};
  double mu_accMix[6][5] = {{0.637, 0.696, 0.629, 0.587, 0.706} , {0.761, 0.698, 0.691, 0.705, 0.741} , {0.811, 0.775, 0.780, 0.772, 0.815} , {0.839, 0.829, 0.810, 0.805, 0.851} , {0.864, 0.858, 0.845, 0.840, 0.847} , {0.895, 0.887, 0.872, 0.868, 0.867}};
 double acc = mu_accMix[njetbin][mt2bin];
  return acc;
}
double Efficiency::accMix_NjetMet(int njetbin, int metbin){
  double mu_accMix[6][4] = {{0.620, 0.677, 0.651, 0.703} , {0.726, 0.719, 0.747, 0.784} , {0.790, 0.784, 0.822, 0.823} , {0.818, 0.831, 0.871, 0.894} , {0.845, 0.870, 0.878, 0.895} , {0.877, 0.882, 0.913, 0.928}};
 double acc = mu_accMix[njetbin][metbin];
  return acc;
}
double Efficiency::accMix_NjetMT2Met(int metbin, int njetbin, int mt2bin){
  //  double mu_accMix[4][6][5] ={{{0.633, 0.687, 0.527, 0.254, 0.721} , {0.748, 0.699, 0.681, 0.730, 0.896} , {0.806, 0.764, 0.796, 0.777, 0.791} , {0.830, 0.819, 0.794, 0.787, 0.769} , {0.859, 0.841, 0.837, 0.802, 0.804} , {0.887, 0.878, 0.859, 0.851, 0.866}}, {{0.671, 0.720, 0.658, 0.666, 0.636} , {0.768, 0.683, 0.694, 0.640, 0.776} , {0.808, 0.784, 0.748, 0.756, 0.851} , {0.852, 0.838, 0.829, 0.777, 0.823} , {0.880, 0.895, 0.856, 0.841, 0.832} , {0.903, 0.887, 0.871, 0.869, 0.831}}, {{0.575, 0.694, 0.710, 0.678, 0.771} , {0.836, 0.700, 0.686, 0.750, 0.698} , {0.878, 0.839, 0.764, 0.800, 0.805} , {0.896, 0.887, 0.862, 0.874, 0.835} , {0.868, 0.913, 0.855, 0.890, 0.862} , {0.938, 0.937, 0.919, 0.874, 0.879}}, {{0.726, 0.688, 0.761, 0.681, 0.690} , {0.809, 0.822, 0.812, 0.834, 0.758} , {0.852, 0.878, 0.828, 0.780, 0.805} , {0.888, 0.909, 0.856, 0.869, 0.912} , {0.913, 0.909, 0.909, 0.957, 0.856} , {0.958, 0.967, 0.927, 0.957, 0.891}}};
    double mu_accMix[4][6][5] ={{{0.640, 0.677, 0.511, 0.271, 0.954} , {0.752, 0.695, 0.673, 0.738, 0.864} , {0.812, 0.764, 0.796, 0.781, 0.782} , {0.832, 0.821, 0.796, 0.786, 0.802} , {0.861, 0.842, 0.840, 0.804, 0.826} , {0.894, 0.880, 0.867, 0.852, 0.895}},{{0.671, 0.708, 0.647, 0.642, 0.717} , {0.770, 0.685, 0.690, 0.650, 0.762} , {0.806, 0.788, 0.754, 0.762, 0.864} , {0.852, 0.836, 0.831, 0.784, 0.824} , {0.877, 0.894, 0.865, 0.839, 0.838} , {0.901, 0.889, 0.874, 0.876, 0.837}},{{0.561, 0.658, 0.694, 0.664, 0.773} , {0.834, 0.705, 0.682, 0.740, 0.688} , {0.883, 0.834, 0.763, 0.804, 0.802} , {0.892, 0.884, 0.863, 0.877, 0.828} , {0.867, 0.908, 0.853, 0.879, 0.850} , {0.936, 0.936, 0.916, 0.877, 0.874}},{{0.714, 0.652, 0.780, 0.696, 0.691} , {0.802, 0.811, 0.804, 0.817, 0.756} , {0.848, 0.880, 0.819, 0.799, 0.808} , {0.884, 0.906, 0.864, 0.861, 0.908} , {0.937, 0.918, 0.908, 0.968, 0.842} , {0.957, 0.964, 0.923, 0.956, 0.884}}};//applying isotrk veto to remaing part of mu CS.


 double acc = mu_accMix[metbin][njetbin][mt2bin];
 return acc;
}

double Efficiency::SBaccMix(int searchbin){
  double mu_acc[45]={0.777, 0.773, 0.835, 0.831, 0.790, 0.765, 0.769, 0.832, 0.857, 0.781, 0.809, 0.771, 0.789, 0.794, 0.825, 0.765, 0.768, 0.756, 0.780, 0.884, 0.877, 0.816, 0.840, 0.843, 0.882, 0.929, 0.941, 0.831, 0.832, 0.878, 0.813, 0.819, 0.895, 0.930, 0.944, 0.826, 0.857, 0.907, 0.852, 0.872, 0.883, 0.879, 0.911, 0.933, 0.954};
  //double mu_acc[37]={0.773, 0.776, 0.838, 0.846, 0.789, 0.767, 0.783, 0.826, 0.860, 0.776, 0.785, 0.780, 0.799, 0.811, 0.840, 0.756, 0.775, 0.773, 0.801, 0.849, 0.807, 0.885, 0.926, 0.951, 0.828, 0.831, 0.880, 0.809, 0.817, 0.901, 0.926, 0.942, 0.825, 0.860, 0.901, 0.800, 0.861};
  double acc = mu_acc[searchbin];
  return acc;
}
double Efficiency::mistag(int ptbin){
  double mistag[18] = {0, 0.0584243, 0.0679155, 0.0877595, 0.0874307, 0.101414, 0.098585, 0.0853752, 0.0818409, 0.0654232, 0.0767843, 0.0524587, 0.0801002, 0.0743795, 0.0597973, 0.0497092, 0.0444659, 0.0444659};//WJet v3

  double rate = mistag[ptbin];
  return rate;
}

double Efficiency::SBmtwTT(int searchbin){
  double SBmtwcorr[45] = {0.904, 0.930, 0.941, 0.941, 0.914, 0.934, 0.953, 0.921, 0.919, 0.948, 0.967, 0.914, 0.930, 0.932, 0.946, 0.908, 0.941, 0.936, 0.953, 0.933, 0.940, 0.893, 0.927, 0.926, 0.908, 0.923, 0.944, 0.904, 0.927, 0.936, 0.897, 0.944, 0.906, 0.925, 0.921, 0.909, 0.922, 0.913, 0.880, 0.937, 0.884, 0.917, 0.897, 0.895, 0.866};//TTbarSingleLep v3 TTtemp

  double SBmtw = SBmtwcorr[searchbin];
  return SBmtw;
}

double Efficiency::SBmtwWJet(int searchbin){

  double SBmtwcorrWJet[45] = {0.879, 0.894, 0.880, 0.904, 0.870, 0.871, 0.896, 0.921, 0.812, 0.865, 0.891, 0.898, 0.899, 0.895, 0.887, 0.898, 0.858, 0.923, 0.868, 0.883, 0.902, 0.811, 0.805, 0.867, 0.879, 0.866, 0.973, 0.855, 0.859, 0.917, 0.886, 0.874, 0.881, 0.897, 0.960, 0.726, 0.857, 0.965, 0.936, 0.908, 0.956, 0.811, 0.899, 0.901, 0.857};//v3 WJettemp

  double SBmtwWJet = SBmtwcorrWJet[searchbin];
  return SBmtwWJet;
}

double Efficiency::SBmtwMix(int searchbin){
  double SBmtwcorrMix[45] = {0.896, 0.916, 0.915, 0.924, 0.892, 0.919, 0.929, 0.925, 0.901, 0.918, 0.930, 0.909, 0.927, 0.926, 0.917, 0.889, 0.925, 0.937, 0.912, 0.933, 0.940, 0.890, 0.917, 0.906, 0.903, 0.916, 0.936, 0.895, 0.929, 0.934, 0.911, 0.911, 0.907, 0.919, 0.916, 0.898, 0.913, 0.918, 0.875, 0.937, 0.883, 0.917, 0.903, 0.890, 0.868};//mix (ttbar singlelep, Wjets, single top)
  double SBmtwMix = SBmtwcorrMix[searchbin];
  return SBmtwMix;
}
double Efficiency::mtwMix(int njetbin, int metbin){
  double mtwcorrMix[6][4] = {{0.892981, 0.875669, 0.88925, 0.943},{0.905773, 0.929934, 0.938666, 0.951932},{0.904658, 0.929749, 0.937562, 0.930164},{0.906773, 0.916012, 0.921934, 0.924404},{0.89874, 0.919228, 0.921755, 0.930223},{0.878255, 0.906088, 0.908429, 0.921293}};
  //double mtwcorrMix[6][7] = {{0.896671, 0.880949, 0.904045, 0.857816, 0.854162, 0.950436, 0.943},{0.900664, 0.921451, 0.927681, 0.931548, 0.930469, 0.953685, 0.951932},{0.900509, 0.917344, 0.929771, 0.929731, 0.939245, 0.934324, 0.930164},{0.904031, 0.916204, 0.914681, 0.917092, 0.919725, 0.926128, 0.924404},{0.895975, 0.908408, 0.910603, 0.92604, 0.917871, 0.929405, 0.930223},{0.875962, 0.885616, 0.901179, 0.910047, 0.906377, 0.911974, 0.921293}};
double mtwMix = mtwcorrMix[njetbin][metbin];
  return mtwMix;
}
double Efficiency::SBisotrkeffTT(int searchbin){
  double SBisotrkeffTT[45] = {0.307, 0.313, 0.324, 0.287, 0.340, 0.342, 0.305, 0.336, 0.362, 0.321, 0.329, 0.289, 0.297, 0.275, 0.340, 0.350, 0.315, 0.299, 0.292, 0.395, 0.338, 0.293, 0.259, 0.297, 0.395, 0.430, 0.401, 0.337, 0.359, 0.378, 0.393, 0.382, 0.390, 0.386, 0.472, 0.324, 0.330, 0.448, 0.372, 0.321, 0.311, 0.374, 0.396, 0.391, 0.204};//ttbarSingleLep v3
 double SBisotrkeff = SBisotrkeffTT[searchbin];
 return SBisotrkeff;
}

double Efficiency::SBisotrkeffWJet(int searchbin){

  double SBisotrkeffWJet[45] = {0.346, 0.330, 0.336, 0.242, 0.309, 0.283, 0.324, 0.297, 0.524, 0.372, 0.487, 0.382, 0.398, 0.448, 0.264, 0.266, 0.301, 0.119, 0.657, 0.629, 0.720, 0.640, 0.078, 0.655, 0.299, 0.415, 0.666, 0.443, 0.256, 0.310, 0.575, 0.333, 0.200, 0.200, 0.200, 0.274, 0.905, 0.467, 0.314, 0.314, 0.311, 0.311, 0.597, 0.597, 0.597};//WJets v3  

double SBisotrkeff = SBisotrkeffWJet[searchbin];
return SBisotrkeff;
}

double Efficiency::SBisotrkeffMix(int searchbin){

  double SBisotrkeffMix[45] = {0.317, 0.324, 0.306, 0.274, 0.367, 0.330, 0.374, 0.300, 0.422, 0.335, 0.397, 0.297, 0.306, 0.300, 0.354, 0.406, 0.333, 0.266, 0.317, 0.384, 0.535, 0.303, 0.259, 0.245, 0.389, 0.416, 0.442, 0.352, 0.354, 0.361, 0.396, 0.362, 0.389, 0.377, 0.489, 0.318, 0.373, 0.506, 0.341, 0.323, 0.303, 0.386, 0.380, 0.395, 0.353};//mix (ttbar singlelep, Wjets, single top)
double SBisotrkeff = SBisotrkeffMix[searchbin];
return SBisotrkeff;
}

double Efficiency::isotrkeffMix(int njetbin){
  double isotrk[6] = {0.185, 0.315, 0.343, 0.369, 0.373, 0.388};//mix (ttbar singlelep, Wjets, single top)
  double eff = isotrk[njetbin];
  return eff;
}
double Efficiency::isotrkeffMix_NjetNbjet(int njetbin, int nbjetbin){
   double isotrk[6][3] = {{0.158, 0.255, 0.333},{0.330, 0.271, 0.279},{0.329, 0.323, 0.247},{0.365, 0.358, 0.299},{0.368, 0.368, 0.305},{0.398, 0.392, 0.347}};//mix (ttbar singlelep, Wjets, single top)
   //double isotrk[6][2] = {{0.159, 0.251},{0.329, 0.268},{0.329, 0.316},{0.366, 0.349},{0.369, 0.356},{0.395, 0.383}};//mix (ttbar singlelep, Wjets, single top)

  double eff = isotrk[njetbin][nbjetbin];
  return eff;
}

double Efficiency::taumucorTT(int searchbin){
   double SBtaumucorr[45]= {0.141, 0.169, 0.194, 0.265, 0.146, 0.165, 0.202, 0.246, 0.153, 0.168, 0.253, 0.140, 0.175, 0.192, 0.253, 0.152, 0.191, 0.242, 0.279, 0.168, 0.265, 0.143, 0.185, 0.204, 0.132, 0.159, 0.148, 0.134, 0.170, 0.182, 0.165, 0.241, 0.136, 0.150, 0.193, 0.135, 0.154, 0.173, 0.200, 0.187, 0.133, 0.186, 0.145, 0.131, 0.120};//SingleLep v3 TTtemp.
  double SBtaumucor = SBtaumucorr[searchbin];
  return SBtaumucor;
}
double Efficiency::taumucorWJet(int searchbin){

  double SBtaumucorr[45]= {0.153, 0.188, 0.200, 0.272, 0.177, 0.154, 0.225, 0.187, 0.130, 0.202, 0.242, 0.133, 0.146, 0.173, 0.245, 0.084, 0.231, 0.167, 0.091, 0.049, 0.077, 0.147, 0.192, 0.106, 0.172, 0.133, 0.233, 0.161, 0.240, 0.155, 0.141, 0.168, 0.179, 0.101, 0.305, 0.256, 0.073, 0.238, 0.414, 0.173, 0.085, 0.090, 0.147, 0.030, 0.003};//v3 WJet temp

  double SBtaumucor = SBtaumucorr[searchbin];
  return SBtaumucor;
}
double Efficiency::SBtaumucorMix(int searchbin){
  double SBtaumucorr[45]= {0.143, 0.166, 0.185, 0.259, 0.151, 0.163, 0.209, 0.258, 0.178, 0.176, 0.223, 0.141, 0.182, 0.186, 0.263, 0.166, 0.203, 0.237, 0.322, 0.139, 0.209, 0.146, 0.202, 0.235, 0.135, 0.158, 0.162, 0.138, 0.175, 0.202, 0.168, 0.268, 0.131, 0.148, 0.197, 0.129, 0.150, 0.191, 0.198, 0.255, 0.133, 0.193, 0.150, 0.134, 0.197};//mix (ttbar singlelep, Wjets, single top)
  double SBtaumucor = SBtaumucorr[searchbin];
  return SBtaumucor;
}
double Efficiency::taumucorMix(int njetbin, int metbin){
  double SBtaumucorr[6][4]= {{0.180573, 0.152611, 0.230756, 0.245318},{0.144058, 0.19431, 0.220997, 0.268273},{0.144386, 0.174184, 0.192034, 0.244345},{0.1407, 0.169858, 0.182826, 0.249126},{0.140436, 0.158769, 0.189431, 0.267826},{0.128646, 0.140671, 0.164106, 0.197228}};
  //double SBtaumucorr[6][7]= {{0.159688, 0.241997, 0.161755, 0.146754, 0.195423, 0.285472, 0.245318},{0.133746, 0.17422, 0.201774, 0.188872, 0.216227, 0.229589, 0.268273},{0.141045, 0.154443, 0.163392, 0.182969, 0.18774, 0.200168, 0.244345},{0.136615, 0.15446, 0.163803, 0.174709, 0.180378, 0.187434, 0.249126},{0.137218, 0.1515, 0.164532, 0.154161, 0.17295, 0.220046, 0.267826},{0.126285, 0.136141, 0.13783, 0.142949, 0.159365, 0.17217, 0.197228}};

  double SBtaumucor = SBtaumucorr[njetbin][metbin];
  return SBtaumucor;
}

double Efficiency::taumucorMix3D(int nbjetbin, int njetbin, int metbin){
  double SBtaumucorr[4][6][4] = {{{0.166, 0.154, 0.294, 0.335},{0.149, 0.204, 0.260, 0.134},{0.151, 0.190, 0.222, 0.208},{0.154, 0.190, 0.135, 0.257},{0.147, 0.212, 0.148, 0.148},{0.140, 0.082, 0.437, 0.437}},{{0.145, 0.124, 0.118, 0.245},{0.150, 0.217, 0.212, 0.269},{0.151, 0.182, 0.223, 0.172},{0.144, 0.185, 0.211, 0.116},{0.150, 0.192, 0.220, 0.184},{0.123, 0.114, 0.122, 0.115}},{{0.315, 0.205, 0.218, 0.220},{0.134, 0.164, 0.177, 0.227},{0.139, 0.163, 0.189, 0.278},{0.147, 0.180, 0.199, 0.277},{0.143, 0.148, 0.199, 0.199},{0.142, 0.169, 0.177, 0.356}},{{0.115, 0.126, 0.264, 0.246},{0.134, 0.176, 0.242, 0.298},{0.130, 0.163, 0.157, 0.241},{0.122, 0.142, 0.165, 0.249},{0.132, 0.151, 0.183, 0.278},{0.124, 0.136, 0.162, 0.191}}};
  double SBtaumucor = SBtaumucorr[nbjetbin][njetbin][metbin];
  return SBtaumucor;
}

double Efficiency::accErr(int searchbin){
   double mu_accMix[45] ={0.005, 0.008, 0.012, 0.024, 0.008, 0.013, 0.025, 0.026, 0.017, 0.022, 0.027, 0.005, 0.007, 0.012, 0.023, 0.017, 0.019, 0.020, 0.031, 0.017, 0.026, 0.009, 0.022, 0.048, 0.005, 0.007, 0.010, 0.008, 0.014, 0.015, 0.023, 0.020, 0.005, 0.007, 0.011, 0.010, 0.013, 0.016, 0.027, 0.021, 0.010, 0.018, 0.016, 0.014, 0.023};
   //double mu_accMix[37] ={0.005, 0.008, 0.011, 0.021, 0.008, 0.012, 0.022, 0.024, 0.016, 0.020, 0.026, 0.004, 0.007, 0.011, 0.020, 0.015, 0.017, 0.018, 0.027, 0.025, 0.045, 0.004, 0.007, 0.009, 0.008, 0.014, 0.015, 0.023, 0.020, 0.004, 0.007, 0.010, 0.009, 0.012, 0.015, 0.027, 0.020};
 double err = mu_accMix[searchbin];
 return err;
}
double Efficiency::accpdfUp(int searchbin){
  double mu_accMix[45] ={0.771, 0.776, 0.845, 0.849, 0.787, 0.764, 0.781, 0.831, 0.853, 0.772, 0.783, 0.771, 0.792, 0.792, 0.805, 0.757, 0.761, 0.769, 0.773, 0.851, 0.865, 0.812, 0.842, 0.839, 0.883, 0.919, 0.945, 0.825, 0.800, 0.850, 0.805, 0.794, 0.896, 0.931, 0.941, 0.805, 0.857, 0.904, 0.847, 0.873, 0.882, 0.884, 0.909, 0.933, 0.955};
  //double mu_accMix[37] ={0.771, 0.776, 0.845, 0.849, 0.787, 0.764, 0.781, 0.831, 0.853, 0.772, 0.783, 0.778, 0.798, 0.809, 0.827, 0.755, 0.773, 0.771, 0.792, 0.848, 0.811, 0.884, 0.917, 0.947, 0.825, 0.801, 0.850, 0.806, 0.796, 0.901, 0.926, 0.942, 0.809, 0.859, 0.902, 0.800, 0.862};

  double err = mu_accMix[searchbin];
  return err;
}
double Efficiency::accpdfDown(int searchbin){
  double mu_accMix[45] ={0.775, 0.777, 0.830, 0.844, 0.791, 0.770, 0.785, 0.820, 0.868, 0.779, 0.788, 0.773, 0.792, 0.799, 0.836, 0.761, 0.765, 0.773, 0.796, 0.855, 0.861, 0.822, 0.841, 0.834, 0.884, 0.938, 0.953, 0.832, 0.869, 0.916, 0.812, 0.842, 0.895, 0.931, 0.944, 0.848, 0.860, 0.903, 0.848, 0.870, 0.880, 0.880, 0.911, 0.933, 0.953};
  //double mu_accMix[37] ={0.775, 0.777, 0.830, 0.844, 0.791, 0.770, 0.785, 0.820, 0.868, 0.779, 0.788, 0.782, 0.799, 0.814, 0.854, 0.758, 0.776, 0.775, 0.810, 0.851, 0.802, 0.885, 0.936, 0.955, 0.832, 0.868, 0.917, 0.813, 0.843, 0.901, 0.926, 0.943, 0.844, 0.861, 0.901, 0.801, 0.860};
  double err = mu_accMix[searchbin];
  return err;
}
double Efficiency::accscaleUp(int searchbin){
  double mu_accMix[45] ={0.772, 0.773, 0.835, 0.839, 0.786, 0.765, 0.781, 0.821, 0.856, 0.775, 0.782, 0.770, 0.791, 0.793, 0.817, 0.762, 0.765, 0.766, 0.782, 0.847, 0.856, 0.815, 0.845, 0.843, 0.883, 0.926, 0.947, 0.827, 0.829, 0.878, 0.807, 0.813, 0.894, 0.930, 0.941, 0.823, 0.855, 0.903, 0.849, 0.869, 0.879, 0.883, 0.911, 0.932, 0.952};
  //double mu_accMix[37] ={0.772, 0.773, 0.835, 0.839, 0.786, 0.765, 0.781, 0.821, 0.856, 0.775, 0.782, 0.778, 0.798, 0.810, 0.838, 0.760, 0.777, 0.769, 0.800, 0.844, 0.803, 0.884, 0.924, 0.949, 0.827, 0.829, 0.878, 0.808, 0.815, 0.900, 0.925, 0.941, 0.824, 0.857, 0.901, 0.801, 0.859};
  double err = mu_accMix[searchbin];
  return err;
}
double Efficiency::accscaleDown(int searchbin){
  double mu_accMix[45] ={0.774, 0.780, 0.842, 0.854, 0.792, 0.769, 0.785, 0.833, 0.865, 0.777, 0.789, 0.773, 0.793, 0.798, 0.821, 0.754, 0.760, 0.775, 0.785, 0.860, 0.871, 0.819, 0.838, 0.832, 0.885, 0.930, 0.950, 0.829, 0.833, 0.883, 0.810, 0.818, 0.897, 0.932, 0.943, 0.828, 0.861, 0.904, 0.847, 0.874, 0.883, 0.882, 0.910, 0.934, 0.956};
  //double mu_accMix[37] ={0.774, 0.780, 0.842, 0.854, 0.792, 0.769, 0.785, 0.833, 0.865, 0.777, 0.789, 0.782, 0.800, 0.814, 0.841, 0.752, 0.772, 0.777, 0.801, 0.855, 0.812, 0.885, 0.928, 0.952, 0.829, 0.833, 0.883, 0.811, 0.820, 0.903, 0.928, 0.943, 0.827, 0.862, 0.902, 0.800, 0.863};
  double err = mu_accMix[searchbin];
  return err;
}
double Efficiency::mtwErr(int njetbin, int metbin){
  double mtwcorrMix[6][4] = {{0.00732952, 0.0137408, 0.0158642, 0.00862698},{0.00197141, 0.00224627, 0.00313586, 0.00343604},{0.0011598, 0.00141399, 0.00204034, 0.00393123},{0.000933405, 0.00171172, 0.00290407, 0.00431023},{0.00111457, 0.00190538, 0.0030743, 0.00302469},{0.00143327, 0.00235669, 0.00309154, 0.0036149}};
 
double err = mtwcorrMix[njetbin][metbin];
  return err;
}
double Efficiency::mtwMetjecUp(int njetbin, int metbin){
  double mtwcorrMix[6][4] = {{0.892211, 0.860436, 0.882465, 0.941232},{0.904211, 0.926509, 0.934947, 0.950093},{0.899434, 0.926471, 0.933147, 0.927752},{0.903825, 0.914324, 0.919297, 0.925284},{0.894584, 0.917324, 0.919833, 0.928528},{0.876036, 0.903039, 0.905927, 0.918944}};
double err = mtwcorrMix[njetbin][metbin];
  return err;
}
double Efficiency::mtwMetjecLow(int njetbin, int metbin){
  double mtwcorrMix[6][4] = {{0.895892, 0.865782, 0.887931, 0.941946},{0.909005, 0.931887, 0.937783, 0.953297},{0.906775, 0.930936, 0.936864, 0.928527},{0.909926, 0.918055, 0.924714, 0.928581},{0.902977, 0.922159, 0.923897, 0.933045},{0.882094, 0.908944, 0.911184, 0.921319}};
  double err = mtwcorrMix[njetbin][metbin];
  return err;
}
double Efficiency::mtwMetjerUp(int njetbin, int metbin){
  double mtwcorrMix[6][4] = {{0.893836, 0.863399, 0.883185, 0.941307},{0.906864, 0.929312, 0.936586, 0.952295},{0.904598, 0.92842, 0.934822, 0.92814},{0.906683, 0.916426, 0.921418, 0.926562},{0.898518, 0.920292, 0.921933, 0.93191},{0.878692, 0.905708, 0.908713, 0.92118}};
  double err = mtwcorrMix[njetbin][metbin];
  return err;
}
double Efficiency::mtwMetjerLow(int njetbin, int metbin){
  double mtwcorrMix[6][4] = {{0.893836, 0.863399, 0.883185, 0.941307},{0.906864, 0.929312, 0.936586, 0.952295},{0.904598, 0.92842, 0.934822, 0.92814},{0.906683, 0.916426, 0.921418, 0.926562},{0.898518, 0.920292, 0.921933, 0.93191},{0.878692, 0.905708, 0.908713, 0.92118}};
  double err = mtwcorrMix[njetbin][metbin];
  return err;
}
double Efficiency::taumucorErr(int njetbin, int metbin){
  double SBtaumucorr[6][4]= {{0.00742852, 0.00708515, 0.0122999, 0.0102427},{0.00191177, 0.0034301, 0.00433467, 0.0078976},{0.00113561, 0.00189638, 0.00311547, 0.005661},{0.000973077, 0.0017064, 0.00293443, 0.00680038},{0.00120194, 0.00200321, 0.00391566, 0.0090098},{0.0011586, 0.00192362, 0.00373763, 0.00740337}};

  double err = SBtaumucorr[njetbin][metbin];
  return err;
}
double Efficiency::isotrkeffErr(int njetbin, int nbjetbin){
  double isotrk[6][3] = {{0.024, 0.040, 0.272},{0.015, 0.012, 0.045},{0.010, 0.010, 0.021},{0.009, 0.010, 0.018},{0.011, 0.012, 0.023},{0.012, 0.013, 0.028}};
  //double isotrk[6][2] = {{0.024, 0.039},{0.014, 0.011},{0.009, 0.009},{0.009, 0.009},{0.011, 0.011},{0.012, 0.012}};
  double err = isotrk[njetbin][nbjetbin];
  return err;
}
#endif
