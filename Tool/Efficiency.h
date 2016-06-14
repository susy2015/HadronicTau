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
  static double mistag(int ptbin);  

  static int Ptbin(double pt); 
  static int Ptbin1(double pt);
  static int Actbin(double act);
  static int Njetbin(int njet);
  static int NJetbin(int nJet);
  static int NBjetbin(int nBjet);
  static int MetFinerbin(double met);
  static int metbin(double met);

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

int Efficiency::metbin(double met){
  int bin = 0;
  if(met>=275 && met<350)bin =1;
  if(met>=350 && met<450)bin =2;
  if(met>=450)bin =3;
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

double Efficiency::SBaccMix(int searchbin){
  double mu_acc[69]={0.775, 0.806, 0.797, 0.784, 0.805, 0.793, 0.786, 0.805, 0.869, 0.796, 0.772, 0.889, 0.774, 0.781, 0.823, 0.700, 0.804, 0.802, 0.848, 0.828, 0.842, 0.781, 0.862, 0.806, 0.827, 0.891, 0.931, 0.781, 0.771, 0.814, 0.881, 0.926, 0.959, 0.866, 0.801, 0.845, 0.873, 0.857, 0.828, 0.743, 0.832, 0.792, 0.885, 0.927, 0.958, 0.809, 0.832, 0.883, 0.942, 0.852, 0.942, 0.845, 0.816, 1.000, 0.889, 0.915, 1.000, 0.800, 0.799, 0.889, 0.901, 0.975, 0.864, 0.938, 0.881, 1.000, 0.943, 0.995, 1.000};
   double acc = mu_acc[searchbin];
  return acc;
}
double Efficiency::mistag(int ptbin){
  double mistag[18] = {0, 0.0584243, 0.0679155, 0.0877595, 0.0874307, 0.101414, 0.098585, 0.0853752, 0.0818409, 0.0654232, 0.0767843, 0.0524587, 0.0801002, 0.0743795, 0.0597973, 0.0497092, 0.0444659, 0.0444659};//WJet v3

  double rate = mistag[ptbin];
  return rate;
}

double Efficiency::mtwMix(int njetbin, int metbin){
  double mtwcorrMix[6][4] = {{0.892981, 0.875669, 0.88925, 0.943},{0.905773, 0.929934, 0.938666, 0.951932},{0.904658, 0.929749, 0.937562, 0.930164},{0.906773, 0.916012, 0.921934, 0.924404},{0.89874, 0.919228, 0.921755, 0.930223},{0.878255, 0.906088, 0.908429, 0.921293}};
double mtwMix = mtwcorrMix[njetbin][metbin];
  return mtwMix;
}

double Efficiency::isotrkeffMix_NjetNbjet(int njetbin, int nbjetbin){
   double isotrk[6][3] = {{0.158, 0.255, 0.333},{0.330, 0.271, 0.279},{0.329, 0.323, 0.247},{0.365, 0.358, 0.299},{0.368, 0.368, 0.305},{0.398, 0.392, 0.347}};//mix (ttbar singlelep, Wjets, single top)
  double eff = isotrk[njetbin][nbjetbin];
  return eff;
}

double Efficiency::taumucorMix(int njetbin, int metbin){
  double SBtaumucorr[6][4]= {{0.180573, 0.152611, 0.230756, 0.245318},{0.144058, 0.19431, 0.220997, 0.268273},{0.144386, 0.174184, 0.192034, 0.244345},{0.1407, 0.169858, 0.182826, 0.249126},{0.140436, 0.158769, 0.189431, 0.267826},{0.128646, 0.140671, 0.164106, 0.197228}};
  double SBtaumucor = SBtaumucorr[njetbin][metbin];
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
