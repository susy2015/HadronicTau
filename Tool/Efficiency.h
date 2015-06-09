#ifndef EFFICIENCY_H
#define EFFICIENCY_H

class Efficiency{
 public:
// The acceptance is from a measurement after baseline
// Stat. uncertainty of the number need to be provided as well
  static double acc(int njetbin);
  static double reco(int ptbin, int actbin);
  static double iso(int ptbin, int actbin);
  static double taumucor(int ptbin);
  static int Ptbin(double pt); 
  static int Ptbin1(double pt);
  static int Actbin(double act);
  static int Njetbin(int njet);
  static double mistag(int ptbin);  
};

int Efficiency::Ptbin(double pt){
  int bin =0;
  if(pt>=10. && pt<20.) bin =1;
  if(pt>=20. && pt<30.) bin =2;  
  if(pt>=30. && pt<40.) bin =3;
  if(pt>=40. && pt<50.) bin =4;
  if(pt>=50. && pt<70.) bin =5;
  if(pt>=70. && pt<100.) bin =6;
  if(pt>=100.) bin =7;
  return bin;
 }
/*
int Efficiency::Ptbin1(double pt){
  int bin =0;
  if(pt>=50. && pt<100.) bin =1;
  if(pt>=100. && pt<150.) bin =2;
  if(pt>=150. && pt<180.) bin =3;
  if(pt>=180. && pt<250.) bin =4;
  if(pt>=250. && pt<300.) bin =5;
  if(pt>=300. && pt<1000.) bin =6;  
  if(pt>=1000.) bin =7;
  return bin;
 }
*/
int Efficiency::Ptbin1(double pt){
  int bin =0;
  if(pt>=25. && pt<50.) bin =1;
  if(pt>=50. && pt<75.) bin =2;
  if(pt>=75. && pt<100.) bin =3;
  if(pt>=100. && pt<125.) bin =4;
  if(pt>=125. && pt<150.) bin =5;
  if(pt>=150. && pt<175.) bin =6;
  if(pt>=175. && pt<200.) bin =7;
  if(pt>=200. && pt<1000.) bin =8;
  if(pt>1000.) bin =9;
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
  if(njet==9) bin =5;
  if(njet==10) bin =6;
  if(njet>10) bin =7;
  return bin;
}


double Efficiency::reco(int ptbin, int actbin){
  double mu_recoeff[8][8] = {{0.952769,0.906542,0.939655,0.93617,0.93038,0.914894,0.92,0.903766},{0.962189,0.966667,0.930147,0.963989,0.963768,0.915254,0.952381,0.95539},{0.971917,0.966019,0.97619,0.934911,0.964126,0.945652,0.951456,0.976048},{0.979784,0.980769,0.949239,0.944206,0.940367,0.960526,0.971429,0.975},{0.980216,0.943089,0.992958,0.986486,0.982955,0.962963,0.9625,0.953757},{0.981664,0.965517,0.949495,0.972125,0.929204,0.960452,0.971223,0.951557},{0.96882,0.990566,0.958042,0.985294,0.954545,0.954545,0.934579,0.95053},{0.925926,0.928571,0.92,0.896552,0.913043,0.9,0.918367,0.932773}};

 double reco = mu_recoeff[ptbin][actbin];
 return reco;
}


double Efficiency::iso(int ptbin, int actbin){
  double mu_isoeff[8][8] = {{0.926544,0.808081,0.77193,0.73743,0.729167,0.752688,0.71831,0.599099},{0.966805,0.884058,0.887597,0.849432,0.833948,0.825,0.762295,0.73622},{0.980271,0.945545,0.936585,0.898773,0.884259,0.8125,0.913462,0.766082},{0.995851,0.993464,0.958115,0.930435,0.925,0.899281,0.913462,0.869792},{0.987544,0.991379,0.962687,0.965686,0.955801,0.935185,0.878378,0.899441},{0.991058,0.982143,0.984127,1,0.990338,0.977528,0.945736,0.930147},{0.997625,1,0.992701,0.995,0.978947,1,0.989796,0.958801},{1,1,1,1,1,1,1,0.996865}};

  double iso = mu_isoeff[ptbin][actbin];
  return iso;
 }

double Efficiency::taumucor(int ptbin){
  double frac = 0.897;
  if(ptbin==1) frac = 0.958;
  if(ptbin==2) frac = 0.982;
  if(ptbin==3) frac = 0.984;
  if(ptbin==4) frac = 1;
  return frac;
 }

double Efficiency::acc(int njetbin){
    double mu_acc[8] = {0.5982, 0.7331, 0.7921, 0.8443, 0.8577, 0.8889, 0.8689, 0.9108};
  //double mu_acc[8] = {0.5815, 0.7284, 0.7898, 0.8431, 0.8591, 0.8882, 0.8663, 0.91};
  double acc = mu_acc[njetbin];
  return acc;
}


double Efficiency::mistag(int ptbin){
  double mistag[10] ={0.003, 0.1251, 0.1131, 0.1602, 0.156, 0.1802, 0.24, 0.4348, 0.4571, 0.4571};
  //double mistag[8] ={0.1131, 0.1273, 0.1628, 0.2758, 0.4375, 0.4286, 0.5, 0.5};
  double rate = mistag[ptbin];
  return rate;
}

#endif
