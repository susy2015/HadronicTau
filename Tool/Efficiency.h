#ifndef EFFICIENCY_H
#define EFFICIENCY_H

class Efficiency{
 public:
// The acceptance is from a measurement after baseline
// Stat. uncertainty of the number need to be provided as well
  static double acc(){return 0.776;}
  static double reco(int ptbin, int actbin);
  static double iso(int ptbin, int actbin);
  static double taumucor(int ptbin);
  static int Ptbin(double pt); 
  static int Ptbin1(double pt);
  static int Actbin(double act);
  
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

int Efficiency::Ptbin1(double pt){
  int bin =0;
  if(pt>40. && pt<=100.) bin =1;
  if(pt>100. && pt<=400.) bin =2;
  if(pt>400. && pt<=600.) bin =3;
  if(pt>600.) bin =4;

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

double Efficiency::reco(int ptbin, int actbin){
  double mu_recoeff[8][8] = {{0.952769,0.906542,0.939655,0.93617,0.93038,0.914894,0.92,0.903766},{0.962189,0.966667,0.930147,0.963989,0.963768,0.915254,0.952381,0.95539},{0.971917,0.966019,0.97619,0.934911,0.964126,0.945652,0.951456,0.976048},{0.979784,0.980769,0.949239,0.944206,0.940367,0.960526,0.971429,0.975},{0.980216,0.943089,0.992958,0.986486,0.982955,0.962963,0.9625,0.953757},{0.981664,0.965517,0.949495,0.972125,0.929204,0.960452,0.971223,0.951557},{0.96882,0.990566,0.958042,0.985294,0.954545,0.954545,0.934579,0.95053},{0.925926,0.928571,0.92,0.896552,0.913043,0.9,0.918367,0.932773}};  

 double reco = mu_recoeff[ptbin][actbin];
 return reco;
}


double Efficiency::iso(int ptbin, int actbin){
  double mu_isoeff[8][8] = {{0.938462,0.804124,0.788991,0.75,0.748299,0.767442,0.724638,0.611111},{0.962771,0.896552,0.897233,0.844828,0.823308,0.820988,0.758333,0.7393},{0.991206,0.959799,0.931707,0.924051,0.902326,0.816092,0.938776,0.773006},{0.993122,0.973856,0.957219,0.918182,0.917073,0.883562,0.901961,0.866667},{0.994495,0.974138,0.950355,0.949772,0.947977,0.942308,0.883117,0.89697},{0.981322,0.982143,0.989362,0.978495,0.985714,0.982353,0.940741,0.923636},{0.981609,1,0.978102,1,0.989418,0.984127,0.96,0.95539},{0.994286,0.980769,0.985507,0.992308,0.986395,1,1,0.984985}};

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


#endif
