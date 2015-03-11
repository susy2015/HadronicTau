#ifndef EFFICIENCY_H
#define EFFICIENCY_H

class Efficiency{
 public:

  static double acc(){return 0.71; }
  static double reco(int ptbin);
  static double iso(int ptbin);
  static double taumucor(int ptbin);
  static int Ptbin(double pt); 
  static int Ptbin1(double pt);
};

int Efficiency::Ptbin(double pt){
  int bin =0;
  if(pt>30. && pt<=40.) bin =1;
  if(pt>40. && pt<=50.) bin =2;
  if(pt>50. && pt<=70.) bin =3;
  if(pt>70. && pt<=100.) bin =4;
  if(pt>100.) bin =5;
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

double Efficiency::reco(int ptbin){
  double reco = 0.95;
  if(ptbin==1) reco = 0.95;
  if(ptbin==2) reco = 0.96;
  if(ptbin==3) reco = 0.95;
  if(ptbin==4) reco = 0.95;
  if(ptbin==5) reco = 0.91;

  return reco;
 }
double Efficiency::iso(int ptbin){
  double iso = 0.77;
  if(ptbin==1) iso = 0.83;
  if(ptbin==2) iso = 0.86;
  if(ptbin==3) iso = 0.88;
  if(ptbin==4) iso = 0.89;
  if(ptbin==5) iso = 0.9;

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
