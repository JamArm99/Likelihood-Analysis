#ifndef unbin_like_H
#define unbin_like_H
#include <TROOT.h>

//Selected parameters for likelihood
std::vector<TString> vars = {"closestPMT","dxPrevx","dyPrevy","dzPrevz","drPrevr","dt_prev_us","n100", "n100_prev","innerPE"};

//Daily rates for tank design
std::vector<double> max_rate = {8.23, 5.71, 1.96e9, 0.379, 0.606, 0.692, 0.698};

//Signal and background files as probability distribution functions
std::vector<TString> pdf_files = {"signal.root","background.root"};
std::vector<TString> pdf = {"signal", "background"};

//Data files for analysis - Order matters
std::vector<TString> data_files = {"big.root","small.root","singles.root","world.root","geo.root","17_N.root","9_Li.root"};

//Colours for stacked hist plot
std::vector<int> colours = {9,38,30,41,8,36,46};

//Legend for stacked hist plot
std::vector<TString> legend_names = {"Big Hartlepool","Small Hartlepool", "Combined Singles", "World Reactors", "Geoneutrinos","17-Nitrogen","9-Lithium"};

//Signal significance function
double sigma(double s, double b) {return sqrt(30)*s/sqrt(s+b);}

//Daily rate calculation from number of events
double rate(double r0, double n, double n0) {return r0*n/n0;}

#endif
