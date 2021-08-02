#include <TROOT.h>
#include "unbin_like.h"

void unbin_like(){

  gROOT->SetBatch(1);//Disbale ROOT canvas output for histogram matrix
  gErrorIgnoreLevel = kWarning;// Turn off saving file messages to terminal

  const char * pwd = gSystem->pwd();//Get current directory for image saving

  //pdf matrix of histograms for signal and background and the number of vars
  TH1F * pdf_arr[pdf.size()][vars.size()];

  //Number of bins for pdfs
  const int n_bins = 200;

  //Loop through signal and background to populate pdf histogram matrix
  for(int f=0; f<pdf.size(); f++){
    TFile * file = new TFile(pdf_files[f]);//Load files from header
    TTree * tree = (TTree*)file->Get("data");//Get the data tree
    gROOT->cd();
    TTree * cut_tree = tree->CopyTree("n100>0 && x>-9999 && dt_prev_us>0 && dt_prev_us<2000 && closestPMT>0");//Make cuts to remove unphysical parameters
    float pdf_size = cut_tree->GetEntries();
    for(int var=0; var<vars.size(); var++){
      TCanvas * c1 = new TCanvas("c1");
      float min = cut_tree->GetMinimum(vars[var]);//Minimum x in pdf
      float max = cut_tree->GetMaximum(vars[var]);//Maximum x in pdf
      float range = max - min;
      float min_hist = min - sqrt(pow((2*range)/n_bins,2));// Minimum x value for histogram
      float max_hist = max + sqrt(pow((2*range)/n_bins,2));// Maximum x value for histogram
      pdf_arr[f][var] = new TH1F(Form("%s_%s",vars[var].Data(),pdf[f].Data()),vars[var],n_bins,min_hist,max_hist); // Creating pdf matrix
      //Filling the histogram will values from data tree
      std::vector<float> vals;
      for(int i=0; i<pdf_size; i++){
        cut_tree->GetEntry(i);
        pdf_arr[f][var]->Fill(cut_tree->GetLeaf(vars[var])->GetValue(0));
      }
      pdf_arr[f][var]->Scale(1/pdf_arr[f][var]->Integral("width"));//Normalising to create probability distribution functions
      pdf_arr[f][var]->Draw("HIST C");//Draw smooth curve
      pdf_arr[f][var]->GetXaxis()->SetTitle(Form("%s",vars[var].Data())); pdf_arr[f][var]->GetXaxis()->SetTitleSize(0.05);
      pdf_arr[f][var]->GetYaxis()->SetTitle("Probability"); pdf_arr[f][var]->GetYaxis()->SetTitleSize(0.05);
      pdf_arr[f][var]->SetLineWidth(2);
      pdf_arr[f][var]->SetTitle("");
      pdf_arr[f][var]->SetStats(0);
      pdf_arr[f][var]->GetYaxis()->SetTitleOffset(1.2);//Ensure pdf y axis fits on canvas
      c1->SetLeftMargin(0.15);
      c1->SaveAs(Form("%s/pdf_images/%s_%s.png",pwd,vars[var].Data(),pdf[f].Data()));//Saving pdf image into folder within pwd
      delete c1; //delete canvas object for loop
    }
    delete tree;
    delete cut_tree;
  }

  //Size array before likelihood analysis i.e. Monte Carlo size
  float size_mc[data_files.size()];

  //Ratio vector array 
  std::vector<double> ratios[data_files.size()];

  for(int f=0; f<data_files.size(); f++){
    TFile * file = new TFile(data_files[f]);//Load files from header
    TTree * tree = (TTree*)file->Get("data");//Get the data tree
    gROOT->cd();
    TTree * cut_tree = tree->CopyTree("n100>0 && x>-9999 && dt_prev_us>0 && dt_prev_us<2000 && closestPMT>0");//Make cuts to remove unphysical parameters

    int size_cut = cut_tree->GetEntries();//Size of cut data tree

    //Number of events data tree
    TTree * run = (TTree*)file->Get("runSummary");
    std::vector<double> events;//Accounting for smaller MC files that make total file input
    //Original MC file size
    for(int i=0; i<run->GetEntries(); i++){
     run->GetEntry(i);
     events.push_back (run->GetLeaf("nEvents")->GetValue(0));
    }

    size_mc[f] = accumulate(events.begin(),events.end(),0.0);//Size of file f

    //Defining local signal and background likelihood arrays
    std::vector<double> ratio_like[vars.size()];

    for(int var=0; var<vars.size(); var++){
      for(int j=0; j<size_cut; j++){
        cut_tree->GetEntry(j);
        double val = cut_tree->GetLeaf(vars[var])->GetValue(0);
        //Finding bin values in the signal and background pdfs
        int x_sig = pdf_arr[0][var]->FindBin(val);
        int x_back = pdf_arr[1][var]->FindBin(val);
        //Finding bin contents of signal and background pdfs
        double y_sig = pdf_arr[0][var]->GetBinContent(x_sig);
        double y_back = pdf_arr[1][var]->GetBinContent(x_back);
        //Adding probabilities to likelihood matrix
        if(y_sig ==0 && y_back != 0){
          ratio_like[var].push_back (log(y_back));
        }
        else if(y_sig !=0 && y_back == 0){
          ratio_like[var].push_back (log(y_sig));
        }
        else if(y_sig != 0 && y_back != 0){
          ratio_like[var].push_back (log(y_sig) - log(y_back));
        }
        else{}//Filter zero terms
      }
    }

    for(int i=0; i<size_cut; i++){
      double ratio_total = 0;
      for(int var=0; var<vars.size(); var++){
        ratio_total += ratio_like[var][i];
      }
      ratios[f].push_back(ratio_total);
    }
  } 

  gROOT->SetBatch(0);// Enable ROOT canvas output

  //Stacked histogram, legend and title for plot
  THStack * stack = new THStack("stack","");
  TLegend *legend = new TLegend(0.6,0.6,0.9,0.9);
  TPaveText * t1 =new TPaveText(0.15,0.91,0.85,0.98,"NBNDC");
  t1->AddText("Likelihood Ratio Analysis");
  t1->SetFillStyle(0);
  t1->SetLineColor(0);

  TH1F * ratio_hist[7];// Histogram array for components
  std::vector<TString> names = {"big", "small","singles","world","geo","N","Li"};//Names for histograms

  //Filling and formatting histograms
  for(int f=0; f<data_files.size(); f++){
    ratio_hist[f] = new TH1F(Form("%s_hist",names[f].Data()), names[f].Data(),200,-80,100);
    for(int i=0; i<ratios[f].size();i++){
      ratio_hist[f]->Fill(ratios[f][i]);//Fill histogram with ratio values
    }
    ratio_hist[f]->SetLineWidth(2);
    ratio_hist[f]->SetLineColor(colours[f]);
    ratio_hist[f]->Scale(1/ratio_hist[f]->Integral("width"));//Normalisation
    stack->Add(ratio_hist[f]);//Adding individual component to stacked hist
    legend->AddEntry(ratio_hist[f],legend_names[f]);
  }

  //Plotting each component
  TCanvas * c1 = new TCanvas("c1");
  stack->Draw("NOSTACK HIST");
  stack->GetXaxis()->SetTitle("Likelihood Ratio"); stack->GetXaxis()->SetTitleSize(0.05);
  stack->GetYaxis()->SetTitle("Normalised Number of Events"); stack->GetYaxis()->SetTitleSize(0.05);
  legend->Draw();
  t1->Draw();

  

  //Output values to a csv file 
  std::ofstream ratio_csv;
  ratio_csv.open ("ratio_threshold.csv");//Open file to appened too
  ratio_csv << "Likelihood Ratio" << "," << "Big" << "," << "Small" << "," << "Combined" << "," << "World" << "," << "Geoneutrinos" << "," << "17 N" << "," << "9 Li" << "," << "Sigma DS1" << "," << "Sigma DS2" << "\n";

  for(int r=0; r<30; r++){
    std::vector<double> ratio_cut[7];//Create vector array for all components
    double rates_arr[7];//Create array for daily rates
    for(int f=0; f<data_files.size(); f++){
      for(int i=0; i<ratios[f].size(); i++){
        if(ratios[f][i] > r){ratio_cut[f].push_back (ratios[f][i]);}else{}//Applying threshold cut
      }
    }
    //Daily rates and 30-day significances
    for(int f=0; f<data_files.size(); f++){
      rates_arr[f] = rate(max_rate[f],ratio_cut[f].size(),size_mc[f]);
    }

    double signal_1 = rates_arr[0] + rates_arr[1];
    double signal_2 = rates_arr[0];

    double background_1 = rates_arr[2] + rates_arr[3] + rates_arr[4] + rates_arr[5] + rates_arr[6];
    double background_2 = rates_arr[1] + rates_arr[2] + rates_arr[3] + rates_arr[4] + rates_arr[5] + rates_arr[6];
    double s1 = sigma(signal_1,background_1);
    double s2 = sigma(signal_2,background_2);

    ratio_csv << r << "," << rates_arr[0] << "," << rates_arr[1]  << "," << rates_arr[2]  << "," <<  rates_arr[3] << "," << rates_arr[4]  << "," << rates_arr[5]  << "," << rates_arr[6] << "," << s1  << "," << s2 <<"\n";

    //Clear the vector array for the next likelihood ratio threshold cut
    for(int f=0; f<data_files.size(); f++){
      ratio_cut[f].clear();
    }
  }

  ratio_csv.close();//close opened csv file

}