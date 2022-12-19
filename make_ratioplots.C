#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TLegendEntry.h"

#include <iostream>

//!! To compute SF = MC/Data
TH1F* DividTGraphs(TGraphAsymmErrors* gr1, TGraphAsymmErrors* gr2){

    int nbins = gr1->GetN();
    printf("%d",nbins);
    double xbins[nbins+1];

    for(int i = 0;  i < nbins; ++i){

        Double_t x = 999; 
        Double_t x_hi = 999; 
        Double_t x_low = 999; 
        Double_t y = 999; 
        gr1->GetPoint(i,x,y);
        x_hi = gr1->GetErrorXhigh(i);
        x_low = gr1->GetErrorXlow(i);
        if(i == nbins-1){
            xbins[i] = x-x_low;
            xbins[i+1] = x+x_hi;
        }else{
            xbins[i] = x-x_low;
        }
    }

    TH1F *h1 = new TH1F("","",nbins,xbins);
    TH1F *h2 = new TH1F("","",nbins,xbins);

    TGraphAsymmErrors* gr[2] = {gr1, gr2};
    TH1F* h[2] = {h1, h2};

    //Loop over bins to do ratio
    //
    for (int k = 0; k < 2; ++k){
        for(int i = 0;  i < nbins+1; ++i){
            //
            //TGraph
            //
            Double_t num_x = 999; 
            Double_t num_y = 999; 
            Double_t num_y_hi = 999; 
            Double_t num_y_low = 999; 

            gr[k]->GetPoint(i,num_x,num_y);
            num_y_hi = gr[k]->GetErrorYhigh(i);
            num_y_low = gr[k]->GetErrorYlow(i);

            double max_error = max(num_y_hi,num_y_low);

            //Convert into TH1D
            h[k]->SetBinContent(h[k]->FindBin(num_x), num_y);
            h[k]->SetBinError(h[k]->FindBin(num_x), max_error);
        }
    }

    //ratio histogram
    h[0]->Divide(h[1]);

    return h[0]; 

}
