#include "TROOT.h"
#include "TList.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDatime.h"
#include "TMath.h"
#include <string.h>
#include <math.h>

TGraphErrors* tg[2][5];

Int_t read_txt(char* folder, Int_t chan, Int_t ipnt, Double_t *wls, Double_t *sigs, Double_t *darks, Double_t *refs, double* sigs_cor)
{
  //chan : 0 or 1
  //ipnt: 0 to 4
  char fnam[256];
  sprintf(fnam,"%s/p%d_%7dU2.TXT", folder, ipnt+1, 1307077+chan);

  FILE* f=fopen(fnam,"r");
  if(!f) {
    printf("read_txt: file %s does not exist\n",fnam);
    return -1;
  }

  Int_t nwl=0;
  char lin[256];

  while(lin==fgets(lin,255,f)){
    for(Int_t i=0; i<strlen(lin); ++i) { if(';'==lin[i]) lin[i]=' '; }
    Double_t wl=0, sig=0, dark=0, ref=0, sig_cor=0;
    Int_t nit=sscanf(lin,"%lf %lf %lf %lf %lf",&wl,&sig,&dark,&ref,&sig_cor);
    if(5==nit){
      wls[nwl]=wl;
      sigs[nwl]=sig;
      darks[nwl]=dark;
      refs[nwl]=ref;
      sigs_cor[nwl]=sig_cor;
      nwl++;
    }
  }

  fclose(f);

  //printf("read_txt: file %s read with %d points\n",fnam, nwl);
  return nwl;
}

void draw2(TPad* lPad, TGraphErrors* tg1, TGraphErrors* tg2, char* tit=0){
  lPad->SetFillColor(0);
  lPad->Clear();

  lPad->SetTitle("Spectra");
  lPad->SetBottomMargin(0.15);
  lPad->SetRightMargin(0.05);
  lPad->SetLeftMargin(0.12);
  lPad->SetTopMargin(0.08);
  lPad->SetGrid();

  tg1->SetLineColor(4);
  tg2->SetLineColor(2);
  tg1->SetMarkerColor(4);
  tg2->SetMarkerColor(2);
  tg1->SetMarkerStyle(29);
  tg2->SetMarkerStyle(29);
  tg1->SetMarkerSize(0.5);
  tg2->SetMarkerSize(0.5);

  tg1->Draw("ALP");

  if(tit)tg1->GetHistogram()->SetTitle(tit);
  tg1->GetXaxis()->SetTitle("wavelength, nm");
  tg1->GetXaxis()->SetNdivisions(515);
  tg1->GetXaxis()->SetLabelSize(0.04);
  tg1->GetXaxis()->SetLabelFont(42);
  tg1->GetXaxis()->SetLabelOffset(0.02);
  tg1->GetXaxis()->SetTitleSize(0.05);
  tg1->GetXaxis()->SetTitleOffset(1.1);
  tg1->GetXaxis()->SetTitleFont(42);

  tg1->GetYaxis()->SetTitle("intensity");
  tg1->GetYaxis()->SetNdivisions(515);
  tg1->GetYaxis()->SetLabelSize(0.04);
  tg1->GetYaxis()->SetLabelOffset(0.02);
  tg1->GetYaxis()->SetLabelFont(42);
  tg1->GetYaxis()->SetTitleSize(0.05);
  tg1->GetYaxis()->SetTitleOffset(1.1);
  tg1->GetYaxis()->SetTitleFont(42);

  tg2->Draw("LP");

  Double_t xmin=1e9,xmax=-1e9,ymin=1e9,ymax=-1e9;

  TObject* obj; 
  TIter nex1(lPad->GetListOfPrimitives());
  while ((obj=nex1())){ 
    if(obj->IsA()==TGraphErrors::Class()){
      TGraphErrors* g=(TGraphErrors*)obj;
      Double_t gxmin,gxmax,gymin,gymax;
      g->ComputeRange( gxmin,gymin,gxmax,gymax);
      if(gxmin<xmin)xmin=gxmin;
      if(gxmax>xmax)xmax=gxmax;
      if(gymin<ymin)ymin=gymin;
      if(gymax>ymax)ymax=gymax;
    }
  }

  //printf("ranges: %f %f %f %f\n",xmin,xmax,ymin,ymax);
  Double_t dx=xmax-xmin;
  Double_t dy=ymax-ymin;

  TIter nex2(lPad->GetListOfPrimitives());
  while ((obj=nex2())){ 
    if(obj->IsA()==TGraphErrors::Class()){
      TGraphErrors* g=(TGraphErrors*)obj;
      g->GetXaxis()->SetLimits(xmin-0.05*dx,xmax+0.05*dx);
      g->GetYaxis()->SetLimits(ymin-0.05*dy,ymax+0.05*dy);
      g->GetXaxis()->UnZoom();
      g->GetYaxis()->UnZoom();
    }
  }

  lPad->Modified();
  lPad->Update();
}

void attlen(Double_t &len, Double_t &slen,Double_t xwl=550, Double_t dxwl=2.5,  Double_t x1=0, Double_t x2=10){
  printf("Length: %f %f \n",len,slen);
  Double_t s0[2][5], s1[2][5], s2[2][5];
  memset(&s0[0][0],0,sizeof(s0));
  memset(&s1[0][0],0,sizeof(s1));
  memset(&s2[0][0],0,sizeof(s2));
  for(Int_t ich=0; ich<2; ++ich){
    for(Int_t ipnt=0; ipnt<5; ++ipnt){
      Int_t n=(tg[ich][ipnt])->GetN();
      Double_t *x=(tg[ich][ipnt])->GetX();
      Double_t *y=(tg[ich][ipnt])->GetY();
      Double_t *ey=(tg[ich][ipnt])->GetEY();
      for(int i=0; i<n; ++i){
        if(x[i]>=xwl-dxwl && x[i]<xwl+dxwl){
          s0[ich][ipnt]+=1;
          s1[ich][ipnt]+=y[i];
          //s2[ich][ipnt]+=y[i]*y[i];
          s2[ich][ipnt]+=ey[i]*ey[i];
        }
      }
      s1[ich][ipnt]/=s0[ich][ipnt];
      //s2[ich][ipnt]=TMath::Sqrt(TMath::Abs(s2[ich][ipnt]/s0[ich][ipnt] - s1[ich][ipnt]*s1[ich][ipnt]));
      s2[ich][ipnt]=TMath::Sqrt(TMath::Abs(s2[ich][ipnt]/s0[ich][ipnt]));
      //printf("ich=%d ipnt=%d s0=%lf s1=%lf s2=%lf\n", ich, ipnt, s0[ich][ipnt], s1[ich][ipnt], s2[ich][ipnt]);
    }
  }
  
  Double_t xs[5]={1,3,5,7,9};
  Double_t ex[5]={0.1,0.1,0.1,0.1,0.1};
  Double_t rr[5], er[5];
  for(int i=0; i<5; ++i){
    if(s0[0][i]>0 && s0[1][i]>0 && s1[0][i]>0 && s1[1][i]>0){
      Double_t a1=s1[0][i];
      Double_t a2=s1[1][i];
      rr[i]=a1/a2;
      Double_t e1=s2[0][i]/TMath::Sqrt(s0[0][i]);
      Double_t e2=s2[1][i]/TMath::Sqrt(s0[1][i]);
      er[i]=rr[i]*(e1/a1+e2/a2);
    }
  }

  TGraphErrors *tgr=new TGraphErrors(5,xs,rr,ex,er);

  TCanvas* cr=(TCanvas*)gROOT->FindObjectAny("cr");
  if(!cr)cr=new TCanvas("cr","att len",600,400);
  cr->cd();
  cr->Clear();
  cr->SetFillColor(0);

  tgr->Draw("ALP");

  cr->SetFillColor(0);
  cr->SetBottomMargin(0.15);
  cr->SetRightMargin(0.05);
  cr->SetLeftMargin(0.12);
  cr->SetTopMargin(0.08);
  cr->SetGrid();

  char tit[256];
  sprintf(tit,"att len, %5.1f+-%5.1f nm",xwl,dxwl);
  tgr->GetHistogram()->SetTitle(tit);

  tgr->SetLineColor(4);
  tgr->SetMarkerColor(4);
  tgr->SetMarkerStyle(29);
  tgr->SetMarkerSize(1.5);

  tgr->GetXaxis()->SetTitle("x, cm");
  tgr->GetXaxis()->SetNdivisions(505);
  tgr->GetXaxis()->SetLabelSize(0.04);
  tgr->GetXaxis()->SetLabelFont(42);
  tgr->GetXaxis()->SetLabelOffset(0.02);
  tgr->GetXaxis()->SetTitleSize(0.05);
  tgr->GetXaxis()->SetTitleOffset(1.0);
  tgr->GetXaxis()->SetTitleFont(42);

  tgr->GetYaxis()->SetTitle("left/right ratio");
  tgr->GetYaxis()->SetNdivisions(505);
  tgr->GetYaxis()->SetLabelSize(0.04);
  tgr->GetYaxis()->SetLabelOffset(0.02);
  tgr->GetYaxis()->SetLabelFont(42);
  tgr->GetYaxis()->SetTitleSize(0.05);
  tgr->GetYaxis()->SetTitleOffset(1.1);
  tgr->GetYaxis()->SetTitleFont(42);

  int ind1=-1, ind2=-1;
  for(Int_t i=0; i<4; ++i){
    if(ind1<0 && xs[i]>=x1)ind1=i;
    if(xs[i]<=x2)ind2=i;
  }

  char funcname[256];
  sprintf(funcname,"flen_%d",(int)xwl);
  TF1* tfr = new TF1(funcname,"([0]*exp(-x*2./[1]))",x1,x2);
  Int_t np=tgr->GetN();
  Double_t *xx=tgr->GetX();
  Double_t *yy=tgr->GetY();
  Double_t p1=2*(xx[ind2]-xx[ind1])/(TMath::Log(TMath::Abs(yy[ind1]))-TMath::Log(TMath::Abs(yy[ind2])));
  Double_t p0=yy[2]*TMath::Exp(5./TMath::Abs(p1));
  //printf("initial: ind1=%d, ind2=%d, p0=%f, p1=%f\n",ind1,ind2,p0,p1);
  tfr->SetParameter(0,p0);
  tfr->SetParameter(1,TMath::Abs(p1));
  tfr->SetBit(kCanDelete);

  tfr->SetLineColor(2);

  //TFitResultPtr r=tgr->Fit(funcname,"S");
  TFitResultPtr r=tgr->Fit(funcname,"QSR");
  Double_t len1=(r->GetParams())[1];
  Double_t err1=(r->GetErrors())[1];
  //printf("For %f+- %f nm att len = %f +- %f\n", xwl, dxwl, len1, err1);

  //tf->Draw("same");
  cr->Modified();
  cr->Update();
  
  len=len1;
  slen=err1;
}

void transm(Double_t xwl=550, Double_t dxwl=2.5, Double_t att=0, Double_t satt=0, Double_t x1=0, Double_t x2=10){
  Double_t s0[2][5], s1[2][5], s2[2][5];
  memset(&s0[0][0],0,sizeof(s0));
  memset(&s1[0][0],0,sizeof(s1));
  memset(&s2[0][0],0,sizeof(s2));
  for(Int_t ich=0; ich<2; ++ich){
    for(Int_t ipnt=0; ipnt<5; ++ipnt){
      Int_t n=(tg[ich][ipnt])->GetN();
      Double_t *x=(tg[ich][ipnt])->GetX();
      Double_t *y=(tg[ich][ipnt])->GetY();
      Double_t *ey=(tg[ich][ipnt])->GetEY();
      for(int i=0; i<n; ++i){
        if(x[i]>=xwl-dxwl && x[i]<xwl+dxwl){
          s0[ich][ipnt]+=1;
          s1[ich][ipnt]+=y[i];
          //s2[ich][ipnt]+=y[i]*y[i];
          s2[ich][ipnt]+=ey[i]*ey[i];
        }
      }
      s1[ich][ipnt]/=s0[ich][ipnt];
      //s2[ich][ipnt]=TMath::Sqrt(TMath::Abs(s2[ich][ipnt]/s0[ich][ipnt] - s1[ich][ipnt]*s1[ich][ipnt]));
      s2[ich][ipnt]=TMath::Sqrt(TMath::Abs(s2[ich][ipnt]/s0[ich][ipnt]));
      //printf("ich=%d ipnt=%d s0=%lf s1=%lf s2=%lf\n", ich, ipnt, s0[ich][ipnt], s1[ich][ipnt], s2[ich][ipnt]);
    }
  }
  
  Double_t xs[5]={1,3,5,7,9};
  Double_t ex[5]={0.1,0.1,0.1,0.1,0.1};
  Double_t rr[5], er[5];
  for(int i=0; i<5; ++i){
    if(s0[0][i]>0 && s0[1][i]>0 && s1[0][i]>0 && s1[1][i]>0){
      Double_t a1=s1[0][i];
      Double_t a2=s1[1][i];
      rr[i]=a1/a2;
      Double_t e1=s2[0][i]/TMath::Sqrt(s0[0][i]);
      Double_t e2=s2[1][i]/TMath::Sqrt(s0[1][i]);
      er[i]=rr[i]*(e1/a1+e2/a2);
    }
  }

  TGraphErrors *tgr=new TGraphErrors(5,xs,rr,ex,er);

  TCanvas* cr=(TCanvas*)gROOT->FindObjectAny("cr");
  if(!cr)cr=new TCanvas("cr","att len",600,400);
  cr->cd();
  cr->Clear();
  cr->SetFillColor(0);

  tgr->Draw("ALP");

  cr->SetFillColor(0);
  cr->SetBottomMargin(0.15);
  cr->SetRightMargin(0.05);
  cr->SetLeftMargin(0.12);
  cr->SetTopMargin(0.08);
  cr->SetGrid();

  char tit[256];
  sprintf(tit,"att len, %5.1f+-%5.1f nm",xwl,dxwl);
  tgr->GetHistogram()->SetTitle(tit);

  tgr->SetLineColor(3);
  tgr->SetMarkerColor(3);
  tgr->SetMarkerStyle(29);
  tgr->SetMarkerSize(1.5);

  tgr->GetXaxis()->SetTitle("x, cm");
  tgr->GetXaxis()->SetNdivisions(505);
  tgr->GetXaxis()->SetLabelSize(0.04);
  tgr->GetXaxis()->SetLabelFont(42);
  tgr->GetXaxis()->SetLabelOffset(0.02);
  tgr->GetXaxis()->SetTitleSize(0.05);
  tgr->GetXaxis()->SetTitleOffset(1.0);
  tgr->GetXaxis()->SetTitleFont(42);

  tgr->GetYaxis()->SetTitle("left/right ratio");
  tgr->GetYaxis()->SetNdivisions(505);
  tgr->GetYaxis()->SetLabelSize(0.04);
  tgr->GetYaxis()->SetLabelOffset(0.02);
  tgr->GetYaxis()->SetLabelFont(42);
  tgr->GetYaxis()->SetTitleSize(0.05);
  tgr->GetYaxis()->SetTitleOffset(1.1);
  tgr->GetYaxis()->SetTitleFont(42);

  char funcname[256];
  sprintf(funcname,"fatt_%d",(int)xwl);
  TF1* tfr = new TF1(funcname,"([0]*pow([1],x/5.))",x1,x2);
  Int_t np=tgr->GetN();
  Double_t *xx=tgr->GetX();
  Double_t *yy=tgr->GetY();
  Double_t p0=yy[0];
  Double_t p1=yy[np-1]/yy[0];
  tfr->SetParameter(0,p0);
  tfr->SetParameter(1,TMath::Abs(p1));
  tfr->SetBit(kCanDelete);

  tfr->SetLineColor(2);

  //TFitResultPtr r=tgr->Fit(funcname,"S");
  TFitResultPtr r=tgr->Fit(funcname,"QSR");
  Double_t att1=(r->GetParams())[1];
  Double_t err1=(r->GetErrors())[1];
  //printf("For %f+- %f nm transmission = %f +- %f\n", xwl, dxwl, att1, err1);
  
  cr->Modified();
  cr->Update();
  
  att=att1;
  satt=err1;
}

void getcorr(Double_t* coef, Double_t xwl=650, Double_t dxwl=50, Double_t len=800){
  Double_t s0[2][5], s1[2][5], s2[2][5];
  memset(&s0[0][0],0,sizeof(s0));
  memset(&s1[0][0],0,sizeof(s1));
  memset(&s2[0][0],0,sizeof(s2));
  for(Int_t ich=0; ich<2; ++ich){
    for(Int_t ipnt=0; ipnt<5; ++ipnt){
      Int_t n=(tg[ich][ipnt])->GetN();
      Double_t *x=(tg[ich][ipnt])->GetX();
      Double_t *y=(tg[ich][ipnt])->GetY();
      Double_t *ey=(tg[ich][ipnt])->GetEY();
      for(int i=0; i<n; ++i){
        if(x[i]>=xwl-dxwl && x[i]<xwl+dxwl){
          s0[ich][ipnt]+=1;
          s1[ich][ipnt]+=y[i];
          //s2[ich][ipnt]+=y[i]*y[i];
          s2[ich][ipnt]+=ey[i]*ey[i];
        }
      }
      s1[ich][ipnt]/=s0[ich][ipnt];
      //s2[ich][ipnt]=TMath::Sqrt(TMath::Abs(s2[ich][ipnt]/s0[ich][ipnt] - s1[ich][ipnt]*s1[ich][ipnt]));
      s2[ich][ipnt]=TMath::Sqrt(TMath::Abs(s2[ich][ipnt]/s0[ich][ipnt]));
      //printf("ich=%d ipnt=%d s0=%lf s1=%lf s2=%lf\n", ich, ipnt, s0[ich][ipnt], s1[ich][ipnt], s2[ich][ipnt]);
    }
  }

  Double_t xs[5]={1,3,5,7,9};
  Double_t rr[5];
  for(int i=0; i<5; ++i){
    rr[i]=TMath::Exp(-2.*(xs[i]-xs[2])/len);
    if(s0[0][i]>0 && s0[1][i]>0 && s1[0][i]>0 && s1[1][i]>0){
      coef[i]=s1[0][i]/s1[1][i]/rr[i];
    }
  }
}

// folder - measurement directory, containing .txt files
// c1, c2 - coefficients (relevant for display of spectra only)
// wl0 - starting wavelength
// dwl - 1/2 of wavelength step
void spectra(char* folder, Double_t c1=1, Double_t c2=1, Double_t wl0=503, Double_t dwl=2.5)
{
  Double_t koef[2]; koef[0]=c1; koef[1]=c2;

  Int_t nwl=0;
  Double_t wls[2048],sigs[2048],darks[2048],refs[2048], sigs_cor[2048];
  Double_t ex[2048], ey[2048];  

  for(Int_t ich=0; ich<2; ++ich){
    for(Int_t ipnt=0; ipnt<5; ++ipnt){
      nwl=read_txt(folder,ich,ipnt,wls,sigs,darks,refs,sigs_cor);
      for(Int_t i=0; i<nwl; ++i){
        sigs_cor[i]*=koef[ich];
        //if(0==ich && 4==ipnt)sigs_cor[i]*=0.92;
        ex[i]=0.1; 
        ey[i]=10*koef[ich];
      }
      tg[ich][ipnt]=new TGraphErrors(nwl,wls,sigs_cor,ex,ey);
    }
  }

  // determine "average" attenuation length for 600...700 nm using points ##2,3,4 (without ##1,5)
  // attlen(650, 50, l_av, sl_av, 2.5, 7.5);
  // wavelength range: 650 +- 50
  // x coordinate range from 2.5 to 7.5, i.e. only points ##2,3,4
  //
  Double_t l_av=0, sl_av=0;
  attlen(l_av, sl_av,650, 50, 2.5, 7.5);
  printf("Average att len = %f +- %f\n",l_av,sl_av);
  if(l_av>1000 || l_av<1){
    printf("weird average att length, setting 1000 cm\n");
    l_av=1000; sl_av=1000;
  }

  // determine correction coefficients for all the 5 points using same 600...700 nm range 
  // and average attenuation length determined above
  Double_t coef[5];
  getcorr(coef,650,50,l_av);
  printf("Corrections: %f %f %f %f %f\n",coef[0],coef[1],coef[2],coef[3],coef[4]);

  // apply corrections
  for(Int_t ipnt=0; ipnt<5; ++ipnt){
    Int_t gn=(tg[1][ipnt])->GetN();
    Double_t *gy=(tg[1][ipnt])->GetY();
    Double_t *gey=(tg[1][ipnt])->GetEY();
    for(Int_t i=0; i<gn; ++i){
      gy[i]*=coef[ipnt];
      gey[i]*=coef[ipnt];
    }
  }

  TCanvas* cs=(TCanvas*)gROOT->FindObjectAny("cs");
  if(!cs)cs=new TCanvas("cs","spectra",900,600);
  cs->Clear();
  cs->SetFillColor(0);
  cs->Divide(3,2);

  for(Int_t ipnt=0; ipnt<5; ++ipnt){
    TPad* p=(TPad*)cs->cd(ipnt+1);
    char titl[256];
    sprintf(titl,"p%d, x=%d cm",ipnt+1,ipnt*2+1);
    draw2(p, tg[0][ipnt], tg[1][ipnt], titl);
  }

  Int_t np=41;

  for(Int_t iwl=0; iwl<np; ++iwl){
    wls[iwl]=wl0+iwl*2*dwl;
    ex[iwl]=dwl/1.7;
    attlen(sigs[iwl], ey[iwl], wls[iwl], dwl );
    if(sigs[iwl]<0){
      sigs[iwl]=5000;
      ey[iwl]=1;
    }
    else if(sigs[iwl]>5000){
      sigs[iwl]=5000;
      ey[iwl]=2500;
    }
    else if(ey[iwl]>5000){
      ey[iwl]=5000;
    }
  }

  TGraphErrors *tgattlen=new TGraphErrors(np,wls,sigs,ex,ey);
  TCanvas* cattlen=(TCanvas*)gROOT->FindObjectAny("cattlen");
  if(!cattlen)cattlen=new TCanvas("cattlen","attenuation length (wavelength)",600,400);
  cattlen->cd();
  cattlen->SetBottomMargin(0.15);
  cattlen->SetRightMargin(0.05);
  cattlen->SetLeftMargin(0.12);
  cattlen->SetTopMargin(0.08);
  cattlen->Clear();
  cattlen->SetFillColor(0);
  cattlen->SetGrid();

  tgattlen->Draw("ALP");
  tgattlen->GetHistogram()->SetTitle("attenuation length");

  tgattlen->SetLineColor(4);
  tgattlen->SetMarkerColor(4);
  tgattlen->SetMarkerStyle(29);
  //tgattlen->SetMarkerSize(0.5);

  tgattlen->GetXaxis()->SetTitle("wavelength, nm");
  tgattlen->GetXaxis()->SetNdivisions(515);
  tgattlen->GetXaxis()->SetLabelSize(0.04);
  tgattlen->GetXaxis()->SetLabelFont(42);
  tgattlen->GetXaxis()->SetLabelOffset(0.02);
  tgattlen->GetXaxis()->SetTitleSize(0.05);
  tgattlen->GetXaxis()->SetTitleOffset(1.0);
  tgattlen->GetXaxis()->SetTitleFont(42);

  tgattlen->GetYaxis()->SetTitle("att length, cm");
  tgattlen->GetYaxis()->SetNdivisions(515);
  tgattlen->GetYaxis()->SetLabelSize(0.04);
  tgattlen->GetYaxis()->SetLabelOffset(0.02);
  tgattlen->GetYaxis()->SetLabelFont(42);
  tgattlen->GetYaxis()->SetTitleSize(0.05);
  tgattlen->GetYaxis()->SetTitleOffset(1.1);
  tgattlen->GetYaxis()->SetTitleFont(42);

  for(Int_t iwl=0; iwl<np; ++iwl){
    wls[iwl]=wl0+iwl*2*dwl;
    ex[iwl]=dwl/7;
    transm(wls[iwl], dwl, sigs[iwl], ey[iwl]);
  }

  TGraphErrors *tgtransm=new TGraphErrors(np,wls,sigs,ex,ey);
  TCanvas* ctransm=(TCanvas*)gROOT->FindObjectAny("ctransm");
  if(!ctransm)ctransm=new TCanvas("ctransm","transmittance (wavelength)",600,400);
  ctransm->cd();
  ctransm->SetBottomMargin(0.15);
  ctransm->SetRightMargin(0.05);
  ctransm->SetLeftMargin(0.12);
  ctransm->SetTopMargin(0.08);
  ctransm->Clear();
  ctransm->SetFillColor(0);
  ctransm->SetGrid();

  tgtransm->Draw("ALP");
  tgtransm->GetHistogram()->SetTitle("transmittance of a 10 cm fiber");

  tgtransm->SetLineColor(4);
  tgtransm->SetMarkerColor(2);
  tgtransm->SetMarkerStyle(29);
  //tgtransm->SetMarkerSize(0.5);

  tgtransm->GetXaxis()->SetTitle("wavelength, nm");
  tgtransm->GetXaxis()->SetNdivisions(515);
  tgtransm->GetXaxis()->SetLabelSize(0.04);
  tgtransm->GetXaxis()->SetLabelFont(42);
  tgtransm->GetXaxis()->SetLabelOffset(0.02);
  tgtransm->GetXaxis()->SetTitleSize(0.05);
  tgtransm->GetXaxis()->SetTitleOffset(1.0);
  tgtransm->GetXaxis()->SetTitleFont(42);

  tgtransm->GetYaxis()->SetTitle("transmittance");
  tgtransm->GetYaxis()->SetNdivisions(515);
  tgtransm->GetYaxis()->SetLabelSize(0.04);
  tgtransm->GetYaxis()->SetLabelOffset(0.02);
  tgtransm->GetYaxis()->SetLabelFont(42);
  tgtransm->GetYaxis()->SetTitleSize(0.05);
  tgtransm->GetYaxis()->SetTitleOffset(1.0);
  tgtransm->GetYaxis()->SetTitleFont(42);

}
