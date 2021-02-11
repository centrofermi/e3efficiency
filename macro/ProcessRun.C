void GetEff(TH1F *h1,TH1F *h2,Float_t &eff,Float_t &err);
void GetEffFit(TH1F *h1,TH1F *h2,Float_t &eff,Float_t &err);

// Parameters
Float_t minZDir = 0.9; // To select vertical tracks
Float_t betaWindow = 0.7; // how large beta is cut
Float_t betapeak = 1.05; // where inverse beta is suppose to be
Float_t betaBackground = 3; // where the residual background is estimated (syst error)
Float_t residualcut = 7; // to decide if matching is accepted (in cm)
TF1 *fgaus2;

void ProcessRun(const char *name="BOLO-01-2017-09-08-16018.root",Int_t ch=0,Int_t hv=16018){

  Float_t timeshift = 0.5;
  const Int_t nbinsX = 45;
  const Int_t nbinsY = 30;
  const Int_t nbinsZ = 6;

  TFile *f = new TFile(name);

  TTree *t = (TTree *) f->Get(Form("EventsCal%i",ch));
  TTree *th = (TTree *) f->Get("Header");

  th->GetEvent(0);
  Float_t distTM = th->GetLeaf("Ztop")->GetValue() - th->GetLeaf("Zmiddle")->GetValue();
  Float_t distMB = th->GetLeaf("Zmiddle")->GetValue() - th->GetLeaf("Zbottom")->GetValue();

  Int_t n = t->GetEntries();

  Float_t x,y,zdir,betainv;

  Int_t nden=0,nnum=0,nnumsyst=0,ndensyst=0,match;

  Float_t eff,erreff,errsyst;

  TH1F *hbeta = new TH1F("hbeta",";#beta_inv",200,-1,4);
  TH1F *hbetaMatch = new TH1F("hbetaMatch",";#beta_inv",200,-1,4);
  TH1F *hbetaNoMatch = new TH1F("hbetaNoMatch",";#beta_inv",200,-1,4);
  TProfile *hEffX = new TProfile("hEffX",";X (cm);#varepsilon",90,-20,160);
  TProfile *hEffY = new TProfile("hEffY",";Y (cm);#varepsilon",60,-10,86);


  TH1F *hbetaMatchPosX[nbinsX];
  TH1F *hbetaNoMatchPosX[nbinsX];
  TH1F *hEffXf = new TH1F("hEffXf","",nbinsX,-20,160);
  for(Int_t i=0;i<nbinsX;i++){
    hbetaMatchPosX[i] = new TH1F(Form("hbetaMatchX%i",i),";#beta_inv",200,-1,4);
    hbetaNoMatchPosX[i] = new TH1F(Form("hbetaNoMatchX%i",i),";#beta_inv",200,-1,4);
  }
  
  TH1F *hbetaMatchPosY[nbinsY];
  TH1F *hbetaNoMatchPosY[nbinsY];
  TH1F *hEffYf = new TH1F("hEffYf","",nbinsY,-10,86);
  for(Int_t i=0;i<nbinsY;i++){
    hbetaMatchPosY[i] = new TH1F(Form("hbetaMatchY%i",i),";#beta_inv",200,-1,4);
    hbetaNoMatchPosY[i] = new TH1F(Form("hbetaNoMatchY%i",i),";#beta_inv",200,-1,4);
  }
  
  TH1F *hbetaMatchPosZ[nbinsZ];
  TH1F *hbetaNoMatchPosZ[nbinsZ];
  TH1F *hEffZf = new TH1F("hEffZf","",nbinsZ,0.7,1);
  for(Int_t i=0;i < nbinsZ;i++){
    hbetaMatchPosZ[i] = new TH1F(Form("hbetaMatchZ%i",i),";#beta_inv",200,-1,4);
    hbetaNoMatchPosZ[i] = new TH1F(Form("hbetaNoMatchZ%i",i),";#beta_inv",200,-1,4);
  }
  
  TH1F *hXdiff = new TH1F("hXdiff","Match in middle chamber;X_{inter} - X_{hit} (cm)",100,-20,20);
  TH1F *hYdiff = new TH1F("hYdiff","Match in middle chamber;Y_{inter} - Y_{hit} (cm)",100,-20,20);
  
  TH1F *htof = new TH1F("htof","Time shift in time of flight;#Deltat (ns)",400,-10,10);
  TF1 *fgaus = new TF1("fgaus","gaus(0)+[3]",-10,10);
  fgaus->SetParameter(0,100);
  fgaus->SetParameter(1,0);
  fgaus->SetParameter(2,1);
  fgaus->SetParameter(3,0);

  fgaus2 = new TF1("fgaus2","gaus(0)+[3]",-1,4);
  fgaus2->SetParameter(0,100);
  //  fgaus2->SetParLimits(0,0,10000000);
  fgaus2->SetParameter(1,betapeak);
  fgaus2->SetParLimits(1,betapeak*0.9,betapeak*1.1);
  fgaus2->SetParameter(2,1);
  fgaus2->SetParLimits(2,0.05,0.5);
  fgaus2->SetParameter(3,0);

  for(Int_t i=0;i<n;i++){
    t->GetEvent(i);

    if(t->GetLeaf("Ntracks")->GetValue() == 0) continue;

    htof->Fill(t->GetLeaf("TrackLength")->GetValue()/30*betapeak - t->GetLeaf("TimeOfFlight")->GetValue()); 

  }

  timeshift = htof->GetXaxis()->GetBinCenter(htof->GetMaximumBin());
  fgaus->SetParLimits(0,1,1000000);
  fgaus->SetParameter(1,timeshift);
  fgaus->SetParLimits(1,timeshift-1,timeshift+1);
  fgaus->SetParLimits(2,0.2,2);

  htof->Fit(fgaus);
  timeshift = fgaus->GetParameter(1);

  printf("time shift = %f +/- %f\n",fgaus->GetParameter(1),fgaus->GetParError(1));

      Int_t iX,iY,iZ;

  for(Int_t i=0;i<n;i++){
    t->GetEvent(i);

    if(t->GetLeaf("Ntracks")->GetValue() == 0) continue;
    
    if(ch==0){
      x = t->GetLeaf("PosXBot")->GetValue();
      y = t->GetLeaf("PosYBot")->GetValue();
    }
    else if(ch==1){
      x = t->GetLeaf("PosXMid")->GetValue();
      y = t->GetLeaf("PosYMid")->GetValue();
    }
    else if(ch==2){
      x = t->GetLeaf("PosXTop")->GetValue();
      y = t->GetLeaf("PosYTop")->GetValue();
    }

    zdir = t->GetLeaf("ZDir")->GetValue();

    betainv = (t->GetLeaf("TimeOfFlight")->GetValue()+timeshift)/t->GetLeaf("TrackLength")->GetValue()*30;

    match = (t->GetLeaf("StripMatched")->GetValue() > -1);

    // if(TMath::Abs(t->GetLeaf("IntersectXMid")->GetValue()-t->GetLeaf("PosXMid")->GetValue()) > residualcut && match) match=0;
    // if(TMath::Abs(t->GetLeaf("IntersectYMid")->GetValue()-t->GetLeaf("PosYMid")->GetValue()) > residualcut && match) match=0;

    iX = (x+20)/180*nbinsX;
    iY = (y+10)/96*nbinsY;
    iZ = (zdir-0.7)/0.3*nbinsZ;

    if(iX < 0 || iX > nbinsX-1) continue;
    if(iY < 0 || iY > nbinsY-1) continue;
    if(iZ < 0 || iZ > nbinsZ-1) continue;

    if(match)
      hbetaMatchPosZ[iZ]->Fill(betainv);
    else
      hbetaNoMatchPosZ[iZ]->Fill(betainv);

    if(zdir < minZDir) continue; // keep only vertical tracks
 
    if(match){
      if(y > 5 && y < 75) hbetaMatchPosX[iX]->Fill(betainv);
      if(x > 40 && y < 100) hbetaMatchPosY[iY]->Fill(betainv);
    }
    else{
      if(y > 5 && y < 75) hbetaNoMatchPosX[iX]->Fill(betainv);
     if(x > 40 && y < 100)  hbetaNoMatchPosY[iY]->Fill(betainv);
    }

    if(y < 5 || y > 75) continue; 
    if(x < 40 || y > 100) continue; 
 
    hbeta->Fill(betainv);

    if(match){
      hbetaMatch->Fill(betainv);
      hXdiff->Fill(t->GetLeaf("IntersectXMid")->GetValue()-t->GetLeaf("PosXMid")->GetValue());
      hYdiff->Fill(t->GetLeaf("IntersectYMid")->GetValue()-t->GetLeaf("PosYMid")->GetValue());
    }
    else{
      hbetaNoMatch->Fill(betainv);
    }
    
    if(TMath::Abs(betainv - betaBackground) < betaWindow){
      ndensyst++;
      if(match) nnumsyst++;
    }

    if(TMath::Abs(betainv - betapeak) > betaWindow) continue; // Skip wrong beta
    
    nden++;
    
    if(match) nnum++;
    
    hEffX->Fill(x,match);
    hEffY->Fill(y,match);
    
  }
  
  nden -= ndensyst;
  nnum -= nnumsyst;

  hbeta->Fit(fgaus2,"");
  fgaus2->FixParameter(1,fgaus2->GetParameter(1));
  fgaus2->FixParameter(2,fgaus2->GetParameter(2));

  Float_t totalError;
  Float_t efffit;

  for(Int_t i=0;i < nbinsX;i++){
    GetEff(hbetaMatchPosX[i],hbetaNoMatchPosX[i],efffit,totalError);
    hEffXf->SetBinContent(i+1,efffit);
    hEffXf->SetBinError(i+1,totalError);
  }
  for(Int_t i=0;i < nbinsY;i++){
    GetEff(hbetaMatchPosY[i],hbetaNoMatchPosY[i],efffit,totalError);
    hEffYf->SetBinContent(i+1,efffit);
    hEffYf->SetBinError(i+1,totalError);
  }
	for(Int_t i=0;i < nbinsZ;i++){
    GetEff(hbetaMatchPosZ[i],hbetaNoMatchPosZ[i],efffit,totalError);
    hEffZf->SetBinContent(i+1,efffit);
    hEffZf->SetBinError(i+1,totalError);
  }


  eff = nnum*1.0/nden;
  erreff = TMath::Sqrt(eff*(1-eff)/nden);
  errsyst = TMath::Sqrt(nnumsyst*1.0)/nden;

  printf("count stats: den=%i num=%i syst=%i\n",nden+nnumsyst,nnum,nnumsyst);
  printf("syst also subtracted to the denominator -> eff = num/(den-syst)\n");

  printf("\nH.V   efficiency    stat_err     syst_err\n");
  printf("%i %f  +/- %f +/- %f\n\n",hv,eff,erreff,errsyst);
  system(Form("echo %i %f %f %f >>results_%i",hv,eff,erreff,errsyst,ch));

  GetEff(hbetaMatch,hbetaNoMatch,efffit,totalError);

  printf("efficiency by integral = %f +/- %f\n",efffit,totalError);

  system(Form("echo %i %f %f >>newresults_%i",hv,efffit,totalError,ch));

  TCanvas *c = new TCanvas();
  hbeta->Draw();
  hbetaMatch->Draw("SAME");
  hbetaNoMatch->Draw("SAME");
  hbetaMatch->SetLineColor(4);
  c->Update();
  TLine *l = new TLine(betapeak-betaWindow,0,betapeak-betaWindow,hbeta->GetMaximum());
  l->SetLineColor(2);
  l->SetLineStyle(2);
  l->SetLineWidth(2);
  l->Draw("SAME");
  TLine *l2 = new TLine(betapeak+betaWindow,0,betapeak+betaWindow,hbeta->GetMaximum());
  l2->SetLineColor(2);
  l2->SetLineStyle(2);
  l2->SetLineWidth(2);
  l2->Draw("SAME");

  TCanvas *c2 = new TCanvas();
  c2->Divide(2,1);
  c2->cd(1);
  hEffX->Draw();
  c2->cd(2);
  hEffY->Draw();

  TCanvas *c3 = new TCanvas();
  c3->Divide(2,1);
  c3->cd(1);
  hXdiff->Draw();
  TLine *l3 = new TLine(-residualcut,0,-residualcut,hXdiff->GetMaximum());
  l3->SetLineColor(2);
  l3->SetLineStyle(2);
  l3->SetLineWidth(2);
  l3->Draw("SAME");
  TLine *l4 = new TLine(residualcut,0,residualcut,hXdiff->GetMaximum());
  l4->SetLineColor(2);
  l4->SetLineStyle(2);
  l4->SetLineWidth(2);
  l4->Draw("SAME");
  c3->cd(2);
  hYdiff->Draw();
  TLine *l5 = new TLine(-residualcut,0,-residualcut,hYdiff->GetMaximum());
  l5->SetLineColor(2);
  l5->SetLineStyle(2);
  l5->SetLineWidth(2);
  l5->Draw("SAME");
  TLine *l6 = new TLine(residualcut,0,residualcut,hYdiff->GetMaximum());
  l6->SetLineColor(2);
  l6->SetLineStyle(2);
  l6->SetLineWidth(2);
  l6->Draw("SAME");

  c2->cd(1);
  hEffXf->Draw();
  c2->cd(2);
  hEffYf->Draw();
  TCanvas *c4 = new TCanvas();
  hEffZf->Draw();
  TLine *l7 = new TLine(minZDir,hEffZf->GetMinimum(),minZDir,hEffZf->GetMaximum());
  l7->Draw("SAME");
  l7->SetLineColor(2);
  l7->SetLineStyle(2);
  l7->SetLineWidth(2);
 
  TFile *fout = new TFile(Form("resultsC%i_%i.root",ch,hv),"RECREATE");
  c->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  fout->Close();
}

void GetEff(TH1F *h1,TH1F *h2,Float_t &eff,Float_t &err){
  Float_t betaWindow = 0.7; // how large beta is cut
  Float_t betapeak = 1.05; // where inverse beta is suppose to be
  Float_t betaBackground = 3; // where the residual background is estimated (syst error)
  
  Float_t num=0,den=0;

  for(Int_t i=1;i <= h1->GetNbinsX();i++){
    Float_t x=h1->GetBinCenter(i);
    if(TMath::Abs(x-betapeak)<betaWindow){
      num+=h1->GetBinContent(i);
      den+=h2->GetBinContent(i);
    }
    else if(TMath::Abs(x-betaBackground)<betaWindow){
      num-=h1->GetBinContent(i);
      den-=h2->GetBinContent(i);
    }
  }
  den += num;

  if(den < 50){
    eff=0;
    err=0;
    return;
  }

  eff = num/den;
  err=TMath::Sqrt(eff*(1-eff)/den);

}

void GetEffFit(TH1F *h1,TH1F *h2,Float_t &eff,Float_t &err){
  if(h1->GetEntries()  + h2->GetEntries() < 100){
    eff=0;
    err=0;
    return;
  }

  for(Int_t i=1;i <= h1->GetNbinsX();i++){
    // set empty bin error to 1
    if(h1->GetBinContent(i) < 1){
      h1->SetBinError(i,1);
    }
    if(h2->GetBinContent(i) < 1){
      h2->SetBinError(i,1);
    }
  }

  h2->Fit(fgaus2,"");
  fgaus2->SetParameter(3,0);
  Float_t intDen = fgaus2->Integral(-1,4)/h2->GetBinWidth(1);
  Float_t errDen = fgaus2->GetParError(0)/fgaus2->GetParameter(0)*intDen;
  h1->Fit(fgaus2,"");
  fgaus2->SetParameter(3,0);
  Float_t intNum = fgaus2->Integral(-1,4)/h1->GetBinWidth(1);
  Float_t errNum = fgaus2->GetParError(0)/fgaus2->GetParameter(0)*intNum;

  Float_t totalError = intDen*intDen*errNum;
  totalError += intNum*intNum*errDen;
  totalError = TMath::Sqrt(totalError);
  totalError /= (intNum+intDen)*(intNum+intDen);

  eff = intNum/(intDen+intNum);
  err = totalError;//TMath::Sqrt(eff*(1-eff)/(intDen+intNum));

  printf("Num from fit = %f -- Den from fit = %f\n",intNum,intNum+intDen);
}
