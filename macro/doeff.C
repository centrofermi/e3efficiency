doeff(const char *ff="BOLO-04-2016-05-26",Int_t tophv=17400,Int_t bottomhv=16280){
	FILE *f = fopen("runs","r");

	Int_t n=0;
	Int_t hv;
	
	TFile *fin,*fin2;
	TH1D *h,*h2,*h3,*hb,*ht,*hb2,*ht2;
	TF1 *fitf=new TF1("fitf","pol0",20,70);
	Float_t x[100],x2[100],ex[100],y[100],ey[100],y2[100],ey2[100];
        Float_t yb[100],eyb[100],yt[100],eyt[100];
	Float_t p,T;

	TF1 *fsig = new TF1("fsig","gaus",0.7,1.3);
        TF1 *fall = new TF1("fall","gaus(0)+pol2(3)",0.6,1.25);
	fsig->SetParameter(0,1000);
        fsig->SetParameter(1,0.93);
        fsig->SetParameter(2,0.05);


	while(fscanf(f,"%i %f %f",&hv,&p,&T)==3){
		fin = TFile::Open(Form("%s-%i_res.root",ff,hv));
		h = (TH1D *) fin->Get("hmEff");
                hh2 = (TH1D *) fin->Get("hmNum");
                hh3 = (TH1D *) fin->Get("hmDen");

                hnum = (TH1D *) fin->Get("hmNum3");
                hden = (TH1D *) fin->Get("hmDen3");
		hden->Add(hnum,-1);
                hnum->Fit(fsig,"WW","",0.7,1.3);
                fall->FixParameter(1,fsig->GetParameter(1));
                fall->FixParameter(2,fsig->GetParameter(2));
                hden->Fit(fall,"WW","",0.6,1.25);
 
                fin2 = TFile::Open(Form("%s-%i_check.root",ff,hv));
                h2 = (TH1D *) fin2->Get("hmNum");
                h3 = (TH1D *) fin2->Get("hmDen");
                hb = (TH1D *) fin2->Get("hbNum");
                hb2 = (TH1D *) fin2->Get("hbDen");
                ht = (TH1D *) fin2->Get("htNum");
                ht2 = (TH1D *) fin2->Get("htDen");
		x[n]=hv;
                x2[n] = hv*(1000/p)*((T+273)/(25+273));
		ex[n]=0;
                y[n]=fsig->GetParameter(0)/(fsig->GetParameter(0)+fall->GetParameter(0));//hh2->GetEntries()/hh3->GetEntries();//fitf->GetParameter(0);
                ey[n]=sqrt(hh2->GetEntries()*(1-y[n]))*y[n]/hh2->GetEntries();
//                h2->Fit(fitf);
                y2[n]=h2->GetEntries()/h3->GetEntries();//fitf->GetParameter(0);
                ey2[n]=sqrt(h2->GetEntries()*(1-y2[n]))/h3->GetEntries();

//                hb->Fit(fitf,"","",2,70);
                yb[n]=hb->GetEntries()/hb2->GetEntries();//fitf->GetParameter(0);
                eyb[n]=sqrt(hb->GetEntries()*(1-yb[n]))/hb2->GetEntries();
//                ht->Fit(fitf,"","",2,70);
                yt[n]=ht->GetEntries()/ht2->GetEntries();//fitf->GetParameter(0);
                cout << n << " " << x[n] << " " << ht->GetEntries() << " " << ht2->GetEntries() << " " << yt[n] << endl; 
                eyt[n]=sqrt(1./(ht->GetEntries()+ht2->GetEntries()))*yt[n];//fitf->GetParError(0);
                eyt[n]=sqrt(ht->GetEntries()*(1-yt[n]))/ht2->GetEntries();

		n++;
	}
	new TCanvas();
        TGraphErrors *gb = new TGraphErrors(n,x,yb,ex,eyb);

        gb->SetMarkerStyle(20);
        gb->Draw("AP");
	gb->SetTitle(Form("Bottom H.V.=%i;middle H.V.",bottomhv));
        gb->Fit(fitf,"","",15000,20000);
	Float_t xtemp[2],ytemp[2],eytemp[2];
	xtemp[0]=bottomhv;
	ytemp[0]=fitf->GetParameter(0);
        eytemp[0]=fitf->GetParError(0);
	TGraphErrors *gbf = new TGraphErrors(1,xtemp,ytemp,ex,eytemp);
	gbf->SetMarkerStyle(20);
        gbf->SetMarkerColor(2);
        gbf->SetLineColor(2);

        new TCanvas();
        TGraphErrors *gt = new TGraphErrors(n,x,yt,ex,eyt);

        gt->SetMarkerStyle(20);
        gt->Draw("AP");
        gt->SetTitle(Form("Top H.V.=%i;middle H.V.",tophv));
        gt->Fit(fitf,"","",15000,20000);
	xtemp[0]=tophv;
        ytemp[0]=fitf->GetParameter(0);
        eytemp[0]=fitf->GetParError(0);
        TGraphErrors *gtf = new TGraphErrors(1,xtemp,ytemp,ex,eytemp);
        gtf->SetMarkerStyle(20);
        gtf->SetMarkerColor(2);
        gtf->SetLineColor(2);


	new TCanvas();
	TGraphErrors *g = new TGraphErrors(n,x,y,ex,ey);

        g->SetMarkerStyle(20);
        g->Draw("AP");

	TGraphErrors *g2 = new TGraphErrors(n,x,y2,ex,ey2);

        TGraphErrors *gcorr = new TGraphErrors(n,x2,y,ex,ey);

        gcorr->SetMarkerStyle(20);
        gcorr->Draw("AP");


        g2->SetMarkerStyle(24);
        g2->Draw("P");


	gbf->Draw("P");
        gtf->Draw("P");

	TFile *fout = new TFile(Form("res%s.root",ff),"RECREATE");
        g->SetName("effMiddle");
	g->Write();
        gcorr->SetName("effMiddleCorr");
        gcorr->Write();
	g2->SetName("effMiddleCheck");
        g2->Write();
        gbf->SetName("effBottom");
        gbf->Write();
        gtf->SetName("effTop");
        gtf->Write();
	fout->Close();

	gt->Draw();

}
