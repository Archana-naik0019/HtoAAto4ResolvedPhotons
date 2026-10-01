void ConvertHistToGraph(TGraph* graph, TH1D* hist) {
    graph->Set(static_cast<int>(hist->GetNbinsX()));

    // Exclude bin 0 and bin N+1 (under/overflow)
    for (Int_t i = 0; i < hist->GetNbinsX(); i++) {
        graph->SetPoint(i, hist->GetBinCenter(i+1), hist->GetBinContent(i+1));
    }
}


//TGraph* SmoothHist(std::string name, TH1D* hist, Float_t smoothness, Float_t span) {
std::tuple<TH1D*, TGraph*> SmoothHist(std::string name, TH1D* hist, Float_t smoothness, Float_t span) {

    TGraph* graph = new TGraph();
    ConvertHistToGraph(graph, hist);

    TGraphSmooth* smoother = new TGraphSmooth();
    TGraph* smooth;
    if(span == -1.0) {
        smooth = smoother->SmoothSuper(graph, "", smoothness);
    }
    else {
        smooth = smoother->SmoothSuper(graph, "", smoothness, span);
    }
    smooth->SetName(name.c_str());
    smooth->SetTitle(name.c_str());

    // Manually convert TGraph to histogram. GetHistogram doesn't do anything: https://root-forum.cern.ch/t/tgraph-gethistogram/8911/2
    TH1D* smoothed_hist = new TH1D(*hist);
    for (int i=0; i < smoothed_hist->GetNbinsX(); ++i) {
        double x,y;
        smooth->GetPoint(i, x, y);
        if(i < smoothed_hist->GetNbinsX() - 2) {
            smoothed_hist->SetBinContent(i, smooth->Eval(x, 0, "S"));
        }
        else {
            double x2,y2;
            smooth->GetPoint(i-1, x2, y2);
            std::cout << y2 << std::endl;
            smoothed_hist->SetBinContent(i, (smooth->Eval(x, 0, "S") + y2) / 2.0);  // Acount for edge case where it uses spline wrong
        }
    }
    for (unsigned short i=1; i < smoothed_hist->GetNbinsX(); ++i) {
        smoothed_hist->SetBinError(i, 0);
    }

    delete graph;
    return std::make_tuple(smoothed_hist, smooth);
}





