#include <TParameter.h>
#include <TFile.h>
#include <string>
#include <TH1D.h>
#include <TGraph2D.h>
#include <vector>
#include <TH2D.h>
#include <TGraph.h>

struct Fxy
{
    double Fx;
    double Fy;
    double xpos;
    double ypos;
};

double gPitch = 300e-6; // 300 um
int gNSteps = 11;
double gXStep = 300e-6 / (gNSteps - 1);
double gYStep = 300e-6 / (gNSteps - 1);
std::pair<double, double> GetPosition(int x_index, int y_index)
{
    double x_pos = gPitch + x_index * gXStep;
    double y_pos = gPitch + y_index * gYStep;
    return std::make_pair(x_pos, y_pos);
}

void process(int x_loop = 11, int y_loop = 11)
{
    auto c1 = new TCanvas("c1", "c1", 800, 600);
    std::vector<Fxy> vFxy;
    for (int x = 0; x < x_loop; x++)
    {
        for (int y = 0; y < y_loop; y++)
        {
            // Perform some processing here
            std::string sFolder = "ana/x_" + std::to_string(x) + "_y_" + std::to_string(y) + "/";
            TFile *file = TFile::Open((sFolder + "ve_bipolar_analysis.root").c_str(), "READ");
            if (file && !file->IsZombie())
            {
                Fxy fxy;
                TParameter<double> *Fx = nullptr;
                file->GetObject("Fx", Fx);
                if (Fx)
                    fxy.Fx = Fx->GetVal();
                TParameter<double> *Fy = nullptr;
                file->GetObject("Fy", Fy);
                if (Fy)
                    fxy.Fy = Fy->GetVal();
                auto pos = GetPosition(x, y);
                fxy.xpos = pos.first;
                fxy.ypos = pos.second;

                std::cout << "x_index: " << x << ", y_index: " << y
                          << ", xpos(um): " << fxy.xpos * 1e6 - 300 << ", ypos(um): " << fxy.ypos * 1e6 - 300
                          << ", Fx: " << fxy.Fx << ", Fy: " << fxy.Fy << std::endl;

                file->Close();
                delete file;
                vFxy.push_back(fxy);
            }
        }
    }

    // Save results to a new ROOT file
    TFile *outFile = TFile::Open("Fxy_results.root", "RECREATE");

    // 作图：Fx与x\y的关系，Fy与x\y的关系
    TGraph2D *grFx = new TGraph2D();
    TGraph2D *grFy = new TGraph2D();
    for (size_t i = 0; i < vFxy.size(); i++)
    {
        grFx->SetPoint(i, vFxy[i].xpos * 1e6 - 300, vFxy[i].ypos * 1e6 - 300, vFxy[i].Fx);
        grFy->SetPoint(i, vFxy[i].xpos * 1e6 - 300, vFxy[i].ypos * 1e6 - 300, vFxy[i].Fy);
    }
    grFx->SetTitle("Fx vs Position;X (um);Y (um);Fx");
    grFy->SetTitle("Fy vs Position;X (um);Y (um);Fy");
    grFx->Write("grFx");
    grFy->Write("grFy");

    TH2D *hFx = new TH2D("hFx", "Fx Heatmap;X (um);Y (um);Fx", 10, 0, 300, 10, 0, 300);
    TH2D *hFy = new TH2D("hFy", "Fy Heatmap;X (um);Y (um);Fy", 10, 0, 300, 10, 0, 300);
    for (size_t i = 0; i < vFxy.size(); i++)
    {
        int x_bin = hFx->GetXaxis()->FindBin(vFxy[i].xpos * 1e6 - 300);
        int y_bin = hFx->GetYaxis()->FindBin(vFxy[i].ypos * 1e6 - 300);
        hFx->SetBinContent(x_bin, y_bin, vFxy[i].Fx);
        hFy->SetBinContent(x_bin, y_bin, vFxy[i].Fy);
    }
    hFx->Write("hFx");
    hFy->Write("hFy");

    // 作图：不同x下，Fy与y的关系，不同y下，Fx与x的关系，画散点图
    // 并且用[0]*TMath::Cos([1]*x + [2])拟合Fy与y的关系, Fx与x的关系
    auto fitFunc = [](Double_t *x, Double_t *par)
    {
        return par[0] * TMath::Cos(par[1] * x[0] + par[2]);
    };
    auto fit = new TF1("fit", fitFunc, -150, 150, 3);
    fit->SetParameters(0.6, 0.01, 0);

    std::vector<TGraph *> vtgxFy;
    for (int x = 0; x < x_loop; x++)
    {
        TGraph *tgFy_vs_y = new TGraph();
        tgFy_vs_y->SetTitle(Form("Fy vs y at x_{index}=%d;x (um);Fy", x));
        tgFy_vs_y->SetName(Form("tgFy_vs_y_xindex_%d", x));
        tgFy_vs_y->SetMarkerColor(kBlue);
        tgFy_vs_y->SetMarkerStyle(20);
        tgFy_vs_y->SetLineColor(kBlue);

        for (int y = 0; y < y_loop; y++)
        {
            for (const auto &fxy : vFxy)
            {
                auto pos = GetPosition(x, y);
                if (std::abs(fxy.xpos - pos.first) < 1e-9 && std::abs(fxy.ypos - pos.second) < 1e-9)
                {
                    tgFy_vs_y->SetPoint(y, fxy.ypos * 1e6 - 300, fxy.Fy);
                    break;
                }
            }
        }
        tgFy_vs_y->Fit(fit, "Q", "", 0, y_loop * 30);

        c1->cd();
        if (x == 5)
            tgFy_vs_y->Draw("AP");

        tgFy_vs_y->Write();
        vtgxFy.push_back(tgFy_vs_y);
    }
    c1->SaveAs("Fy_vs_y_all_x.png");
    std::vector<TGraph *> vtgxFx;
    for (int y = 0; y < y_loop; y++)
    {
        TGraph *tgFx_vs_x = new TGraph();
        tgFx_vs_x->SetTitle(Form("Fx vs x at y_{index}=%d;x (um);Fx", y));
        tgFx_vs_x->SetName(Form("tgFx_vs_x_yindex_%d", y));
        tgFx_vs_x->SetMarkerColor(kRed);
        tgFx_vs_x->SetMarkerStyle(20);
        tgFx_vs_x->SetLineColor(kRed);
        for (int x = 0; x < x_loop; x++)
        {
            for (const auto &fxy : vFxy)
            {
                auto pos = GetPosition(x, y);
                if (std::abs(fxy.xpos - pos.first) < 1e-9 && std::abs(fxy.ypos - pos.second) < 1e-9)
                {
                    tgFx_vs_x->SetPoint(x, fxy.xpos * 1e6 - 300, fxy.Fx);
                    break;
                }
            }
        }
        fit->SetLineColor(kBlue);
        tgFx_vs_x->Fit(fit, "Q", "", 0, x_loop * 30);
        tgFx_vs_x->Write();

        if (y == 5)
            tgFx_vs_x->Draw("AP");
        vtgxFx.push_back(tgFx_vs_x);
    }
    c1->SaveAs("Fx_vs_x_all_y.png");

    outFile->Close();
    delete outFile;
}