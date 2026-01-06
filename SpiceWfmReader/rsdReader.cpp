// rsd_reader.C - ROOT macro to read LTspice waveform data
// Usage: root -l rsd_reader.C

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

void rsdReader()
{
    const char *filename = "rsd_4x4_model.txt";
    const char *rootfile = "rsd_4x4_model.root";

    std::ifstream infile(filename);
    if (!infile.is_open())
    {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }

    std::string line;

    // Skip header lines until we reach "Values:"
    while (std::getline(infile, line))
    {
        if (line.find("Values:") != std::string::npos)
        {
            break;
        }
    }

    // Variables for TTree
    double time, V_inj_node, I_Iinj;
    double I_VE[16]; // VE1 to VE16

    // Create ROOT file and TTree
    TFile *f = new TFile(rootfile, "RECREATE");
    TTree *tree = new TTree("waveform", "LTspice Waveform Data");

    // Create branches
    tree->Branch("time", &time, "time/D");
    tree->Branch("V_inj_node", &V_inj_node, "V_inj_node/D");
    tree->Branch("I_Iinj", &I_Iinj, "I_Iinj/D");

    // Create branches for VE currents
    for (int i = 0; i < 16; i++)
    {
        std::string branch_name = Form("I_VE%d", i + 1);
        tree->Branch(branch_name.c_str(), &I_VE[i], Form("%s/D", branch_name.c_str()));
    }

    // Vectors for TGraphs
    std::vector<double> v_time;
    std::vector<double> v_V_inj_node;
    std::vector<double> v_I_Iinj;
    std::vector<std::vector<double>> v_I_VE(16); // 16 vectors for VE currents

    int line_count = 0;
    int data_count = 0;

    // Read data
    while (std::getline(infile, line))
    {
        if (line.empty())
            continue;

        // Check if line starts with a number (data line)
        if (isdigit(line[0]) || line[0] == '-')
        {
            std::istringstream iss(line);
            double value;
            int col = 0;

            while (iss >> value)
            {
                if (col == 0)
                {
                    time = value;
                    v_time.push_back(value);
                }
                else if (col == 1)
                {
                    V_inj_node = value;
                    v_V_inj_node.push_back(value);
                }
                else if (col == 2)
                {
                    I_Iinj = value;
                    v_I_Iinj.push_back(value);
                }
                else if (col >= 3 && col <= 18)
                {
                    int ve_index = col - 3;
                    I_VE[ve_index] = value;
                    v_I_VE[ve_index].push_back(value);
                }
                col++;
            }

            if (col >= 19)
            { // Should have 19 columns
                tree->Fill();
                data_count++;
            }
        }
        line_count++;
    }

    std::cout << "Read " << data_count << " data points from " << filename << std::endl;

    // Create TGraphs
    int npoints = v_time.size();

    TGraph *g_V_inj = new TGraph(npoints, &v_time[0], &v_V_inj_node[0]);
    g_V_inj->SetName("g_V_inj_node");
    g_V_inj->SetTitle("Voltage at inj_node;Time (s);Voltage (V)");
    g_V_inj->SetLineColor(kBlue);
    g_V_inj->SetLineWidth(2);

    TGraph *g_I_inj = new TGraph(npoints, &v_time[0], &v_I_Iinj[0]);
    g_I_inj->SetName("g_I_Iinj");
    g_I_inj->SetTitle("Injection Current;Time (s);Current (A)");
    g_I_inj->SetLineColor(kRed);
    g_I_inj->SetLineWidth(2);

    // Create TGraphs for VE currents
    std::vector<TGraph *> g_VE_currents;
    int colors[] = {kBlack, kRed, kBlue, kGreen + 2, kMagenta, kCyan + 1,
                    kOrange + 1, kViolet, kSpring, kTeal, kAzure + 1,
                    kPink + 1, kYellow + 2, kGray, kOrange + 7, kSpring + 5};

    TMultiGraph *mg_VE = new TMultiGraph("mg_VE_currents", "VE Electrode Currents;Time (s);Current (A)");

    for (int i = 0; i < 16; i++)
    {
        TGraph *g = new TGraph(npoints, &v_time[0], &v_I_VE[i][0]);
        std::string name = Form("g_I_VE%d", i + 1);
        std::string title = Form("Current through VE%d;Time (s);Current (A)", i + 1);
        g->SetName(name.c_str());
        g->SetTitle(title.c_str());
        g->SetLineColor(colors[i % 16]);
        g->SetLineWidth(1);

        g_VE_currents.push_back(g);
        mg_VE->Add(g);
    }

    // Write everything to file
    tree->Write();
    g_V_inj->Write();
    g_I_inj->Write();
    mg_VE->Write();

    // Also write individual VE graphs
    for (auto *g : g_VE_currents)
    {
        g->Write();
    }

    // Create some plots for visualization
    TCanvas *c1 = new TCanvas("c1", "Waveforms", 1200, 800);
    c1->Divide(2, 2);

    // Plot 1: Injection voltage
    c1->cd(1);
    gPad->SetGrid(1, 1);
    g_V_inj->Draw("AL");

    // Plot 2: Injection current
    c1->cd(2);
    gPad->SetGrid(1, 1);
    g_I_inj->Draw("AL");

    // Plot 3: All VE currents
    c1->cd(3);
    gPad->SetGrid(1, 1);
    mg_VE->Draw("AL");

    // Plot 4: Zoomed in view of VE currents (first 10 data points)
    c1->cd(4);
    gPad->SetGrid(1, 1);

    // Create a zoomed graph (first 100 points)
    int zoom_points = (npoints > 100) ? 100 : npoints;
    TMultiGraph *mg_zoom = new TMultiGraph("mg_VE_zoom", "VE Currents (Zoomed);Time (s);Current (A)");

    for (int i = 0; i < 8; i++)
    { // Show only first 8 for clarity
        double *x_zoom = new double[zoom_points];
        double *y_zoom = new double[zoom_points];

        for (int j = 0; j < zoom_points; j++)
        {
            x_zoom[j] = v_time[j];
            y_zoom[j] = v_I_VE[i][j];
        }

        TGraph *g_zoom = new TGraph(zoom_points, x_zoom, y_zoom);
        g_zoom->SetLineColor(colors[i]);
        g_zoom->SetLineWidth(2);
        mg_zoom->Add(g_zoom);
    }
    mg_zoom->Draw("AL");

    c1->Write();

    // Clean up
    delete c1;
    f->Close();

    std::cout << "Data saved to " << rootfile << std::endl;
    std::cout << "TTree 'waveform' with " << tree->GetEntries() << " entries" << std::endl;
    std::cout << "TGraphs created for all variables" << std::endl;
}