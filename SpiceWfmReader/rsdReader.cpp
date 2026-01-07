// rsd_full_reader.C - ROOT macro to read full LTspice waveform data
// Usage: root -l rsd_full_reader.C

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TString.h>
#include <TList.h>
#include <TObjString.h>
#include <TNamed.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>

// Structure to hold variable information
struct VariableInfo
{
    int index;
    TString name;
    TString type;
    TString comment;
};

struct HeaderInfo
{
    TString title;
    TString date;
    TString plotname;
    TString flags;
    int nVariables;
    int nPoints;
    double offset;
    TString command;
    std::vector<TString> toVector() const
    {
        std::vector<TString> vec;
        vec.push_back(Form("Title: %s", title.Data()));
        vec.push_back(Form("Date: %s", date.Data()));
        vec.push_back(Form("Plotname: %s", plotname.Data()));
        vec.push_back(Form("Flags: %s", flags.Data()));
        vec.push_back(Form("No. Variables: %d", nVariables));
        vec.push_back(Form("No. Points: %d", nPoints));
        vec.push_back(Form("Offset: %.6e", offset));
        vec.push_back(Form("Command: %s", command.Data()));
        return vec;
    }
};

int rsdReader(std::string filename_input = "rsd_4x4_model.raw", std::string rootfile_input = "rsd_4x4_model_full.root")
{
    std::string filename = filename_input;
    std::string rootfile = rootfile_input;

    std::ifstream infile(filename);
    if (!infile.is_open())
    {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return 1;
    }

    std::string line;
    HeaderInfo header_info;
    std::vector<VariableInfo> variables;
    bool in_header = true;
    bool in_variables = false;
    bool in_data = false;

    // ============================================
    // STEP 1: Parse header and variable information
    // ============================================

    std::cout << "\n=== Parsing LTspice Waveform File ===\n"
              << std::endl;

    while (std::getline(infile, line))
    {
        // Skip empty lines
        if (line.empty())
            continue;

        // Check for section transitions
        if (line.find("Title:") == 0)
        {
            header_info.title = line.c_str();
            std::cout << "File: " << line.substr(7) << std::endl;
            continue;
        }
        else if (line.find("Date:") == 0)
        {
            header_info.date = line.c_str();
            std::cout << "Date: " << line.substr(6) << std::endl;
            continue;
        }
        else if (line.find("Plotname:") == 0)
        {
            header_info.plotname = line.c_str();
            std::cout << "Analysis: " << line.substr(10) << std::endl;
            continue;
        }
        else if (line.find("Flags:") == 0)
        {
            header_info.flags = line.c_str();
            std::cout << "Flags: " << line.substr(7) << std::endl;
            continue;
        }
        else if (line.find("No. Variables:") == 0)
        {
            header_info.nVariables = std::stoi(line.substr(15));
            std::cout << "Variables: " << line.substr(15) << std::endl;
            continue;
        }
        else if (line.find("No. Points:") == 0)
        {
            header_info.nPoints = std::stoi(line.substr(12));
            std::cout << "Data points: " << line.substr(12) << std::endl;
            continue;
        }
        else if (line.find("Offset:") == 0)
        {
            header_info.offset = std::stod(line.substr(8));
            continue;
        }
        else if (line.find("Command:") == 0)
        {
            header_info.command = line.c_str();
            std::cout << "Command: " << line.substr(9) << std::endl;
            continue;
        }
        else if (line.find("Variables:") == 0)
        {
            in_variables = true;
            in_header = false;
            std::cout << "\n=== Variable List ===" << std::endl;
            continue;
        }
        else if (line.find("Values:") == 0)
        {
            in_variables = false;
            in_data = true;
            std::cout << "\n=== Reading Data ===" << std::endl;
            break;
        }

        // Parse variable information
        if (in_variables)
        {
            std::istringstream iss(line);
            int index;
            std::string name, comment;

            if (iss >> index >> name >> comment)
            {
                VariableInfo var;
                var.index = index;
                var.name = name.c_str();
                var.comment = comment.c_str();
                // 拆分name，格式：I(VE1) -> type=I, name=VE1
                size_t paren_pos = name.find('(');
                if (paren_pos != TString::kNPOS && name.back() == ')')
                {
                    var.type = name.substr(0, paren_pos).c_str();
                    var.name = name.substr(paren_pos + 1, name.length() - paren_pos - 2).c_str();
                }
                else if (name == "time")
                    var.type = "TIME";
                else
                    var.type = "UNKNOWN";

                variables.push_back(var);

                std::cout << "  [" << index << "] " << var.type << "-" << var.name << " (" << comment << ")" << std::endl;
            }
        }
    }

    std::cout << "Variables in header: " << header_info.nVariables
              << ", Parsed variables: " << variables.size() << std::endl
              << std::endl;
    if (variables.size() != static_cast<size_t>(header_info.nVariables))
    {
        std::cerr << "Error: Mismatch in number of variables!" << std::endl;
        return 1;
    }
    // ============================================
    // STEP 2: Read data
    // ============================================

    // Create ROOT file
    TFile *f = new TFile(rootfile.c_str(), "RECREATE");

    // Store header information in a TNamed object
    auto header_info_vec = header_info.toVector();
    for (size_t i = 0; i < header_info_vec.size(); i++)
    {
        TNamed *header_line = new TNamed(Form("header_%02d", (int)i), header_info_vec[i].Data());
        header_line->Write();
    }

    // Store variable information in a TTree
    TTree *var_tree = new TTree("var", "Variable Information");
    int var_index;
    char var_name[100];
    char var_type[50];
    char var_comment[200];

    var_tree->Branch("index", &var_index, "index/I");
    var_tree->Branch("name", var_name, "name/C");
    var_tree->Branch("type", var_type, "type/C");
    var_tree->Branch("comment", var_comment, "comment/C");

    for (const auto &var : variables)
    {
        var_index = var.index;
        strncpy(var_name, var.name.Data(), 99);
        strncpy(var_type, var.type.Data(), 49);
        strncpy(var_comment, var.comment.Data(), 199);
        var_tree->Fill();
    }
    var_tree->Write();

    // Create data TTree with dynamic branches based on variables
    TTree *data_tree = new TTree("data", "Waveform Data");

    // Arrays to hold data (simplified: we know there are 19 variables)
    const int nVars = variables.size();
    const int nPoints = header_info.nPoints;
    double *data = new double[nVars]; // data[0] = time, data[1] = V(inj_node), etc.

    // Create branches for each variable
    for (int i = 0; i < nVars; i++)
    {
        // Create valid ROOT branch name (remove special characters)
        TString branch_name = variables[i].name;
        branch_name.ReplaceAll("(", "_");
        branch_name.ReplaceAll(")", "");
        branch_name.ReplaceAll(")", "");

        std::string leaf_name = Form("%s/D", branch_name.Data());
        data_tree->Branch(branch_name.Data(), &data[i], leaf_name.c_str());

        std::cout << "Created branch: " << branch_name << " for variable "
                  << variables[i].name << std::endl;
    }

    // Vectors for TGraphs (for key variables)
    std::vector<double> v_time;
    std::vector<double> v_voltage;
    std::vector<double> v_current_inj;
    std::map<int, std::vector<double>> v_currents; // VE currents

    int data_count = 0;

    // Read data lines
    while (std::getline(infile, line))
    {
        if (line.empty())
            continue;

        int read_point;
        int row;
        double value;
        // Check if line starts with data (number or negative sign)
        if (isdigit(line[0]) || line[0] == '-')
        {
            std::istringstream iss(line);
            iss >> read_point >> value;
            if (data_count != read_point)
                std::cerr << "Warning: Expected data point " << data_count << " but got " << read_point << std::endl;

            row = 0;
            data[row++] = value;
        }
        else if (line[0] == ' ' || line[0] == '\t')
        {
            std::istringstream iss(line);
            iss >> value;
            data[row++] = value;
        }

        // Fill tree if row == nVars
        if (row == nVars)
        {
            data_tree->Fill();
            data_count++;
            if (data_count % 1000 == 0)
                std::cout << "\r  Read " << data_count << " data points..." << std::flush;
            else if (data_count == nPoints)
                std::cout << "\r  Read " << data_count << " data points." << std::endl;

            // Store in vectors for TGraphs
            v_time.push_back(data[0]);        // time
            v_voltage.push_back(data[1]);     // V(inj_node)
            v_current_inj.push_back(data[2]); // I(Iinj)

            // Store VE currents (indices 3-18)
            for (int i = 3; i < nVars; i++)
                v_currents[i].push_back(data[i]);
        }
        else if (row > nVars)
            std::cerr << "Warning: Line " << data_count << " has extra data values!" << std::endl;
    }

    std::cout << "\nSuccessfully read " << data_count << " data points" << std::endl;

    // ============================================
    // STEP 3: Create TGraphs for visualization
    // ============================================

    int npoints = v_time.size();

    // Create TGraph for injection voltage
    TGraph *g_voltage = new TGraph(npoints, &v_time[0], &v_voltage[0]);
    g_voltage->SetName("g_voltage_inj_node");
    g_voltage->SetTitle(Form("Voltage at inj_node;Time (s);Voltage (V) - %d points", npoints));
    g_voltage->SetLineColor(kBlue);
    g_voltage->SetLineWidth(2);
    g_voltage->Write();

    // Create TGraph for injection current
    TGraph *g_current_inj = new TGraph(npoints, &v_time[0], &v_current_inj[0]);
    g_current_inj->SetName("g_current_Iinj");
    g_current_inj->SetTitle(Form("Injection Current;Time (s);Current (A) - %d points", npoints));
    g_current_inj->SetLineColor(kRed);
    g_current_inj->SetLineWidth(2);
    g_current_inj->Write();

    // Create TGraphs for VE currents
    TMultiGraph *mg_VE_currents = new TMultiGraph("mg_VE_currents",
                                                  "Currents through VE Electrodes;Time (s);Current (A)");

    int colors[] = {kBlack, kRed, kBlue, kGreen + 2, kMagenta, kCyan + 1,
                    kOrange + 1, kViolet, kSpring, kTeal, kAzure + 1,
                    kPink + 1, kYellow + 2, kGray, kOrange + 7, kSpring + 5};

    for (int i = 3; i < nVars; i++)
    {
        int ve_index = i - 2; // VE1 is index 3 in data array
        if (v_currents.find(i) != v_currents.end() && !v_currents[i].empty())
        {
            TGraph *g = new TGraph(npoints, &v_time[0], v_currents[i].data());
            TString name = Form("g_current_VE%d", ve_index);
            TString title = Form("Current through %s;Time (s);Current (A)",
                                 variables[i].name.Data());
            g->SetName(name);
            g->SetTitle(title);
            g->SetLineColor(colors[(i - 3) % 16]);
            g->SetLineWidth(1);
            g->Write();
            mg_VE_currents->Add(g);
        }
    }
    mg_VE_currents->Write();

    // ============================================
    // STEP 4: Create summary TNamed objects
    // ============================================

    TNamed *summary = new TNamed("summary", Form(
                                                "LTspice Waveform Data\n"
                                                "File: %s\n"
                                                "Variables: %d\n"
                                                "Data points: %d\n"
                                                "Time range: %.2e to %.2e seconds\n"
                                                "Voltage range: %.2e to %.2e V\n"
                                                "Injection current range: %.2e to %.2e A",
                                                filename.c_str(),
                                                nVars,
                                                data_count,
                                                v_time.front(),
                                                v_time.back(),
                                                *std::min_element(v_voltage.begin(), v_voltage.end()),
                                                *std::max_element(v_voltage.begin(), v_voltage.end()),
                                                *std::min_element(v_current_inj.begin(), v_current_inj.end()),
                                                *std::max_element(v_current_inj.begin(), v_current_inj.end())));
    summary->Write();

    // ============================================
    // STEP 5: Create visualization canvas
    // ============================================

    TCanvas *c_main = new TCanvas("c_main", "LTspice Waveform Analysis", 1400, 900);
    c_main->Divide(3, 2);

    // Plot 1: Voltage waveform
    c_main->cd(1);
    gPad->SetGrid(1, 1);
    gPad->SetLogx(1);
    g_voltage->Draw("AL");

    // Plot 2: Current waveform
    c_main->cd(2);
    gPad->SetGrid(1, 1);
    gPad->SetLogx(1);
    g_current_inj->Draw("AL");

    // Plot 3: Combined VE currents
    c_main->cd(3);
    gPad->SetGrid(1, 1);
    gPad->SetLogx(1);
    mg_VE_currents->Draw("AL");

    // Plot 4: Voltage histogram
    c_main->cd(4);
    TH1D *h_voltage = new TH1D("h_voltage", "Voltage Distribution;Voltage (V);Counts",
                               100, v_voltage.front(), v_voltage.back());
    for (auto v : v_voltage)
        h_voltage->Fill(v);
    h_voltage->SetFillColor(kBlue);
    h_voltage->SetLineColor(kBlue);
    h_voltage->Draw();

    // Plot 5: Current histogram
    c_main->cd(5);
    TH1D *h_current = new TH1D("h_current", "Injection Current Distribution;Current (A);Counts",
                               100, v_current_inj.front(), v_current_inj.back());
    for (auto i : v_current_inj)
        h_current->Fill(i);
    h_current->SetFillColor(kRed);
    h_current->SetLineColor(kRed);
    h_current->Draw();

    // Plot 6: Text summary
    c_main->cd(6);
    TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9);
    pt->AddText("LTspice Data Summary");
    pt->AddText("====================");
    pt->AddText(Form("File: %s", filename.c_str()));
    pt->AddText(Form("Total variables: %d", nVars));
    pt->AddText(Form("Data points: %d", data_count));
    pt->AddText(Form("Time range: %.2e to %.2e s", v_time.front(), v_time.back()));
    pt->AddText(Form("Voltage: min=%.2e, max=%.2e V",
                     *std::min_element(v_voltage.begin(), v_voltage.end()),
                     *std::max_element(v_voltage.begin(), v_voltage.end())));
    pt->AddText(Form("Inj. current: min=%.2e, max=%.2e A",
                     *std::min_element(v_current_inj.begin(), v_current_inj.end()),
                     *std::max_element(v_current_inj.begin(), v_current_inj.end())));
    pt->SetFillColor(0);
    pt->SetTextAlign(12);
    pt->Draw();

    c_main->Write();

    // ============================================
    // STEP 6: Cleanup and finalize
    // ============================================

    data_tree->Write();

    // Cleanup
    delete[] data;
    delete c_main;
    f->Close();

    std::cout << "\n=== Processing Complete ===" << std::endl;
    std::cout << "Output saved to: " << rootfile << std::endl;
    std::cout << "Contents:" << std::endl;
    std::cout << "  1. Header information (TNamed objects)" << std::endl;
    std::cout << "  2. 'var' TTree - variable metadata" << std::endl;
    std::cout << "  3. 'data' TTree - " << data_count << " data points" << std::endl;
    std::cout << "  4. TGraph objects for key waveforms" << std::endl;
    std::cout << "  5. Summary canvas with plots" << std::endl;
    std::cout << "\nTo access data in ROOT:" << std::endl;
    std::cout << "  TFile f(\"" << rootfile << "\");" << std::endl;
    std::cout << "  f.ls();  // List contents" << std::endl;
    std::cout << "  TTree* t = (TTree*)f.Get(\"data\");" << std::endl;
    std::cout << "  t->Print();  // Show branches" << std::endl;
    std::cout << "  t->Draw(\"V_inj_node:time\");  // Plot voltage vs time" << std::endl;

    return 0;
}