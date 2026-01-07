// ve_waveform_analyzer_bipolar.C - Analyze bipolar VE waveforms and Iinj injection current
// Usage: root -l ve_waveform_analyzer_bipolar.C

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TString.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TMath.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TText.h>
#include <TProfile.h>
#include <TMarker.h>
#include <TParameter.h>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>

// Structure to store channel analysis results
struct ChannelAnalysis
{
    int channel;

    // 基线信息
    double baseline;     // 基线值（接近0）
    double baseline_std; // 基线标准差

    // 负信号特征（主信号）
    double signal_neg_max;  // 负信号最大幅度（绝对值）
    double time_at_neg_max; // 负信号最大值对应的时间

    // 正信号特征（次要信号）
    double signal_pos_max;  // 正信号最大幅度
    double time_at_pos_max; // 正信号最大值对应的时间

    // 总体特征
    double peak_to_peak; // 峰峰值（正最大 - 负最大）
    double signal_ratio; // 正负信号幅度比 |正/负|

    // 负信号时间参数
    double neg_rise_time_10_90; // 负信号10%-90%上升时间
    double neg_rise_time_20_80; // 负信号20%-80%上升时间
    double neg_fall_time_90_10; // 负信号90%-10%下降时间
    double neg_fwhm;            // 负信号半高宽

    // 信号积分
    double signal_integral_neg; // 负信号积分（负电荷）
    double signal_integral_pos; // 正信号积分（正电荷）
    double signal_integral_net; // 净积分

    // 特征时间点（针对负信号）
    double time_10_percent_neg;
    double time_20_percent_neg;
    double time_80_percent_neg;
    double time_90_percent_neg;
    double time_50_percent_rising_neg;
    double time_50_percent_falling_neg;

    // 与Iinj的关系
    double correlation_with_Iinj; // 与Iinj波形的相关系数
    double time_delay_to_Iinj;    // 相对于Iinj的时间延迟

    // 其他
    double noise_level;
    bool is_anomalous;
    TString anomaly_type;

    // 波形
    TGraph *waveform;
};

// Iinj注入电流分析结构
struct IinjAnalysis
{
    double baseline;
    double baseline_std;
    double signal_pos_max;  // 正信号最大幅度
    double time_at_pos_max; // 正信号最大值时间
    double signal_integral; // 信号积分（总注入电荷）
    double rise_time_10_90; // 10%-90%上升时间
    double fwhm;            // 半高宽

    // 特征时间点
    double time_10_percent;
    double time_90_percent;
    double time_50_percent_rising;

    TGraph *waveform;
};

class VEWaveformAnalyzerBipolar
{
private:
    TFile *input_file;
    TTree *data_tree;
    TString output_filename;
    TString output_reportname;

    // Analysis parameters
    double baseline_start_time; // 基线开始时间（从后面取）
    double baseline_end_time;   // 基线结束时间

    double signal_search_window_start; // 信号搜索窗口开始
    double signal_search_window_end;   // 信号搜索窗口结束

    // Results storage
    std::map<int, ChannelAnalysis> ve_results;
    IinjAnalysis iinj_result;
    std::vector<int> anomalous_channels;

    // Calculated position estimators
    double fFx;
    double fFy;

public:
    VEWaveformAnalyzerBipolar(const std::string &input_file_name, const std::string &output_file_name = "ve_bipolar_analysis.root")
        : output_filename(output_file_name)
    {
        output_reportname = output_filename;
        output_reportname.ReplaceAll(".root", "_report.txt");

        // 打开输入文件
        input_file = new TFile(input_file_name.c_str(), "READ");
        if (!input_file || input_file->IsZombie())
        {
            std::cerr << "Error: Cannot open input file " << input_file_name << std::endl;
            return;
        }

        // 获取数据树
        data_tree = (TTree *)input_file->Get("data");
        if (!data_tree)
        {
            std::cerr << "Error: Cannot find TTree 'waveform' in file" << std::endl;
            return;
        }

        // 设置分析参数 - 从后面取基线
        baseline_start_time = 4.0e-10; // 400 ps（假设信号在更早时间）
        baseline_end_time = 1.0e-9;    // 1 ns

        // 信号搜索窗口（主要信号应该在这个窗口内）
        signal_search_window_start = 0.0;   // 从t=0开始
        signal_search_window_end = 30.0e-9; // 30 ns
        // signal_search_window_end = 3.0e-10; // 300 ps

        std::cout << "VEWaveformAnalyzerBipolar initialized" << std::endl;
        std::cout << "  Baseline region (from end): " << baseline_start_time << " to " << baseline_end_time << " s" << std::endl;
        std::cout << "  Signal search window: " << signal_search_window_start << " to " << signal_search_window_end << " s" << std::endl;
    }

    ~VEWaveformAnalyzerBipolar()
    {
        if (input_file)
            input_file->Close();
    }

    // 主分析函数
    void AnalyzeAll()
    {
        std::cout << "\n=== Starting Analysis of Bipolar Waveforms ===" << std::endl;

        // 1. 首先分析Iinj注入电流
        std::cout << "\n1. Analyzing Iinj Injection Current..." << std::endl;
        AnalyzeIinj();

        // 2. 分析所有VE通道
        std::cout << "\n2. Analyzing VE Channels..." << std::endl;

        // 获取时间数据
        std::vector<double> time_data = GetTimeData();
        if (time_data.empty())
        {
            std::cerr << "Error: No time data found" << std::endl;
            return;
        }

        int npoints = time_data.size();
        std::cout << "  Number of time points: " << npoints << std::endl;

        // 获取Iinj数据用于相关性分析
        std::vector<double> iinj_data = GetIinjData();

        // 分析每个VE通道
        for (int channel = 1; channel <= 16; channel++)
        {
            std::cout << "\n  Analyzing VE Channel " << channel << "..." << std::endl;

            // 获取该通道的数据
            std::vector<double> current_data = GetChannelData(channel);
            if (current_data.empty())
            {
                std::cerr << "    Warning: No data for VE" << channel << std::endl;
                continue;
            }

            // 创建TGraph
            TGraph *waveform = new TGraph(npoints, &time_data[0], &current_data[0]);
            waveform->SetName(Form("waveform_VE%d", channel));
            waveform->SetTitle(Form("Bipolar Waveform VE%d;Time (s);Current (A)", channel));

            // 执行通道分析（包括与Iinj的关系）
            ChannelAnalysis result = AnalyzeVEChannel(channel, waveform, time_data, current_data, iinj_data);
            ve_results[channel] = result;

            // 检查是否异常
            if (result.is_anomalous)
            {
                anomalous_channels.push_back(channel);
                std::cout << "    * ANOMALY DETECTED: " << result.anomaly_type << std::endl;
            }

            std::cout << "    Negative peak: " << result.signal_neg_max << " A at t=" << result.time_at_neg_max << " s" << std::endl;
            std::cout << "    Positive peak: " << result.signal_pos_max << " A at t=" << result.time_at_pos_max << " s" << std::endl;
            std::cout << "    Neg rise time (10-90%): " << result.neg_rise_time_10_90 << " s" << std::endl;
            std::cout << "    Delay to Iinj: " << result.time_delay_to_Iinj << " s" << std::endl;
        }

        // 计算Fx=(VE6+VE7-VE10-VE11)/(VE6+VE7+VE10+VE11), Fy=(VE6+VE10-VE7-VE11)/(VE6+VE7+VE10+VE11)
        std::cout << "\n3. Calculating Position Estimators Fx and Fy..." << std::endl;
        auto result_ve6 = ve_results.find(6);
        auto result_ve7 = ve_results.find(7);
        auto result_ve10 = ve_results.find(10);
        auto result_ve11 = ve_results.find(11);
        if (result_ve6 != ve_results.end() && result_ve7 != ve_results.end() &&
            result_ve10 != ve_results.end() && result_ve11 != ve_results.end())
        {
            double integral_ve6 = TMath::Abs(result_ve6->second.signal_integral_neg);
            double integral_ve7 = TMath::Abs(result_ve7->second.signal_integral_neg);
            double integral_ve10 = TMath::Abs(result_ve10->second.signal_integral_neg);
            double integral_ve11 = TMath::Abs(result_ve11->second.signal_integral_neg);

            double Fx_numerator = (integral_ve6 + integral_ve7) - (integral_ve10 + integral_ve11);
            double Fx_denominator = (integral_ve6 + integral_ve7) + (integral_ve10 + integral_ve11);
            double Fx = Fx_denominator != 0 ? Fx_numerator / Fx_denominator : 0.0;

            double Fy_numerator = (integral_ve6 + integral_ve10) - (integral_ve7 + integral_ve11);
            double Fy_denominator = (integral_ve6 + integral_ve7) + (integral_ve10 + integral_ve11);
            double Fy = Fy_denominator != 0 ? Fy_numerator / Fy_denominator : 0.0;

            std::cout << "    VE6 integral: " << integral_ve6 << std::endl;
            std::cout << "    VE7 integral: " << integral_ve7 << std::endl;
            std::cout << "    VE10 integral: " << integral_ve10 << std::endl;
            std::cout << "    VE11 integral: " << integral_ve11 << std::endl;
            std::cout << "    Fx: " << Fx << ", Fy: " << Fy << std::endl;
            fFx = Fx;
            fFy = Fy;
        }
        else
            std::cout << "  Warning: Cannot compute Fx and Fy - missing VE6, VE7, VE10, or VE11 data." << std::endl;

        std::cout << "\n=== Analysis Complete ===" << std::endl;
        std::cout << "  Total VE channels analyzed: " << ve_results.size() << std::endl;
        std::cout << "  Anomalous VE channels: " << anomalous_channels.size() << std::endl;
        std::cout << "  Iinj injection charge: " << iinj_result.signal_integral << " C" << std::endl;
    }

    // 获取时间数据
    std::vector<double> GetTimeData()
    {
        std::vector<double> time_data;

        double time_value;
        data_tree->SetBranchAddress("time", &time_value);

        Long64_t nentries = data_tree->GetEntries();
        for (Long64_t i = 0; i < nentries; i++)
        {
            data_tree->GetEntry(i);
            time_data.push_back(time_value);
        }

        return time_data;
    }

    // 获取Iinj数据
    std::vector<double> GetIinjData()
    {
        std::vector<double> iinj_data;

        double iinj_value;
        data_tree->SetBranchAddress("Iinj", &iinj_value);

        Long64_t nentries = data_tree->GetEntries();
        for (Long64_t i = 0; i < nentries; i++)
        {
            data_tree->GetEntry(i);
            iinj_data.push_back(iinj_value);
        }

        return iinj_data;
    }

    // 获取指定VE通道的数据
    std::vector<double> GetChannelData(int channel)
    {
        std::vector<double> channel_data;

        TString branch_name = Form("VE%d", channel);

        double current_value;
        auto branch = data_tree->GetBranch(branch_name.Data());
        if (branch)
        {
            branch->SetAddress(&current_value);
            Long64_t nentries = data_tree->GetEntries();
            for (Long64_t i = 0; i < nentries; i++)
            {
                data_tree->GetEntry(i);
                channel_data.push_back(current_value);
            }
            branch->ResetAddress();
        }

        return channel_data;
    }

    // 分析Iinj注入电流
    void AnalyzeIinj()
    {
        std::vector<double> time_data = GetTimeData();
        std::vector<double> iinj_data = GetIinjData();

        if (time_data.empty() || iinj_data.empty())
        {
            std::cerr << "Error: Cannot get Iinj data" << std::endl;
            return;
        }

        // 创建Iinj波形
        TGraph *iinj_waveform = new TGraph(time_data.size(), &time_data[0], &iinj_data[0]);
        iinj_waveform->SetName("waveform_Iinj");
        iinj_waveform->SetTitle("Iinj Injection Current;Time (s);Current (A)");

        iinj_result.waveform = iinj_waveform;

        int npoints = time_data.size();
        double *x = time_data.data();
        double *y = iinj_data.data();

        // 1. 计算基线（从后面取）
        CalculateBaselineIinj(x, y, npoints);

        // // 2. 基线校正
        // std::vector<double> y_corrected(npoints);
        // for (int i = 0; i < npoints; i++)
        //     y_corrected[i] = y[i] - iinj_result.baseline;

        // 3. 寻找正信号最大值（Iinj是正信号）
        FindPositivePeakIinj(x, y, npoints);

        // 4. 计算特征时间点
        CalculateCharacteristicTimesIinj(x, y, npoints);

        // 5. 计算上升时间
        iinj_result.rise_time_10_90 = iinj_result.time_90_percent - iinj_result.time_10_percent;

        // 6. 计算信号积分（总注入电荷）
        CalculateIntegralIinj(x, y, npoints);

        // 7. 计算FWHM
        CalculateFWHMIinj(x, y, npoints);

        std::cout << "  Iinj positive peak: " << iinj_result.signal_pos_max << " A" << std::endl;
        std::cout << "  Peak time: " << iinj_result.time_at_pos_max << " s" << std::endl;
        std::cout << "  Rise time (10-90%): " << iinj_result.rise_time_10_90 << " s" << std::endl;
        std::cout << "  Total injection charge: " << iinj_result.signal_integral << " C" << std::endl;
        std::cout << "  FWHM: " << iinj_result.fwhm << " s" << std::endl;
    }

    void CalculateBaselineIinj(double *x, double *y, int npoints)
    {
        iinj_result.baseline = 0;
        iinj_result.baseline_std = 0;
    }

    void FindPositivePeakIinj(double *x, double *y, int npoints)
    {
        double max_val = -1e100;
        int max_idx = -1;

        // 在信号搜索窗口内寻找最大值
        for (int i = 0; i < npoints; i++)
            if (x[i] >= signal_search_window_start && x[i] <= signal_search_window_end)
                if (y[i] > max_val)
                {
                    max_val = y[i];
                    max_idx = i;
                }

        if (max_idx >= 0)
        {
            iinj_result.signal_pos_max = max_val;
            iinj_result.time_at_pos_max = x[max_idx];
        }
        else
        {
            // 如果没有在窗口内找到，搜索整个范围
            for (int i = 0; i < npoints; i++)
            {
                if (y[i] > max_val)
                {
                    max_val = y[i];
                    max_idx = i;
                }
            }
            iinj_result.signal_pos_max = max_val;
            iinj_result.time_at_pos_max = x[max_idx];
        }
    }

    void CalculateCharacteristicTimesIinj(double *x, double *y, int npoints)
    {
        double threshold_10 = 0.1 * iinj_result.signal_pos_max;
        double threshold_50 = 0.5 * iinj_result.signal_pos_max;
        double threshold_90 = 0.9 * iinj_result.signal_pos_max;

        bool found_10 = false, found_50 = false, found_90 = false;

        // 寻找上升沿（从开始到最大值）
        for (int i = 0; i < npoints && x[i] <= iinj_result.time_at_pos_max; i++)
        {
            if (!found_10 && y[i] >= threshold_10)
            {
                iinj_result.time_10_percent = x[i];
                found_10 = true;
            }
            if (!found_50 && y[i] >= threshold_50)
            {
                iinj_result.time_50_percent_rising = x[i];
                found_50 = true;
            }
            if (!found_90 && y[i] >= threshold_90)
            {
                iinj_result.time_90_percent = x[i];
                found_90 = true;
            }
        }

        // 如果没有找到，设置默认值
        if (!found_10)
            iinj_result.time_10_percent = 0.0;
        if (!found_50)
            iinj_result.time_50_percent_rising = iinj_result.time_at_pos_max / 2.0;
        if (!found_90)
            iinj_result.time_90_percent = iinj_result.time_at_pos_max * 0.9;
    }

    void CalculateIntegralIinj(double *x, double *y, int npoints)
    {
        double integral = 0.0;
        for (int i = 1; i < npoints; i++)
        {
            double dt = x[i] - x[i - 1];
            double avg_y = 0.5 * (y[i] + y[i - 1]);
            integral += avg_y * dt;
        }
        iinj_result.signal_integral = integral;
    }

    void CalculateFWHMIinj(double *x, double *y, int npoints)
    {
        double half_max = 0.5 * iinj_result.signal_pos_max;
        double left_time = 0.0, right_time = 0.0;

        // 寻找左半高点
        for (int i = 0; i < npoints - 1; i++)
        {
            if (y[i] <= half_max && y[i + 1] >= half_max)
            {
                double x1 = x[i], x2 = x[i + 1];
                double y1 = y[i], y2 = y[i + 1];
                left_time = x1 + (x2 - x1) * (half_max - y1) / (y2 - y1);
                break;
            }
        }

        // 寻找右半高点
        for (int i = npoints - 1; i > 0; i--)
        {
            if (y[i] <= half_max && y[i - 1] >= half_max)
            {
                double x1 = x[i - 1], x2 = x[i];
                double y1 = y[i - 1], y2 = y[i];
                right_time = x1 + (x2 - x1) * (half_max - y1) / (y2 - y1);
                break;
            }
        }

        iinj_result.fwhm = (right_time > left_time) ? (right_time - left_time) : 0.0;
    }

    // 分析VE通道（双极性信号）
    ChannelAnalysis AnalyzeVEChannel(int channel, TGraph *waveform,
                                     const std::vector<double> &time_data,
                                     const std::vector<double> &current_data,
                                     const std::vector<double> &iinj_data)
    {
        ChannelAnalysis result;
        result.channel = channel;
        result.waveform = waveform;

        int npoints = time_data.size();
        const double *x = time_data.data();
        const double *y = current_data.data();

        // 1. 计算基线（从后面取），这里强制设成0
        CalculateBaselineVE(result, x, y, npoints);

        // 3. 寻找负信号最大值（主信号）
        FindNegativePeak(result, x, y, npoints);

        // 4. 寻找正信号最大值（次要信号）
        FindPositivePeakVE(result, x, y, npoints);

        // 5. 计算峰峰值和信号比
        result.peak_to_peak = result.signal_pos_max - result.signal_neg_max; // 注意：neg_max是负值
        result.signal_ratio = (result.signal_pos_max != 0 && result.signal_neg_max != 0) ? TMath::Abs(result.signal_pos_max / result.signal_neg_max) : 0.0;

        // 6. 计算负信号的特征时间点
        CalculateNegativeSignalTimes(result, x, y, npoints);

        // 7. 计算负信号的上升时间和FWHM
        CalculateNegativeSignalParameters(result);

        // 8. 计算信号积分
        CalculateIntegralsVE(result, x, y, npoints);

        // 9. 分析与Iinj的关系
        std::vector<double> y_corrected(npoints);
        for (int i = 0; i < npoints; i++)
            y_corrected[i] = y[i] - result.baseline;
        AnalyzeCorrelationWithIinj(result, time_data, y_corrected, iinj_data);

        // 10. 检查异常
        CheckForAnomaliesVE(result);

        // 11. 计算噪声水平
        result.noise_level = result.baseline_std;

        return result;
    }

    void CalculateBaselineVE(ChannelAnalysis &result, const double *x, const double *y, int npoints)
    {
        result.baseline = 0;
        result.baseline_std = 0;
    }

    void FindNegativePeak(ChannelAnalysis &result, const double *x, const double *y, int npoints)
    {
        double min_val = 1e100; // 寻找最小值（负信号）
        int min_idx = -1;

        // 在信号搜索窗口内寻找最小值
        for (int i = 0; i < npoints; i++)
        {
            if (x[i] >= signal_search_window_start && x[i] <= signal_search_window_end)
            {
                if (y[i] < min_val)
                {
                    min_val = y[i];
                    min_idx = i;
                }
            }
        }

        if (min_idx >= 0)
        {
            result.signal_neg_max = min_val; // 这是负值
            result.time_at_neg_max = x[min_idx];
        }
        else
        {
            // 如果没有在窗口内找到，搜索整个范围
            for (int i = 0; i < npoints; i++)
            {
                if (y[i] < min_val)
                {
                    min_val = y[i];
                    min_idx = i;
                }
            }
            result.signal_neg_max = min_val;
            result.time_at_neg_max = x[min_idx];
        }
    }

    void FindPositivePeakVE(ChannelAnalysis &result, const double *x, const double *y, int npoints)
    {
        double max_val = -1e100;
        int max_idx = -1;

        // 在负峰值之后寻找正峰值
        for (int i = 0; i < npoints; i++)
        {
            if (x[i] > result.time_at_neg_max)
            { // 在负峰值之后
                if (y[i] > max_val)
                {
                    max_val = y[i];
                    max_idx = i;
                }
            }
        }

        if (max_idx >= 0)
        {
            result.signal_pos_max = max_val;
            result.time_at_pos_max = x[max_idx];
        }
        else
        {
            result.signal_pos_max = 0.0;
            result.time_at_pos_max = 0.0;
        }
    }

    void CalculateNegativeSignalTimes(ChannelAnalysis &result, const double *x, const double *y, int npoints)
    {
        // 使用负信号的绝对值来计算百分比
        double signal_amplitude = TMath::Abs(result.signal_neg_max);

        double threshold_10 = -0.1 * signal_amplitude; // 负值
        double threshold_20 = -0.2 * signal_amplitude;
        double threshold_50 = -0.5 * signal_amplitude;
        double threshold_80 = -0.8 * signal_amplitude;
        double threshold_90 = -0.9 * signal_amplitude;

        bool found_10 = false, found_20 = false, found_50 = false;
        bool found_80 = false, found_90 = false;

        // 寻找上升沿（从开始到负峰值）
        for (int i = 0; i < npoints && x[i] <= result.time_at_neg_max; i++)
        {
            if (!found_10 && y[i] <= threshold_10)
            {
                result.time_10_percent_neg = x[i];
                found_10 = true;
            }
            if (!found_20 && y[i] <= threshold_20)
            {
                result.time_20_percent_neg = x[i];
                found_20 = true;
            }
            if (!found_50 && y[i] <= threshold_50)
            {
                result.time_50_percent_rising_neg = x[i];
                found_50 = true;
            }
            if (!found_80 && y[i] <= threshold_80)
            {
                result.time_80_percent_neg = x[i];
                found_80 = true;
            }
            if (!found_90 && y[i] <= threshold_90)
            {
                result.time_90_percent_neg = x[i];
                found_90 = true;
            }
        }

        // 寻找下降沿（从负峰值开始）
        bool found_50_fall = false;
        for (int i = 0; i < npoints; i++)
        {
            if (x[i] >= result.time_at_neg_max)
            {
                if (!found_50_fall && y[i] >= threshold_50)
                {
                    // 线性插值
                    if (i > 0)
                    {
                        double x1 = x[i - 1], x2 = x[i];
                        double y1 = y[i - 1], y2 = y[i];
                        result.time_50_percent_falling_neg = x1 + (x2 - x1) * (threshold_50 - y1) / (y2 - y1);
                    }
                    else
                    {
                        result.time_50_percent_falling_neg = x[i];
                    }
                    found_50_fall = true;
                    break;
                }
            }
        }

        // 设置默认值
        if (!found_10)
            result.time_10_percent_neg = 0.0;
        if (!found_20)
            result.time_20_percent_neg = 0.0;
        if (!found_50)
            result.time_50_percent_rising_neg = result.time_at_neg_max / 2.0;
        if (!found_80)
            result.time_80_percent_neg = result.time_at_neg_max * 0.8;
        if (!found_90)
            result.time_90_percent_neg = result.time_at_neg_max * 0.9;
        if (!found_50_fall)
            result.time_50_percent_falling_neg = result.time_at_neg_max + 1e-11;
    }

    void CalculateNegativeSignalParameters(ChannelAnalysis &result)
    {
        // 上升时间（从10%到90%）
        result.neg_rise_time_10_90 = result.time_90_percent_neg - result.time_10_percent_neg;
        result.neg_rise_time_20_80 = result.time_80_percent_neg - result.time_20_percent_neg;

        // 下降时间（从90%到10%）
        if (result.time_50_percent_falling_neg > result.time_at_neg_max)
            result.neg_fall_time_90_10 = 2.0 * (result.time_50_percent_falling_neg - result.time_at_neg_max);
        else
            result.neg_fall_time_90_10 = 0.0;

        // FWHM
        result.neg_fwhm = result.time_50_percent_falling_neg - result.time_50_percent_rising_neg;
    }

    void CalculateIntegralsVE(ChannelAnalysis &result, const double *x, const double *y, int npoints)
    {
        double integral_neg = 0.0, integral_pos = 0.0;

        for (int i = 1; i < npoints; i++)
        {
            double dt = x[i] - x[i - 1];
            double avg_y = 0.5 * (y[i] + y[i - 1]);

            if (avg_y < 0)
                integral_neg += avg_y * dt; // 负面积
            else
                integral_pos += avg_y * dt; // 正面积
        }

        result.signal_integral_neg = integral_neg;
        result.signal_integral_pos = integral_pos;
        result.signal_integral_net = integral_neg + integral_pos; // 净电荷
    }

    void AnalyzeCorrelationWithIinj(ChannelAnalysis &result,
                                    const std::vector<double> &time_data,
                                    const std::vector<double> &ve_data,
                                    const std::vector<double> &iinj_data)
    {
        if (iinj_data.empty() || ve_data.empty())
            return;

        int n = std::min(time_data.size(), std::min(ve_data.size(), iinj_data.size()));

        // 计算相关系数
        double sum_ve = 0.0, sum_iinj = 0.0;
        double sum_ve2 = 0.0, sum_iinj2 = 0.0, sum_ve_iinj = 0.0;

        for (int i = 0; i < n; i++)
        {
            sum_ve += ve_data[i];
            sum_iinj += iinj_data[i];
            sum_ve2 += ve_data[i] * ve_data[i];
            sum_iinj2 += iinj_data[i] * iinj_data[i];
            sum_ve_iinj += ve_data[i] * iinj_data[i];
        }

        double mean_ve = sum_ve / n;
        double mean_iinj = sum_iinj / n;

        double cov = sum_ve_iinj / n - mean_ve * mean_iinj;
        double std_ve = TMath::Sqrt(sum_ve2 / n - mean_ve * mean_ve);
        double std_iinj = TMath::Sqrt(sum_iinj2 / n - mean_iinj * mean_iinj);

        if (std_ve > 0 && std_iinj > 0)
            result.correlation_with_Iinj = cov / (std_ve * std_iinj);
        else
            result.correlation_with_Iinj = 0.0;

        // 计算时间延迟（VE负峰值相对于Iinj正峰值）
        result.time_delay_to_Iinj = result.time_at_neg_max - iinj_result.time_at_pos_max;
    }

    void CheckForAnomaliesVE(ChannelAnalysis &result)
    {
        result.is_anomalous = false;
        result.anomaly_type = "Normal";

        // 1. 检查负信号幅度是否异常
        double avg_neg_peak = CalculateAverageNegativePeak();
        if (avg_neg_peak < 0)
        { // avg_neg_peak是负值
            double ratio = TMath::Abs(result.signal_neg_max / avg_neg_peak);
            if (ratio > 5.0)
            {
                result.is_anomalous = true;
                result.anomaly_type = "Negative peak too large";
            }
            else if (ratio < 0.2)
            {
                result.is_anomalous = true;
                result.anomaly_type = "Negative peak too small";
            }
        }

        // 2. 检查正负信号比例是否异常
        if (result.signal_ratio > 0.5)
        { // 正信号大于负信号的50%
            result.is_anomalous = true;
            result.anomaly_type = "Positive component too large";
        }

        // 3. 检查上升时间是否异常
        double avg_rise_time = CalculateAverageRiseTime();
        if (avg_rise_time > 0)
        {
            double ratio = result.neg_rise_time_10_90 / avg_rise_time;
            if (ratio > 3.0 || ratio < 0.33)
            {
                result.is_anomalous = true;
                if (ratio > 3.0)
                    result.anomaly_type = "Rise time too long";
                else
                    result.anomaly_type = "Rise time too short";
            }
        }

        // 4. 检查噪声水平
        if (result.noise_level > 1e-9)
        { // 1 nA噪声阈值
            result.is_anomalous = true;
            result.anomaly_type = "High noise";
        }

        // 5. 检查与Iinj的时间延迟是否异常
        double avg_delay = CalculateAverageDelay();
        if (TMath::Abs(result.time_delay_to_Iinj - avg_delay) > 5e-11)
        { // 50 ps
            result.is_anomalous = true;
            result.anomaly_type = "Abnormal timing delay";
        }
    }

    double CalculateAverageNegativePeak()
    {
        double sum = 0.0;
        int count = 0;

        for (const auto &pair : ve_results)
        {
            sum += pair.second.signal_neg_max;
            count++;
        }

        return (count > 0) ? sum / count : 0.0;
    }

    double CalculateAverageRiseTime()
    {
        double sum = 0.0;
        int count = 0;

        for (const auto &pair : ve_results)
        {
            if (pair.second.neg_rise_time_10_90 > 0)
            {
                sum += pair.second.neg_rise_time_10_90;
                count++;
            }
        }

        return (count > 0) ? sum / count : 0.0;
    }

    double CalculateAverageDelay()
    {
        double sum = 0.0;
        int count = 0;

        for (const auto &pair : ve_results)
        {
            sum += pair.second.time_delay_to_Iinj;
            count++;
        }

        return (count > 0) ? sum / count : 0.0;
    }

    // 保存分析结果
    void SaveResults()
    {
        TFile *output_file = new TFile(output_filename.Data(), "RECREATE");

        // 1. 保存Iinj分析结果
        if (iinj_result.waveform)
            iinj_result.waveform->Write();

        // 2. 保存VE波形
        for (auto &pair : ve_results)
        {
            ChannelAnalysis &result = pair.second;
            if (result.waveform)
                result.waveform->Write();
        }

        // 3. 创建汇总TTree
        TTree *summary_tree = new TTree("ve_channel_summary", "VE Channel Bipolar Analysis");

        int ch;
        double baseline, baseline_std;
        double signal_neg_max, time_at_neg_max;
        double signal_pos_max, time_at_pos_max;
        double peak_to_peak, signal_ratio;
        double neg_rise_time_10_90, neg_rise_time_20_80;
        double neg_fall_time_90_10, neg_fwhm;
        double signal_integral_neg, signal_integral_pos, signal_integral_net;
        double correlation_with_Iinj, time_delay_to_Iinj;
        double noise_level;
        bool is_anomalous;
        char anomaly_type[100];

        summary_tree->Branch("channel", &ch, "channel/I");
        summary_tree->Branch("baseline", &baseline, "baseline/D");
        summary_tree->Branch("baseline_std", &baseline_std, "baseline_std/D");
        summary_tree->Branch("signal_neg_max", &signal_neg_max, "signal_neg_max/D");
        summary_tree->Branch("time_at_neg_max", &time_at_neg_max, "time_at_neg_max/D");
        summary_tree->Branch("signal_pos_max", &signal_pos_max, "signal_pos_max/D");
        summary_tree->Branch("time_at_pos_max", &time_at_pos_max, "time_at_pos_max/D");
        summary_tree->Branch("peak_to_peak", &peak_to_peak, "peak_to_peak/D");
        summary_tree->Branch("signal_ratio", &signal_ratio, "signal_ratio/D");
        summary_tree->Branch("neg_rise_time_10_90", &neg_rise_time_10_90, "neg_rise_time_10_90/D");
        summary_tree->Branch("neg_rise_time_20_80", &neg_rise_time_20_80, "neg_rise_time_20_80/D");
        summary_tree->Branch("neg_fall_time_90_10", &neg_fall_time_90_10, "neg_fall_time_90_10/D");
        summary_tree->Branch("neg_fwhm", &neg_fwhm, "neg_fwhm/D");
        summary_tree->Branch("signal_integral_neg", &signal_integral_neg, "signal_integral_neg/D");
        summary_tree->Branch("signal_integral_pos", &signal_integral_pos, "signal_integral_pos/D");
        summary_tree->Branch("signal_integral_net", &signal_integral_net, "signal_integral_net/D");
        summary_tree->Branch("correlation_with_Iinj", &correlation_with_Iinj, "correlation_with_Iinj/D");
        summary_tree->Branch("time_delay_to_Iinj", &time_delay_to_Iinj, "time_delay_to_Iinj/D");
        summary_tree->Branch("noise_level", &noise_level, "noise_level/D");
        summary_tree->Branch("is_anomalous", &is_anomalous, "is_anomalous/O");
        summary_tree->Branch("anomaly_type", anomaly_type, "anomaly_type/C");

        for (const auto &pair : ve_results)
        {
            const ChannelAnalysis &result = pair.second;

            ch = result.channel;
            baseline = result.baseline;
            baseline_std = result.baseline_std;
            signal_neg_max = result.signal_neg_max;
            time_at_neg_max = result.time_at_neg_max;
            signal_pos_max = result.signal_pos_max;
            time_at_pos_max = result.time_at_pos_max;
            peak_to_peak = result.peak_to_peak;
            signal_ratio = result.signal_ratio;
            neg_rise_time_10_90 = result.neg_rise_time_10_90;
            neg_rise_time_20_80 = result.neg_rise_time_20_80;
            neg_fall_time_90_10 = result.neg_fall_time_90_10;
            neg_fwhm = result.neg_fwhm;
            signal_integral_neg = result.signal_integral_neg;
            signal_integral_pos = result.signal_integral_pos;
            signal_integral_net = result.signal_integral_net;
            correlation_with_Iinj = result.correlation_with_Iinj;
            time_delay_to_Iinj = result.time_delay_to_Iinj;
            noise_level = result.noise_level;
            is_anomalous = result.is_anomalous;
            strncpy(anomaly_type, result.anomaly_type.Data(), 99);

            summary_tree->Fill();
        }
        summary_tree->Write();

        // 4. 创建Iinj分析结果的TTree
        TTree *iinj_tree = new TTree("iinj_analysis", "Iinj Injection Current Analysis");
        double iinj_baseline, iinj_baseline_std, iinj_signal_pos_max;
        double iinj_time_at_pos_max, iinj_signal_integral;
        double iinj_rise_time_10_90, iinj_fwhm;

        iinj_tree->Branch("baseline", &iinj_baseline, "baseline/D");
        iinj_tree->Branch("baseline_std", &iinj_baseline_std, "baseline_std/D");
        iinj_tree->Branch("signal_pos_max", &iinj_signal_pos_max, "signal_pos_max/D");
        iinj_tree->Branch("time_at_pos_max", &iinj_time_at_pos_max, "time_at_pos_max/D");
        iinj_tree->Branch("signal_integral", &iinj_signal_integral, "signal_integral/D");
        iinj_tree->Branch("rise_time_10_90", &iinj_rise_time_10_90, "rise_time_10_90/D");
        iinj_tree->Branch("fwhm", &iinj_fwhm, "fwhm/D");

        iinj_baseline = iinj_result.baseline;
        iinj_baseline_std = iinj_result.baseline_std;
        iinj_signal_pos_max = iinj_result.signal_pos_max;
        iinj_time_at_pos_max = iinj_result.time_at_pos_max;
        iinj_signal_integral = iinj_result.signal_integral;
        iinj_rise_time_10_90 = iinj_result.rise_time_10_90;
        iinj_fwhm = iinj_result.fwhm;

        iinj_tree->Fill();
        iinj_tree->Write();

        // 5. 创建直方图和画布
        CreateVisualization(output_file);

        // 6. 保存fFx, fFy
        TParameter<double> _fx("Fx", fFx);
        _fx.Write();
        TParameter<double> _fy("Fy", fFy);
        _fy.Write();

        // 7. 保存文本报告
        SaveTextReport(output_reportname.Data());

        output_file->Close();
        std::cout << "\nResults saved to: " << output_filename << std::endl;
    }

    void CreateVisualization(TFile *output_file)
    {
        output_file->cd();

        // 画布1：所有VE波形叠加（基线校正后）
        TCanvas *c1 = new TCanvas("c_all_ve_waveforms", "All VE Waveforms (Baseline Corrected)", 1400, 900);
        c1->Divide(4, 4);

        for (int ch = 1; ch <= 16; ch++)
        {
            if (ve_results.find(ch) != ve_results.end())
            {
                c1->cd(ch);
                gPad->SetGrid(1, 1);

                TGraph *wf = ve_results[ch].waveform;
                if (wf)
                {
                    wf->SetLineColor(kBlue);
                    wf->Draw("AL");

                    // 标记负峰值
                    TMarker *m_neg = new TMarker(ve_results[ch].time_at_neg_max,
                                                 ve_results[ch].signal_neg_max, 20);
                    m_neg->SetMarkerColor(kRed);
                    m_neg->Draw();

                    // 标记正峰值（如果存在）
                    if (ve_results[ch].signal_pos_max > 0)
                    {
                        TMarker *m_pos = new TMarker(ve_results[ch].time_at_pos_max,
                                                     ve_results[ch].signal_pos_max, 21);
                        m_pos->SetMarkerColor(kGreen);
                        m_pos->Draw();
                    }
                }
            }
        }
        c1->Write();

        // 画布2：Iinj和典型VE波形对比
        TCanvas *c2 = new TCanvas("c_iinj_vs_ve", "Iinj vs VE Waveforms", 1200, 800);
        c2->Divide(2, 2);

        // 子图1：Iinj波形
        c2->cd(1);
        gPad->SetGrid(1, 1);
        if (iinj_result.waveform)
        {
            iinj_result.waveform->SetLineColor(kRed);
            iinj_result.waveform->SetLineWidth(2);
            iinj_result.waveform->Draw("AL");

            TMarker *m_iinj = new TMarker(iinj_result.time_at_pos_max,
                                          iinj_result.signal_pos_max, 29);
            m_iinj->SetMarkerColor(kRed);
            m_iinj->Draw();
        }

        // 子图2：VE负峰值幅度分布
        c2->cd(2);
        std::vector<double> ve_channels;
        std::vector<double> ve_peak_values;
        for (const auto &pair : ve_results)
        {
            ve_channels.push_back(pair.first);
            ve_peak_values.push_back(TMath::Abs(pair.second.signal_neg_max)); // 取绝对值
        }

        Int_t nPoints = ve_channels.size();  // 数据点数
        Double_t *x = ve_channels.data();    // 通道号数组
        Double_t *y = ve_peak_values.data(); // 峰值数组

        // 创建TGraph散点图
        TGraph *scatterPlot = new TGraph(nPoints, x, y);
        scatterPlot->SetTitle("VE Negative Peak Distribution");
        scatterPlot->GetXaxis()->SetTitle("VE Channel Number");
        scatterPlot->GetYaxis()->SetTitle("Peak Amplitude (a.u.)");

        // 设置散点样式
        scatterPlot->SetMarkerStyle(20);    // 实心圆点
        scatterPlot->SetMarkerSize(0.8);    // 点大小
        scatterPlot->SetMarkerColor(kBlue); // 点颜色
        scatterPlot->SetLineColor(kBlue);   // 线颜色（如果有线）

        // 可选：添加Iinj线
        double y_Iinj = TMath::Abs(iinj_result.signal_pos_max);
        TLine *iinjLine = new TLine(x[0], y_Iinj, x[nPoints - 1], y_Iinj);
        iinjLine->SetLineColor(kGreen + 2);
        iinjLine->SetLineStyle(2); // 虚线
        iinjLine->SetLineWidth(2);

        // 设置绘图坐标轴范围，从0到y_Injj的1.2倍
        auto hFrame = new TH1F("hFrame", "VE Negative Peak Distribution", 1, 0, ve_channels.size() + 1);
        hFrame->SetMinimum(-y_Iinj * 0.1);
        hFrame->SetMaximum(y_Iinj * 1.2);
        hFrame->GetXaxis()->SetTitle("VE Channel Number");
        hFrame->GetYaxis()->SetTitle("Peak Amplitude (a.u.)");

        hFrame->SetStats(kFALSE);     // 不显示统计信息
        hFrame->Draw();               // 先绘制坐标轴框架
        hFrame->SetLineColor(kWhite); // 隐藏框架线

        // 绘制散点图
        iinjLine->Draw();
        scatterPlot->Draw("P"); // "A"=绘制坐标轴, "P"=绘制点

        // 可选：添加图例
        TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
        leg->AddEntry(scatterPlot, "VE Peaks", "p");
        leg->AddEntry(iinjLine, "I_{inj} Peak", "l");
        leg->Draw();

        // 子图3：上升时间分布
        c2->cd(3);
        TH1D *h_rise_time = new TH1D("h_rise_time",
                                     "Rise Time Distribution (10-90%);Rise Time (s);Count",
                                     50, 0, 0);

        double rise_time_max = 0.0;
        for (const auto &pair : ve_results)
            if (pair.second.neg_rise_time_10_90 > rise_time_max)
                rise_time_max = pair.second.neg_rise_time_10_90;
        h_rise_time->SetBins(50, 0, rise_time_max * 1.5);

        for (const auto &pair : ve_results)
            if (pair.second.neg_rise_time_10_90 > 0)
                h_rise_time->Fill(pair.second.neg_rise_time_10_90);

        h_rise_time->SetFillColor(kGreen);
        h_rise_time->Draw();

        // 子图4：时间延迟分布
        c2->cd(4);
        TH1D *h_delay = new TH1D("h_delay",
                                 "Time Delay to Iinj;Delay (s);Count",
                                 100, -100e-9, 100e-9);

        for (const auto &pair : ve_results)
            h_delay->Fill(pair.second.time_delay_to_Iinj);

        h_delay->SetFillColor(kMagenta);
        h_delay->Draw();

        c2->Write();

        // 画布3：2D分布图
        TCanvas *c3 = new TCanvas("c_2d_distributions", "2D Distributions", 1200, 800);
        c3->Divide(2, 2);

        // 子图1：负峰值 vs 上升时间
        c3->cd(1);
        TH2D *h2_peak_vs_rise = new TH2D("h2_peak_vs_rise",
                                         "Negative Peak vs Rise Time;Negative Peak (A);Rise Time 10-90% (s)",
                                         50, 0, 0, 50, 0, 0);

        for (const auto &pair : ve_results)
        {
            h2_peak_vs_rise->Fill(TMath::Abs(pair.second.signal_neg_max),
                                  pair.second.neg_rise_time_10_90);
        }
        h2_peak_vs_rise->Draw("COLZ");

        // 子图2：负峰值 vs 时间延迟
        c3->cd(2);
        TH2D *h2_peak_vs_delay = new TH2D("h2_peak_vs_delay",
                                          "Negative Peak vs Time Delay;Negative Peak (A);Delay to Iinj (s)",
                                          50, 0, 0, 50, -1e-10, 1e-10);

        for (const auto &pair : ve_results)
        {
            h2_peak_vs_delay->Fill(TMath::Abs(pair.second.signal_neg_max),
                                   pair.second.time_delay_to_Iinj);
        }
        h2_peak_vs_delay->Draw("COLZ");

        // 子图3：信号积分分布
        c3->cd(3);
        TH1D *h_integral = new TH1D("h_integral",
                                    "Net Signal Integral Distribution;Net Charge (C);Count",
                                    50, 0, 0);

        for (const auto &pair : ve_results)
        {
            h_integral->Fill(TMath::Abs(pair.second.signal_integral_net));
        }

        if (h_integral->GetEntries() > 0)
        {
            double max_val = h_integral->GetMaximum();
            h_integral->SetBins(50, 0, max_val * 1.2);
        }

        h_integral->SetFillColor(kOrange);
        h_integral->Draw();

        // 子图4：正负信号比例
        c3->cd(4);
        TH1D *h_ratio = new TH1D("h_ratio",
                                 "Positive/Negative Signal Ratio;|Pos/Neg| Ratio;Count",
                                 50, 0, 1.0);

        for (const auto &pair : ve_results)
        {
            h_ratio->Fill(pair.second.signal_ratio);
        }

        h_ratio->SetFillColor(kCyan);
        h_ratio->Draw();

        c3->Write();

        // 画布4：异常通道分析
        if (!anomalous_channels.empty())
        {
            TCanvas *c4 = new TCanvas("c_anomalous", "Anomalous Channels", 1200, 800);
            c4->Divide(3, 3);

            int plot_count = 0;
            for (int ch : anomalous_channels)
            {
                if (plot_count >= 9)
                    break; // 最多显示9个

                c4->cd(plot_count + 1);
                gPad->SetGrid(1, 1);

                TGraph *wf = ve_results[ch].waveform;
                if (wf)
                {
                    wf->SetLineColor(kRed);
                    wf->SetLineWidth(2);
                    wf->SetTitle(Form("VE%d - %s", ch, ve_results[ch].anomaly_type.Data()));
                    wf->Draw("AL");
                }

                plot_count++;
            }
            c4->Write();
        }
    }

    void SaveTextReport(std::string report_filename = "ve_bipolar_analysis_report.txt")
    {
        std::ofstream report(report_filename);

        report << "VE Channel Bipolar Waveform Analysis Report\n";
        report << "===========================================\n\n";

        report << "I. Iinj Injection Current Analysis\n";
        report << "   ------------------------------\n";
        report << "   Positive peak: " << iinj_result.signal_pos_max << " A\n";
        report << "   Peak time: " << iinj_result.time_at_pos_max << " s\n";
        report << "   Rise time (10-90%): " << iinj_result.rise_time_10_90 << " s\n";
        report << "   FWHM: " << iinj_result.fwhm << " s\n";
        report << "   Total injection charge: " << iinj_result.signal_integral << " C\n";
        report << "   Baseline: " << iinj_result.baseline << " A (std=" << iinj_result.baseline_std << ")\n\n";

        report << "II. VE Channel Analysis Summary\n";
        report << "   ---------------------------\n";
        report << "   Total channels analyzed: " << ve_results.size() << "\n";
        report << "   Anomalous channels: " << anomalous_channels.size() << "\n";
        report << "   Average negative peak: " << CalculateAverageNegativePeak() << " A\n";
        report << "   Average rise time: " << CalculateAverageRiseTime() << " s\n";
        report << "   Average delay to Iinj: " << CalculateAverageDelay() << " s\n\n";

        report << "III. Channel-by-Channel Results\n";
        report << "   ----------------------------------------------------------------------------------------\n";
        report << "   Ch  NegPeak(A)    Time@Neg(s)   PosPeak(A)    Rise10-90(s)  Delay(s)      NetCharge(C)   Status\n";
        report << "   ----------------------------------------------------------------------------------------\n";

        for (int ch = 1; ch <= 16; ch++)
        {
            if (ve_results.find(ch) != ve_results.end())
            {
                const ChannelAnalysis &r = ve_results[ch];
                report << Form("   %2d  %-12.2e  %-12.2e  %-12.2e  %-12.2e  %-12.2e  %-12.2e  %s\n",
                               ch, r.signal_neg_max, r.time_at_neg_max, r.signal_pos_max,
                               r.neg_rise_time_10_90, r.time_delay_to_Iinj, r.signal_integral_net,
                               r.is_anomalous ? r.anomaly_type.Data() : "Normal");
            }
        }

        report << "\nIV. Anomalous Channels Details\n";
        if (anomalous_channels.empty())
        {
            report << "   No anomalous channels detected.\n";
        }
        else
        {
            for (int ch : anomalous_channels)
            {
                const ChannelAnalysis &r = ve_results[ch];
                report << Form("   VE%d: %s\n", ch, r.anomaly_type.Data());
                report << Form("        Negative peak: %.2e A at %.2e s\n", r.signal_neg_max, r.time_at_neg_max);
                report << Form("        Rise time: %.2e s\n", r.neg_rise_time_10_90);
                report << Form("        Delay to Iinj: %.2e s\n", r.time_delay_to_Iinj);
                report << Form("        Net charge: %.2e C\n\n", r.signal_integral_net);
            }
        }

        report << "\n V. Position Esitimator Calculation\n";
        report << "   --------------------------------\n";
        report << "   Fx = (VE6 + VE7 - VE11 - VE10) / (VE6 + VE7 + VE11 + VE10)\n";
        report << "   Fy = (VE6 + VE10 - VE7 - VE11) / (VE6 + VE7 + VE11 + VE10)\n";
        report << "   Fx = " << fFx << '\t' << "Fy = " << fFy << '\n';

        report.close();
        std::cout << "Text report saved to: ve_bipolar_analysis_report.txt" << std::endl;
    }

    std::string GetOutputReportName() const
    {
        return output_reportname.Data();
    }

    std::string GetOutputRootName() const
    {
        return output_filename.Data();
    }
};

// 主函数
int analyze(std::string input_file = "rsd_4x4_model_full.root",
            std::string output_file = "ve_bipolar_analysis.root")
{
    // 创建分析器
    VEWaveformAnalyzerBipolar analyzer(input_file, output_file);

    // 执行分析
    analyzer.AnalyzeAll();

    // 保存结果
    analyzer.SaveResults();

    std::cout << "\n=== Bipolar Analysis Complete ===" << std::endl;
    std::cout << "Output files created:" << std::endl;
    std::cout << "  1. " << analyzer.GetOutputRootName() << " - ROOT file with all results" << std::endl;
    std::cout << "  2. " << analyzer.GetOutputReportName() << " - Text summary report" << std::endl;
    std::cout << "\nKey parameters analyzed:" << std::endl;
    std::cout << "  - VE negative peak amplitude and timing" << std::endl;
    std::cout << "  - VE positive peak amplitude" << std::endl;
    std::cout << "  - VE rise time (10-90%, 20-80%)" << std::endl;
    std::cout << "  - VE signal integrals (negative, positive, net)" << std::endl;
    std::cout << "  - Correlation and time delay relative to Iinj" << std::endl;
    std::cout << "  - Iinj positive peak and injection charge" << std::endl;
    std::cout << "  - Position estimator Fx,Fy" << std::endl;

    return 0;
}

// 单独分析特定通道
void analyze_single_channel_bipolar(int channel)
{
    const char *input_file = "rsd_4x4_model_full.root";

    TFile *f = new TFile(input_file, "READ");
    if (!f || f->IsZombie())
    {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    TTree *t = (TTree *)f->Get("data");
    if (!t)
    {
        std::cerr << "Error: Cannot find TTree 'waveform'" << std::endl;
        return;
    }

    // 获取时间数据
    std::vector<double> time_data;
    double time_val;
    t->SetBranchAddress("time", &time_val);

    Long64_t nentries = t->GetEntries();
    for (Long64_t i = 0; i < nentries; i++)
    {
        t->GetEntry(i);
        time_data.push_back(time_val);
    }

    // 获取VE通道数据
    TString branch_name = Form("VE%d", channel);
    std::vector<double> ve_data;
    double ve_val;
    t->SetBranchAddress(branch_name.Data(), &ve_val);

    for (Long64_t i = 0; i < nentries; i++)
    {
        t->GetEntry(i);
        ve_data.push_back(ve_val);
    }

    // 获取Iinj数据
    std::vector<double> iinj_data;
    double iinj_val;
    t->SetBranchAddress("Iinj", &iinj_val);

    for (Long64_t i = 0; i < nentries; i++)
    {
        t->GetEntry(i);
        iinj_data.push_back(iinj_val);
    }

    // 创建分析器
    VEWaveformAnalyzerBipolar analyzer(input_file);

    // 创建波形
    TGraph *waveform = new TGraph(time_data.size(), &time_data[0], &ve_data[0]);
    waveform->SetName(Form("VE%d_waveform", channel));

    // 分析该通道（简化版本）
    ChannelAnalysis result = analyzer.AnalyzeVEChannel(channel, waveform, time_data, ve_data, iinj_data);

    // 打印详细结果
    std::cout << "\n=== Detailed Bipolar Analysis for VE" << channel << " ===" << std::endl;
    std::cout << "Negative peak (main signal): " << result.signal_neg_max << " A at t=" << result.time_at_neg_max << " s" << std::endl;
    std::cout << "Positive peak: " << result.signal_pos_max << " A at t=" << result.time_at_pos_max << " s" << std::endl;
    std::cout << "Peak-to-peak: " << result.peak_to_peak << " A" << std::endl;
    std::cout << "Positive/Negative ratio: " << result.signal_ratio << std::endl;
    std::cout << "Rise time (10-90%): " << result.neg_rise_time_10_90 << " s" << std::endl;
    std::cout << "Rise time (20-80%): " << result.neg_rise_time_20_80 << " s" << std::endl;
    std::cout << "FWHM: " << result.neg_fwhm << " s" << std::endl;
    std::cout << "Negative charge: " << result.signal_integral_neg << " C" << std::endl;
    std::cout << "Positive charge: " << result.signal_integral_pos << " C" << std::endl;
    std::cout << "Net charge: " << result.signal_integral_net << " C" << std::endl;
    std::cout << "Correlation with Iinj: " << result.correlation_with_Iinj << std::endl;
    std::cout << "Time delay to Iinj: " << result.time_delay_to_Iinj << " s" << std::endl;
    std::cout << "Baseline: " << result.baseline << " A (std=" << result.baseline_std << ")" << std::endl;
    std::cout << "Status: " << (result.is_anomalous ? result.anomaly_type.Data() : "Normal") << std::endl;

    f->Close();
}