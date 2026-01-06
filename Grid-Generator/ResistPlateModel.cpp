// ResistivePlateModel.cpp (关键部分)
#include "ResistPlateModel.hh"
#include <cmath>
#include <iostream>
#include <TROOT.h>

// 构造函数
ResistivePlateModel::ResistivePlateModel(double pitchX, double pitchY, double rSheet, double cJunction, double gridSize)
    : fPitchX(pitchX), fPitchY(pitchY), fRSheet(rSheet), fCJunction(cJunction), fGridSize(gridSize)
{

    // 计算电路参数
    fRUnit = fRSheet; // 简化近似
    fCUnit = fCJunction * fGridSize * fGridSize;

    // 默认4x4电极，计算平板总尺寸
    fLx = 3 * fPitchX; // 4个电极有3个间隙
    fLy = 3 * fPitchY;
}

ResistivePlateModel::~ResistivePlateModel()
{
    for (auto ele : fElectrodes)
        delete ele;
    fElectrodes.clear();
}

// 构建网格和基本RC网络
void ResistivePlateModel::BuildModel(double electrodeBuffer, int extraBufferCells)
{
    // 全局缓冲偏移
    fElectrodeBuffer = electrodeBuffer;
    fGlobalBuf = extraBufferCells;

    // 计算总平板尺寸
    // 假设4×4电极，间距500μm
    fLx = 3 * fPitchX + 2 * fElectrodeBuffer; // 3个间隙 + 两边缓冲
    fLy = 3 * fPitchY + 2 * fElectrodeBuffer;

    // 计算总网格数
    fNx = int(fLx / fGridSize) + 1 + 2 * fGlobalBuf;
    fNy = int(fLy / fGridSize) + 1 + 2 * fGlobalBuf;

    // 创建节点
    fNodes.clear();
    int nodeId = 0;
    for (int i = 0; i < fNx; i++)
    {
        for (int j = 0; j < fNy; j++)
        {
            GridNode node;
            node.id = nodeId++;

            // 计算实际坐标（考虑所有缓冲）
            node.x = (i - fGlobalBuf) * fGridSize - fElectrodeBuffer;
            node.y = (j - fGlobalBuf) * fGridSize - fElectrodeBuffer;

            node.name = Form("N_%d_%d", i, j);

            // 标记边界节点（最外圈接地）
            node.isBoundary = (i == 0 || i == fNx - 1 ||
                               j == 0 || j == fNy - 1);
            fNodes.push_back(node);
        }
    }

    // 创建电阻和电容
    fResistors.clear();
    fCapacitors.clear();

    // 添加电容和电阻连接
    for (int i = 0; i < fNx; i++)
    {
        for (int j = 0; j < fNy; j++)
        {
            int nodeId = i * fNy + j;

            // 每个节点对地电容
            Capacitor cap;
            cap.name = Form("C_%s", fNodes[nodeId].name.c_str());
            cap.node = nodeId;
            cap.value = fCUnit;
            fCapacitors.push_back(cap);

            // 水平方向电阻（向右连接）
            if (i < fNx - 1)
            {
                int rightNodeId = (i + 1) * fNy + j;
                ConnectNodes(nodeId, rightNodeId, fRUnit, "R_h");
            }

            // 垂直方向电阻（向上连接）
            if (j < fNy - 1)
            {
                int upNodeId = i * fNy + (j + 1);
                ConnectNodes(nodeId, upNodeId, fRUnit, "R_v");
            }
        }
    }
}

// 设置注入点
void ResistivePlateModel::SetInjectionPoint(double relX, double relY, double rSpread, double rDiffuseFactor)
{
    // 计算绝对坐标，注意这里起始点为(-fElectrodeBuffer, -fElectrodeBuffer)
    double absX = relX * fLx - fElectrodeBuffer;
    double absY = relY * fLy - fElectrodeBuffer;
    SetInjectionPointAbsolute(absX, absY, rSpread, rDiffuseFactor);
}

void ResistivePlateModel::SetInjectionPointAbsolute(double absX, double absY, double rSpread, double rDiffuseFactor)
{
    // 找到中心节点
    fInjection.centerNodeId = FindNearestNode(absX, absY);
    fInjection.rSpread = rSpread;
    fInjection.rDiffuse = fRUnit * rDiffuseFactor;

    // 找到最近的4个邻居节点用于扩散电阻
    auto centerNode = fNodes[fInjection.centerNodeId];
    fInjection.neighborNodes.clear();

    // 搜索最近的几个节点（排除自己）
    std::vector<std::pair<double, int>> distances;
    for (const auto &node : fNodes)
    {
        if (node.id == fInjection.centerNodeId)
            continue;
        double dx = node.x - centerNode.x;
        double dy = node.y - centerNode.y;
        double dist = sqrt(dx * dx + dy * dy);
        if (dist < fGridSize * 1.5)
        { // 1.5倍网格内的节点
            distances.push_back({dist, node.id});
        }
    }

    // 取最近的4个
    std::sort(distances.begin(), distances.end());
    for (int i = 0; i < std::min(4, (int)distances.size()); i++)
    {
        fInjection.neighborNodes.push_back(distances[i].second);
    }
}

void ResistivePlateModel::SetDeltaInjection(double totalCharge, double pulseWidth)
{
    fIsPWL = false;
    fInjectionCharge = totalCharge;
    fInjectionPulseWidth = pulseWidth;
}

void ResistivePlateModel::GenerateSPICEInjectionDelta(std::ostream &out, double totalCharge, double pulseWidth)
{
    out << "* 三角形δ脉冲模型 *\n";
    // Generate a pulse current source for injection
    // double totalCharge = 1.6e-19 * 5000;      // 总电荷量 (库仑), 5000 e-h对
    // double pulseWidth = 1e-12;                // 脉冲宽度 (秒)，建议1ps
    double I_peak = totalCharge / pulseWidth; // 峰值电流 (安培)
    double t_rise = pulseWidth / 2.0;
    double t_fall = pulseWidth / 2.0;
    out << Form("Iinj INJ_NODE 0 PWL(0 0 %.3e %.3e %.3e 0) ; 三角形δ脉冲\n",
                t_rise, I_peak,
                t_rise + t_fall);
}

void ResistivePlateModel::SetPWLInjection(const std::string &csvFile)
{
    fIsPWL = true;
    sfPWLFileName = csvFile;
}

void ResistivePlateModel::GenerateSPICEInjectionPWL(std::ostream &out, const std::string &filename)
{
    out << "* PWL注入模型 *\n";
    out << Form("Iinj INJ_NODE 0 PWL FILE=\"%s\" ; PWL 电流源\n", filename.c_str());
}

void ResistivePlateModel::SetSimulationParameters(double totalTime, double timeStep)
{
    fTotalTime = totalTime;
    fTimeStep = timeStep;
}

void ResistivePlateModel::GenerateSPICENetlist(const std::string &filename)
{
    std::ofstream out(filename);

    out << "* RSD 4x4电极阵列阻性平板模型 (ROOT C++生成)\n";
    out << "* 参数: Rsheet=" << fRSheet << " Cjunc=" << fCJunction << " Grid=" << fGridSize << "\n\n";

    // 1. 输出电容
    out << "** 对地电容 **\n";
    for (const auto &cap : fCapacitors)
    {
        out << Form("C%s %s 0 %.3e\n", cap.name.c_str(),
                    fNodes[cap.node].name.c_str(), cap.value);
    }
    out << "\n";

    // 2. 输出电阻
    out << "** 网格电阻 **\n";
    for (const auto &res : fResistors)
    {
        out << Form("%s %s %s %.3e\n", res.name.c_str(),
                    fNodes[res.node1].name.c_str(),
                    fNodes[res.node2].name.c_str(),
                    res.value);
    }
    out << "\n";

    // 3. 边界接地
    out << "** 四周接地边界 **\n";
    for (const auto &node : fNodes)
    {
        if (node.isBoundary)
        {
            out << Form("Vgnd_%s %s 0 0\n", node.name.c_str(), node.name.c_str());
        }
    }
    out << "\n";

    // 4. 电极连接， Rsheet-SiN-电极 模型
    out << "** 电极连接 **\n";
    for (const auto &ele : fElectrodes)
    {
        out << Form("* 电极 %s 在 (%.1f, %.1f)um\n",
                    ele->name.c_str(), ele->x * 1e6, ele->y * 1e6);

        // 4.1. 电极虚拟节点（连接平板）
        if (fVirtualVoltageSource)
        {
            out << Form("* 电极虚拟节点\n");
            out << Form("V%s %s 0 0  ; 用于测量电极总电流\n",
                        ele->name.c_str(), ele->nodePad.c_str());
        }

        // 4.2. 所有SiN耦合电容（连接到同一个Pad节点）
        for (size_t i = 0; i < ele->coupledNodes.size(); i++)
        {
            int nodeId = ele->coupledNodes[i].first;
            double cSiN = ele->coupledNodes[i].second;

            out << Form("CSiN_%s_%d %s %s %.3e\n",
                        ele->name.c_str(), (int)i + 1,
                        fNodes[nodeId].name.c_str(), // 连接到阻性层网格
                        ele->nodePad.c_str(),        // 连接到电极pad
                        cSiN);
        }

        // 4.3. 放大器输入电阻（输出节点到地）
        out << Form("* 放大器输入阻抗\n");
        out << Form("Ramp_%s %s 0 %.3e\n",
                    ele->name.c_str(), ele->nodePad.c_str(), ele->rAmp);

        out << "\n";
    }

    // 5. 电流注入点模型
    out << "** 电流注入点 (面源模型) **\n";

    if (!fIsPWL)
        GenerateSPICEInjectionDelta(out, fInjectionCharge, fInjectionPulseWidth);
    else
        GenerateSPICEInjectionPWL(out, sfPWLFileName);

    // out << Form("Iinj INJ_NODE 0 PULSE(0 1e-6 0 10e-12 10e-12 1e-9 5e-9)\n");
    out << Form("R_spread INJ_NODE %s %.3e\n",
                fNodes[fInjection.centerNodeId].name.c_str(), fInjection.rSpread);

    for (int i = 0; i < fInjection.neighborNodes.size(); i++)
    {
        out << Form("R_diff%d %s %s %.3e\n", i + 1,
                    fNodes[fInjection.centerNodeId].name.c_str(),
                    fNodes[fInjection.neighborNodes[i]].name.c_str(),
                    fInjection.rDiffuse);
    }
    out << "\n";

    // 6. 仿真设置
    out << "** 仿真设置 **\n";
    out << Form(".tran 0 %.3e 0 %.3e\n", fTotalTime, fTimeStep);

    if (!fSaveAll)
    {
        for (const auto &nodeName : fSaveNodeNames)
            if (fVirtualVoltageSource)
                out << Form(".save I(V%s)\n", nodeName.c_str());
            else
                out << Form(".save V(%s)\n", nodeName.c_str());
        out << Form(".save I(Iinj)\n");
        out << Form(".save V(INJ_NODE)\n");
    }

    out << ".backanno\n";
    out << ".end\n";

    // 8. 运行指令
    out << "** 运行指令 **\n";
    // out << ".control\n";
    // out << "run\n";
    // out << "write waveforms.raw all\n";
    // out << ".endc\n";

    out.close();
    std::cout << "网表已生成至: " << filename << std::endl;
    std::cout << "总节点数: " << fNodes.size() << std::endl;
    std::cout << "总电阻数: " << fResistors.size() << std::endl;
    std::cout << "总电容数: " << fCapacitors.size() << std::endl;
}

int ResistivePlateModel::FindNearestNode(double x, double y)
{
    // 将物理坐标转换为网格索引（考虑缓冲偏移）
    int i = int(round((x + fElectrodeBuffer + fGlobalBuf * fGridSize) / fGridSize));
    int j = int(round((y + fElectrodeBuffer + fGlobalBuf * fGridSize) / fGridSize));

    // 边界检查
    i = std::max(0, std::min(i, fNx - 1));
    j = std::max(0, std::min(j, fNy - 1));

    // 返回节点ID
    return i * fNy + j;
}

std::vector<std::pair<int, double>> ResistivePlateModel::FindWeightedNodes(
    double x, double y, int maxNodes)
{

    std::vector<std::pair<int, double>> result;

    // 1. 找到电极所在的网格单元
    // 将坐标调整到考虑缓冲的网格坐标系
    double gridX = (x + fElectrodeBuffer + fGlobalBuf * fGridSize) / fGridSize;
    double gridY = (y + fElectrodeBuffer + fGlobalBuf * fGridSize) / fGridSize;

    int i0 = int(floor(gridX)); // 单元左下角索引
    int j0 = int(floor(gridY));

    // 2. 检查边界
    i0 = std::max(0, std::min(i0, fNx - 2));
    j0 = std::max(0, std::min(j0, fNy - 2));

    // 3. 获取单元四个角节点
    int nodeIndices[4];
    nodeIndices[0] = i0 * fNy + j0;             // 左下
    nodeIndices[1] = (i0 + 1) * fNy + j0;       // 右下
    nodeIndices[2] = i0 * fNy + (j0 + 1);       // 左上
    nodeIndices[3] = (i0 + 1) * fNy + (j0 + 1); // 右上

    // 4. 计算双线性插值权重（与距离的倒数相关）
    double dx = gridX - i0; // 在单元内的相对位置 [0,1]
    double dy = gridY - j0;

    // 四个角的权重（与距离成反比）
    double weights[4];
    weights[0] = (1 - dx) * (1 - dy); // 左下
    weights[1] = dx * (1 - dy);       // 右下
    weights[2] = (1 - dx) * dy;       // 左上
    weights[3] = dx * dy;             // 右上

    // 5. 归一化权重并构建结果（排除权重过小的节点）
    double sumWeight = 0;
    for (int k = 0; k < 4; k++)
    {
        if (weights[k] > 0.01)
        { // 忽略权重小于1%的连接
            sumWeight += weights[k];
        }
    }

    for (int k = 0; k < 4; k++)
    {
        if (weights[k] > 0.01 && result.size() < maxNodes)
        {
            // 权重与距离成反比，所以权重越大，电阻应该越小
            // 电阻值 = 接触电阻基值 / 归一化权重
            double normalizedWeight = weights[k] / sumWeight;
            // 注意：这里返回的是权重因子，实际电阻值将在电极连接时计算
            result.push_back({nodeIndices[k], 1.0 / normalizedWeight});
        }
    }

    // 6. 如果没有找到合适的节点（特殊情况），返回最近节点
    if (result.empty())
    {
        int nearest = FindNearestNode(x, y);
        result.push_back({nearest, 1.0});
    }

    return result;
}

void ResistivePlateModel::ConnectNodes(int id1, int id2, double value, const std::string &prefix)
{
    // 检查是否已存在连接（避免重复）
    for (const auto &res : fResistors)
    {
        if ((res.node1 == id1 && res.node2 == id2) ||
            (res.node1 == id2 && res.node2 == id1))
        {
            return; // 已存在连接
        }
    }

    Resistor res;
    res.name = Form("%s_%s_to_%s", prefix.c_str(),
                    fNodes[id1].name.c_str(), fNodes[id2].name.c_str());
    res.node1 = id1;
    res.node2 = id2;
    res.value = value;

    fResistors.push_back(res);
}

void ResistivePlateModel::AddRectangularElectrode(double centerX, double centerY, double sizeX, double sizeY, double cSiN_per_area, double rAmp, double rContact)
{
    auto ele_pointer = new Electrode(); // 在这里不拥有指针
    Electrode &ele = *ele_pointer;
    ele.x = centerX;
    ele.y = centerY;
    ele.sizeX = sizeX;
    ele.sizeY = sizeY;
    ele.cSiN_per_area = cSiN_per_area;
    ele.rAmp = rAmp;

    // 生成节点名称
    int eleIndex = fElectrodes.size() + 1;
    ele.name = Form("E%d", eleIndex);
    ele.nodePad = Form("PAD_%s", ele.name.c_str());

    // 计算电极覆盖的矩形区域
    double x_min = centerX - sizeX / 2.0;
    double x_max = centerX + sizeX / 2.0;
    double y_min = centerY - sizeY / 2.0;
    double y_max = centerY + sizeY / 2.0;

    double totalCSiN = 0.0;

    // 找到区域内所有网格节点
    std::vector<int> nodes_in_area;
    for (const auto &node : fNodes)
    {
        if (node.x >= x_min && node.x <= x_max &&
            node.y >= y_min && node.y <= y_max)
        {
            double node_x_min = node.x - fGridSize / 2.0;
            double node_x_max = node.x + fGridSize / 2.0;
            double node_y_min = node.y - fGridSize / 2.0;
            double node_y_max = node.y + fGridSize / 2.0;
            // 计算重叠面积
            double overlap_x = std::max(0.0, std::min(x_max, node_x_max) - std::max(x_min, node_x_min));
            double overlap_y = std::max(0.0, std::min(y_max, node_y_max) - std::max(y_min, node_y_min));
            if (overlap_x <= 0 || overlap_y <= 0)
                continue;

            double overlap_area = overlap_x * overlap_y;
            double cSiN = overlap_area * ele.cSiN_per_area;
            totalCSiN += cSiN;
            ele.coupledNodes.push_back({node.id, cSiN});
            nodes_in_area.push_back(node.id);
        }
    }
    ele.coveredGridNodes = nodes_in_area;

    fElectrodes.push_back(ele_pointer);
    if (fVirtualVoltageSource)
        AddWaveformNode(ele.name); // 用虚拟电压源监测电流
    else
        AddWaveformNode(ele.nodePad); // 监测电极pad电压

    // 输出信息
    std::cout << "电极 " << ele.name
              << " (" << sizeX * 1e6 << "×" << sizeY * 1e6 << " μm²)"
              << " 耦合到 " << ele.coupledNodes.size()
              << " 个网格节点" << std::endl;
    std::cout << "  总SiN电容: " << totalCSiN * 1e12 << " pF" << std::endl;
    std::cout << "  电极节点: " << ele.nodePad << std::endl;
}

// main.cpp
#include "ResistPlateModel.hh"
#include <TROOT.h>

void generateRSDModel(std::string sFolder = "./", double xInj = 0.5, double yInj = 0.5, bool kIsRelative = true)
{
    // 1. 定义物理参数 (示例值，根据实际情况修改)
    double pitchX = 300e-6;          // 300 um, Unit: m
    double pitchY = 300e-6;          // 300 um, Unit: m
    double rSheet = 5e3;             // 5 kOhm/sq
    double gridSize = pitchY / 30.0; // 10 um, Unit: m
    // double gridSize = pitchY / 2.0;               // 150 um, Unit: m

    double thickness_Si = 50e-6;                  // 50 um
    double epsilon_Si = 11.7 * 8.854e-12;         // Si介电常数
    double cJunction = epsilon_Si / thickness_Si; // F/m^2

    // 2. 读出电路典型参数
    double rAmp = 50e6;                                 // 50 MΩ (放大器输入阻抗)
    double thickness_SiN = 100e-9;                      // 100 nm
    double epsilon_SiN = 7.5 * 8.854e-12;               // SiN介电常数
    double cSiN_per_area = epsilon_SiN / thickness_SiN; // F/m^2

    // 3. 创建模型
    ResistivePlateModel model(pitchX, pitchY, rSheet, cJunction, gridSize);

    // 4. 构建基本网格（添加1层缓冲网格）
    // model.BuildModel(1);
    double electrodeBuffer = 100e-6;                     // 电极缓冲100um
    int extraBufferCells = 1;                            // 额外1层网格缓冲
    model.BuildModel(electrodeBuffer, extraBufferCells); // 电极缓冲100um，额外1层网格缓冲

    // 5. 添加16个电极（4x4阵列）
    model.SetVirtualVoltageSource(true); // 使用虚拟电压源监测电极电流
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            double centerX = i * pitchX;
            double centerY = j * pitchY;
            // model.AddElectrode(x, y, 0.001); // 1 mOhm接触电阻
            model.AddRectangularElectrode(centerX, centerY, electrodeBuffer, electrodeBuffer, cSiN_per_area, rAmp, 0.001);
        }
    }

    // 6. 设置注入点（在中间4个电极的中心）
    model.SetPWLInjection("LGAD-wfm.csv");
    // model.SetInjectionPoint(0.5, 0.5, 50.0, 0.1);
    if (kIsRelative)
        model.SetInjectionPoint(xInj, yInj, 50.0, 0.1); // 相对坐标
    else
        model.SetInjectionPointAbsolute(xInj, yInj, 50.0, 0.1); // 绝对坐标

    // 7. 仿真设置
    // 7.1 仿真时间与步长
    // model.SetSimulationParameters(20e-9, 0.1e-12); // 总时间20ns，时间步长0.1ps
    model.SetSimulationParameters(150e-9, 10e-12); // 总时间150ns，时间步长10ps
    // 7.2 设置保存波形
    model.SaveAllNodeWaveforms(false); // 不保存所有节点
    // model.SaveAllNodeWaveforms(true); // 保存所有节点

    // 7. 生成SPICE网表
    model.GenerateSPICENetlist(sFolder + "rsd_4x4_model.cir");

    std::cout << "\n模型生成完成！" << std::endl;
    std::cout << "使用命令在终端中运行LTspice仿真:" << std::endl;
    std::cout << "  'wine ~/.wine/drive_c/.../LTspice.exe rsd_4x4_model.cir'" << std::endl;
}