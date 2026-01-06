
// ResistivePlateModel.h
#ifndef RESISTIVE_PLATE_MODEL_H
#define RESISTIVE_PLATE_MODEL_H

#include <vector>
#include <string>
#include <fstream>

/*
 * ============================================================================
 * RSD（Resistive Silicon Detector）阻性硅探测器SPICE模型生成器
 * 作者: [John]
 * 版本: 1.0
 * 日期: [2025-Dec-30]
 * ============================================================================
 *
 * 模型概述:
 * 本代码生成用于LTspice仿真的RSD探测器SPICE网表文件。
 * 模型模拟4×4电极阵列的阻性平板，包含完整的读出电路和δ电流注入。
 *
 * 物理模型:
 * 1. 阻性平板: 离散化为RC传输线网络
 * 2. 电极: 有限尺寸(100×100 μm²)，通过加权电阻连接多个网格节点
 * 3. 读出链: GND - RSD Plate - 电极 - C_couple - R_amp - GND
 * 4. 注入: δ电流源模拟粒子入射产生的瞬时电荷
 *
 * ============================================================================
 * 坐标系统说明（关键！）:
 * ============================================================================
 *
 * 全局坐标系定义:
 *
 * 实际实现:
 *
 *  整个平板物理范围:
 *      X: [-fElectrodeBuffer, 3*fPitchX + fElectrodeBuffer]
 *      Y: [-fElectrodeBuffer, 3*fPitchY + fElectrodeBuffer]
 *      示例: 当pitch=300μm, buffer=100μm时
 *          X范围: [-100μm, 3×300+100 = 1000μm]
 *          Y范围: [-100μm, 1000μm]
 *
 *  电极阵列布局（4×4）:
 *      电极中心位于网格交点，采用"左上角原点"约定:
 *          - 第(0,0)个电极中心在 (0, 0)μm
 *          - 第(i,j)个电极中心在 (i×pitchX, j×pitchY)μm
 *
 *      电极实际覆盖区域（100×100 μm²）:
 *          对于电极(i,j)，覆盖区域为:
 *              X: [i×pitchX - 50μm, i×pitchX + 50μm]
 *              Y: [j×pitchY - 50μm, j×pitchY + 50μm]
 *
 *  重要坐标点示例（单位μm）:
 *      电极(0,0): 中心(0,0), 范围[-50,50]×[-50,50]
 *      电极(1,0): 中心(300,0), 范围[250,350]×[-50,50]
 *      电极(0,1): 中心(0,300), 范围[-50,50]×[250,350]
 *      电极(3,3): 中心(900,900), 范围[850,950]×[850,950]
 *
 *  注入点坐标:
 *      相对坐标(0-1)对应电极区域范围[0, 3×pitch]
 *      例如(0.5,0.5)对应(450,450)μm（中间4电极中心）
 *
 * ============================================================================
 * 网格离散化:
 * ============================================================================
 *
 *  网格尺寸: fGridSize (例如15μm)
 *  总网格数: (fNx, fNy)
 *
 *  节点坐标计算:
 *      node.x = (i - fGlobalBuf) * fGridSize - fElectrodeBuffer
 *      node.y = (j - fGlobalBuf) * fGridSize - fElectrodeBuffer
 *
 *  其中:
 *      fGlobalBuf: 边界缓冲层数（数值稳定）
 *      fElectrodeBuffer: 电极外物理缓冲距离（100μm）
 *
 *  电极到网格的映射:
 *      每个电极覆盖多个网格节点
 *      使用加权连接（电阻值与距离成反比）
 *      实现从连续电极到离散网格的物理精确映射
 *
 * ============================================================================
 * 电路元件说明:
 * ============================================================================
 *
 *  1. 平板RC网络:
 *      R_unit ≈ R_sheet (单位电阻)
 *      C_unit = C_junction × fGridSize² (单位电容)
 *
 *  2. 电极连接:
 *      Vmeas_E*: 0V电压源，用于测量电极电流
 *      R{E*}_*: 加权连接电阻（μΩ级）
 *      Ccouple_E*: 交流耦合电容（pF级）
 *      Ramp_E*: 放大器输入阻抗（MΩ级）
 *
 *  3. 电流注入:
 *      Iinj: δ电流脉冲源
 *      R_spread: 垂直扩展电阻（模拟电荷"涌入"阻抗）
 *      R_diff*: 局域横向扩散电阻
 *
 * ============================================================================
 * 物理参数典型值:
 * ============================================================================
 *
 *  几何参数:
 *      pitchX = pitchY = 300e-6 m (300 μm)
 *      gridSize = pitch/20 = 15e-6 m (15 μm)
 *      electrodeBuffer = 100e-6 m (100 μm)
 *
 *  材料参数:
 *      R_sheet = 5e3 Ω/□ (5 kΩ/□)
 *      C_junction = 2e-6 F/m² (0.2 pF/μm²)
 *
 *  信号参数:
 *      δ注入电荷 = 1.6e-19 × 5000 ≈ 8e-16 C (5000个电子)
 *      脉冲宽度 = 1e-12 s (1 ps)
 *      峰值电流 ≈ 0.8 mA
 *
 *  读出参数:
 *      C_couple = 3e-12 F (3 pF)
 *      R_amp = 50e6 Ω (50 MΩ)
 *
 * ============================================================================
 * 使用流程:
 * ============================================================================
 *
 *  1. 初始化模型: ResistivePlateModel(pitchX, pitchY, R_sheet, C_junction, gridSize)
 *  2. 构建网格: BuildModel(electrodeBuffer, boundaryLayers)
 *  3. 添加电极: AddRectangularElectrode(centerX, centerY, sizeX, sizeY, ...)
 *  4. 设置注入: SetInjectionPoint(x, y, ...)
 *  5. 生成网表: GenerateSPICENetlist("filename.cir")
 *  6. LTspice仿真: 导入生成的.cir文件
 *
 * ============================================================================
 * 输出文件结构:
 * ============================================================================
 *
 *  生成的文件包含:
 *      - 文件头和信息
 *      - 平板RC网格定义
 *      - 边界接地条件
 *      - 16个电极的完整读出电路
 *      - δ电流注入模型
 *      - 仿真设置和测量指令
 *
 *  关键测量信号:
 *      I(Vmeas_E*): 各电极总电流
 *      V(OUT_E*): 各电极输出信号电压
 *      .meas指令自动计算电荷量等参数
 *
 * ============================================================================
 * 验证方法:
 * ============================================================================
 *
 *  1. 对称性检查: 中心注入时，中间4电极信号应近似相等
 *  2. 电荷守恒: 所有电极电流之和 ≈ 注入电流
 *  3. 网格收敛: 减小gridSize，结果变化应小于1%
 *  4. 物理合理性: 信号幅度、时间常数应符合预期
 *
 * ============================================================================
 * 注意事项:
 * ============================================================================
 *
 *  1. 坐标一致性: 所有物理坐标使用相同参考系（本文件中明确定义）
 *  2. 单位统一: 所有计算使用SI单位（米，秒，法拉，欧姆）
 *  3. 加权连接: 电极与网格节点的连接电阻反比于几何权重
 *  4. 边界处理: 四周接地，通过缓冲层减小边界效应
 *  5. 时间步长: 仿真需使用足够小时同步长（如0.1ps）以分辨ps级脉冲
 *
 * ============================================================================
 * 修订历史:
 * ============================================================================
 *
 *  v1.0 [日期] 初始版本
 *      - 实现基本RC网格
 *      - 添加有限尺寸电极支持
 *      - 包含完整读出电路
 *      - 支持δ电流注入
 *
 *  // 在此添加后续修订说明
 *
 * ============================================================================
 */

struct GridNode
{
    int id;
    double x, y;
    bool isBoundary;
    std::string name;
};

struct Resistor
{
    std::string name;
    int node1, node2;
    double value;
};

struct Capacitor
{
    std::string name;
    int node;
    double value;
};

struct Electrode
{
    std::string name;
    double x, y;
    double sizeX, sizeY; // 电极尺寸
    double rContact;
    // std::vector<std::pair<int, double>> connectedNodes; //  连接到平板的节点， 节点ID + 权重因子
    std::vector<std::pair<int, double>> coupledNodes; //  NodeId, cSiN

    // double cCouple; // 可选：耦合电容
    double cSiN_per_area; // SiN 单位面积电容(F/m^2)
    double rAmp;          // 可选：放大器输入电阻

    // 连接到多个网格节点（整个电极区域）
    std::vector<int> coveredGridNodes; // 电极覆盖的所有网格节点

    std::string nodePad; // 电极节点（C_couple和R_amp之间）
    // std::string nodeOutput;    // 输出信号节点（C_couple和R_amp之间）
};

// 注入点信息
struct InjectionPoint
{
    int centerNodeId;
    double rSpread, rDiffuse;
    std::vector<int> neighborNodes;
};

class ResistivePlateModel
{
public:
    // 构造函数：初始化物理参数
    ResistivePlateModel(double pitchX, double pitchY, double rSheet, double cJunction, double gridSize);
    ~ResistivePlateModel();

    // 构建模型核心方法
    void BuildModel(double electrodeBuffer = 100e-6, // 电极外缓冲尺寸 (100μm)
                    int extraBufferCells = 2);       // 额外的边界缓冲层数
    void GenerateSPICENetlist(const std::string &filename);

    // 注入点设置
    /// 相对坐标，0-1之间
    void SetInjectionPoint(double relX, double relY, double rSpread = 50.0, double rDiffuseFactor = 0.1);
    /// 绝对坐标，单位米
    void SetInjectionPointAbsolute(double absX, double absY, double rSpread = 50.0, double rDiffuseFactor = 0.1);

    // 电流源设置
    /// @param totalCharge 总电荷量 (库仑), 5000 e-h对
    /// @param pulseWidth 脉冲宽度 (秒)，建议1ps
    void SetDeltaInjection(double totalCharge = 1.6e-19 * 5000, double pulseWidth = 1e-12);
    void GenerateSPICEInjectionDelta(std::ostream &out, double totalCharge = 1.6e-19 * 5000, double pulseWidth = 1e-12);

    void SetPWLInjection(const std::string &csvFile);
    void GenerateSPICEInjectionPWL(std::ostream &out, const std::string &filename);

    // 仿真设置
    void SetSimulationParameters(double totalTime = 20e-9, double timeStep = 0.1e-12);
    void SaveAllNodeWaveforms(bool saveAll = true) { fSaveAll = saveAll; }
    void AddWaveformNode(const std::string &nodeName) { fSaveNodeNames.push_back(nodeName); }

    // 电极定义（添加电极并自动计算最近节点）
    // void AddElectrode(double absX, double absY, double cSiN_per_area, double rAmp, double rContact = 0.001);
    void AddRectangularElectrode(double centerX, double centerY, // 电极中心
                                 double sizeX, double sizeY,     // 电极尺寸 (100e-6 m)
                                 double cSiN_per_area, double rAmp,
                                 double rContact = 0.001);
    void SetVirtualVoltageSource(bool useVirtual) { fVirtualVoltageSource = useVirtual; }

private:
    // 物理参数
    double fPitchX, fPitchY; // 电极间距
    double fRSheet;          // 面电阻
    double fCJunction;       // 单位面积结电容 (F/um^2)
    double fGridSize;        // 网格尺寸

    // 计算出的电路参数
    double fRUnit; // 单位电阻
    double fCUnit; // 单位电容

    // 模型几何
    int fNx, fNy;    // X,Y方向网格总数
    double fLx, fLy; // 平板总尺寸

    double fElectrodeBuffer;            // 电极外缓冲尺寸
    int fGlobalBuf;                     // 全局缓冲网格数
    bool fVirtualVoltageSource = false; // 是否使用虚拟电压源测量电流

    // 电流源设置
    bool fIsPWL = false;         // 是否使用PWL电流源，默认不使用，使用瞬时三角脉冲
    double fInjectionCharge;     // 注入电荷量
    double fInjectionPulseWidth; // 注入脉冲宽度

    std::string sfPWLFileName = "injection_pwl.txt"; // PWL文件名

    // 仿真设置
    double fTotalTime; // 总仿真时间
    double fTimeStep;  // 时间步长

    bool fSaveAll = true;                    // 是否保存所有节点波形
    std::vector<std::string> fSaveNodeNames; // 需要单独保存波形的节点名称

    // 节点与元件存储

    std::vector<GridNode> fNodes;

    std::vector<Resistor> fResistors;

    std::vector<Capacitor> fCapacitors;

    std::vector<Electrode *> fElectrodes;

    // 注入点信息
    InjectionPoint fInjection;

    // 内部方法
    int FindNearestNode(double x, double y);
    std::vector<std::pair<int, double>> FindWeightedNodes(double x, double y, int maxNodes = 4);
    void ConnectNodes(int id1, int id2, double value, const std::string &prefix);
};
#endif