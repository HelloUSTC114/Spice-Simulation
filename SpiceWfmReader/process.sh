#!/bin/bash

# 脚本名称: batch_process_rsd.sh
# 功能: 批量处理x=0到X_MAX, y=0到Y_MAX的所有组合
# 用法: ./batch_process_rsd.sh <workdir> <savedir> <x_max> <y_max>

# 设置参数
workdir=${1:-"/mnt/e/Project-LGAD2/SPICE/Spice-Simulation/Grid-Generator/simu/"}  # 第一个参数：工作目录
savedir=${2:-"/mnt/e/Project-LGAD2/SPICE/Spice-Simulation/SpiceWfmReader/ana/"}  # 第二个参数：保存目录
x_max=${3:-10}                    # 第三个参数：x最大值（默认1）
y_max=${4:-10}                    # 第四个参数：y最大值（默认10）

# 检查必要的参数是否提供
if [ -z "$workdir" ] || [ -z "$savedir" ]; then
    echo "错误: 请提供工作目录和保存目录"
    echo "用法: $0 <workdir> <savedir> [x_max] [y_max]"
    echo "示例: $0 /data/work /data/output 5 8"
    exit 1
fi

# 检查x_max和y_max是否为数字
re='^[0-9]+$'
if ! [[ $x_max =~ $re ]] ; then
    echo "错误: x_max必须是整数"
    exit 1
fi

if ! [[ $y_max =~ $re ]] ; then
    echo "错误: y_max必须是整数"
    exit 1
fi

echo "========== 批量处理开始 =========="
echo "工作目录: $workdir"
echo "保存目录: $savedir"
echo "x范围: 0 到 $x_max"
echo "y范围: 0 到 $y_max"
echo "总组合数: $(($x_max + 1)) × $(($y_max + 1)) = $(( ($x_max + 1) * ($y_max + 1) ))"
echo "开始时间: $(date)"
echo "========================================"

# 成功和失败计数器
success_count=0
fail_count=0
skip_count=0

# 循环x从0到x_max
for ((x=0; x<=x_max; x++)); do
    echo "---------- 处理x=$x ----------"
    
    # 循环y从0到y_max
    for ((y=0; y<=y_max; y++)); do
        echo "处理 x=$x, y=$y"
        
        # 构建输入文件路径
        input_file="$workdir/x_${x}_y_${y}/rsd_4x4_model.raw"
        echo "  输入文件: $input_file"
        
        # 检查输入文件是否存在
        if [ ! -f "$input_file" ]; then
            echo "  警告: 输入文件不存在，跳过"
            ((skip_count++))
            continue
        fi
        
        # 构建输出目录路径
        output_dir="$savedir/x_${x}_y_${y}/"
        if [ -d "$output_dir" ]; then
            echo "  警告: 输出文件夹已存在，跳过: $output_dir"
            ((skip_count++))
            continue  # 跳过当前循环
        fi
        
        # 创建输出目录（如果不存在）
        mkdir -p "$output_dir"
        
        if [ ! -d "$output_dir" ]; then
            echo "  错误: 无法创建输出目录，跳过"
            ((fail_count++))
            continue
        fi
        
        # 第一次运行ROOT程序
        echo "  第一次运行ROOT分析..."
        output_root1="$output_dir/rsd_4x4_model_full.root"
        root -q -b "rsdReader.cpp(\"$input_file\", \"$output_root1\")"
        
        # 检查第一次运行是否成功
        if [ $? -ne 0 ]; then
            echo "  警告: 第一次ROOT运行可能失败，继续尝试第二次..."
        fi
        
        # 第二次运行ROOT程序
        echo "  第二次运行ROOT分析..."
        output_root2="$output_dir/ve_bipolar_analysis.root"
        root -q -b "analyze.cpp(\"$output_root1\", \"$output_root2\")"
        
        # 检查第二次运行是否成功
        if [ $? -eq 0 ]; then
            echo "  ✓ 处理成功"
            ((success_count++))
        else
            echo "  ✗ 处理失败"
            ((fail_count++))
        fi
    done
    
    echo ""
done

echo "========== 批量处理完成 =========="
echo "结束时间: $(date)"
echo "----------------------------------------"
echo "统计信息:"
echo "  成功: $success_count"
echo "  失败: $fail_count"
echo "  跳过: $skip_count"
echo "  总计: $((success_count + fail_count + skip_count))"
echo "========================================"

# 如果所有文件都失败，则退出状态为1
if [ $success_count -eq 0 ] && [ $fail_count -gt 0 ]; then
    exit 1
fi