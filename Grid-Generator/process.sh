# 运行root -q -b process.cpp(x,y,xstep,ystep)
# 其中x,y为像素坐标，xstep,ystep为注入位置步长，每次并行16个进程
#!/bin/bash


pitch="300e-6"
nSteps=11
# xstep=$(echo "$pitch / ($nSteps - 1)" | bc -l)
xstep=$(echo "$pitch / ($nSteps - 1)" | root -l | sed 's/.* //' | sed 's/)//')
ystep=$xstep


# 使用seq生成序列
for ix in $(seq 0 1 $((nSteps-1))); do
    for iy in $(seq 0 1 $((nSteps-1))); do
        # 判断是否存在文件夹
        if [ ! -d "simu/x_${ix}_y_${iy}" ]; then
            mkdir "simu/x_${ix}_y_${iy}"
            root -q -b -l process.cpp\($ix,$iy,$xstep,$ystep\) &
        fi
        if (( (ix * nSteps + iy + 1) % 4 == 0 )); then
            wait
        fi
    done
done
wait