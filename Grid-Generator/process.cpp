#include <TSystem.h>
#include "ResistPlateModel.cpp"
void process(int x, int y, double xstep, double ystep)
{
    std::string sFolder = Form("simu/x_%d_y_%d/", x, y);
    gSystem->mkdir(sFolder.c_str(), true);
    gSystem->Exec(Form("cp LGAD-wfm.csv %s 2> /dev/null", sFolder.c_str()));
    // generateRSDModel("simu/");
    // generateRSDModel(sFolder.c_str());
    double pitch = 300e-6;
    double xInj = 0.5 * 2 * pitch + xstep * x;
    double yInj = 0.5 * 2 * pitch + ystep * y;
    generateRSDModel(sFolder, xInj, yInj, false);

    // 运行LTspice仿真
    // 获取当前PWD
    std::string pwd = gSystem->GetWorkingDirectory();
    // 切换到仿真文件夹
    gSystem->ChangeDirectory(sFolder.c_str());
    // gSystem->Exec(Form("powershell.exe \" ltspice.exe -b rsd_4x4_model.cir\" "));
    gSystem->Exec(Form("powershell.exe -command \" ltspice.exe -b -ascii rsd_4x4_model.cir\" "));    // batch mode
    // gSystem->Exec(Form("powershell.exe -command \" ltspice.exe  -ascii rsd_4x4_model.cir\" "));
    // 切换回原PWD
    gSystem->ChangeDirectory(pwd.c_str());
}

void process()
{
    double pitch = 300e-6;
    int nSteps = 2;
    double xstep = pitch / (nSteps - 1);
    double ystep = xstep;
    for (int ix = 0; ix < nSteps; ix++)
        for (int iy = 0; iy < nSteps; iy++)
            process(ix, iy, xstep, ystep);
}