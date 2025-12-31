#include <TSystem.h>
#include "ResistPlateModel.cpp"
void process()
{
    gSystem->mkdir("simu", true);
    // generateRSDModel("simu/");
    generateRSDModel("simu/");
}