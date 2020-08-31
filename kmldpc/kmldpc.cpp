#include "randnum.h"
#include "ldpc_linear_system.h"
#include <iostream>

CLCRandNum rndGen0;
CWHRandNum rndGen1;

int main(int argc, char* argv[])
{
    int flag0, flag1;

    flag0 = 0;
    rndGen0.SetSeed(flag0);

    flag1 = 0;
    rndGen1.SetSeed(flag1);

    std::cout << "[Info] Start simulation" << std::endl;

    LDPC_Linear_System sim_ldpc;
    sim_ldpc.Simulator();

    std::cout << "[Info] Simulation done" << std::endl;

    system("pause");

    return 0;
}