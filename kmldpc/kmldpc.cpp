﻿#include "randnum.h"
#include "ldpc_linear_system.h"
#include <iostream>
#include "Log.h"
#include <matio.h>

CLCRandNum rndGen0;
CWHRandNum rndGen1;

int main(int argc, char* argv[])
{
    std::ofstream logFile("kmldpc.log");
    kmldpc::TeeStream logTee (logFile, std::cout);

    if (logFile.is_open() && logFile.good())
        kmldpc::Log::get().setLogStream(logTee);
    else
        kmldpc::Log::get().setLogStream(std::cout);

    int flag0, flag1;

    flag0 = 0;
    rndGen0.SetSeed(flag0);

    flag1 = 0;
    rndGen1.SetSeed(flag1);

    std::cout << "[Info] Start simulation" << std::endl;

    LDPC_Linear_System sim_ldpc;
    sim_ldpc.Simulator();

    std::cout << "[Info] Simulation done" << std::endl;

    return 0;
}