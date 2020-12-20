#include <iostream>

#include "randnum.h"
#include "log.h"
#include "toml.hpp"

#include "ldpclinearsystem.h"

int main()
{
    std::string dirname = "logs";
    std::string logFileName = lab::logger::Log::GetCurrentSystemTime() + "-kmldpc.logger";
    std::ofstream logFile(dirname + "/" + logFileName);
    lab::logger::TeeStream logTee (logFile, std::cout);

    if (logFile.is_open() && logFile.good())
        lab::logger::Log::Get().SetLogStream(logTee);
    else
        lab::logger::Log::Get().SetLogStream(std::cout);

    int flag0, flag1;

    flag0 = 0;
    lab::CLCRandNum::Get().SetSeed(flag0);

    flag1 = 0;
    lab::CWHRandNum::Get().SetSeed(flag1);
    lab::logger::Log::Get().SetLevel(lab::logger::Info);
    LOG(lab::logger::Info, true) << "Start simulation" << std::endl;

    std::ifstream ifs("config.toml", std::ios_base::binary);
    if (ifs.good()) {
        auto arguments = toml::parse(ifs);
        LDPCLinearSystem simulator(arguments);
        simulator.Simulator();
    } else {
        LOG(lab::logger::Info, true) << "Encouter error while opening config.toml" << std::endl;
    }

    LOG(lab::logger::Info, true) << "Simulation done" << std::endl;

    return 0;
}