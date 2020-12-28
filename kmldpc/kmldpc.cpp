#include <iostream>
#include <list>
#include <numeric>
#include "randnum.h"
#include "log.h"
#include "thread_pool.h"
#include "toml.hpp"
#include "ldpclinearsystem.h"

int main() {
    std::string dirname = "logs";
    std::string logFileName = lab::logger::Log::get_time() + "-kmldpc.logger";
    std::ofstream logFile(dirname + "/" + logFileName);
    lab::logger::TeeStream logTee(logFile, std::cout);
    if (logFile.is_open() && logFile.good())
        lab::logger::Log::get().set_log_stream(logTee);
    else
        lab::logger::Log::get().set_log_stream(std::cout);
    int flag0, flag1;
    flag0 = 0;
    lab::CLCRandNum::Get().SetSeed(flag0);
    flag1 = 0;
    lab::CWHRandNum::Get().SetSeed(flag1);
    lab::logger::Log::get().set_log_level(lab::logger::Info);
    LOG(lab::logger::Info, true) << "Start simulation" << std::endl;
    std::ifstream ifs("config.toml", std::ios_base::binary);
    if (ifs.is_open()) {
        auto arguments = toml::parse(ifs);
        LDPCLinearSystem simulator(arguments);
        simulator.Simulator();
        ifs.close();
    } else {
        LOG(lab::logger::Info, true) << "Encouter error while opening config.toml" << std::endl;
    }
    LOG(lab::logger::Info, true) << "Simulation done" << std::endl;
    return 0;
}