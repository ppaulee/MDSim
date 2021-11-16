//
// Created by ethan on 11/16/2021.
//

#include "Log.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace MolSim{

    std::shared_ptr<spdlog::logger> Log::s_Logger;

    void Log::Init() {
        spdlog::set_pattern("%^[%T] %n [%l]: %v%$");
        s_Logger = spdlog::stdout_color_mt("MolSim");
        s_Logger->set_level(spdlog::level::trace);
    }

}