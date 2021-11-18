//
// Created by ethan on 11/16/2021.
//

#include "Log.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace MolSim{

    std::shared_ptr<spdlog::logger> Log::consoleLogger;
    std::shared_ptr<spdlog::logger> Log::fileLogger;

    void Log::Init() {
        spdlog::set_pattern("%^[%T] %n [%l]: %v%$");
        consoleLogger = spdlog::stdout_color_mt("MolSimConsole");
        consoleLogger->set_level(spdlog::level::trace);

        fileLogger = spdlog::stdout_color_mt("MolSimFileReader");
        fileLogger->set_level(spdlog::level::trace);
    }

}