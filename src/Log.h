//
// Created by ethan on 11/16/2021.
//

#ifndef PSEMOLDYN_GROUPB_LOG_H
#define PSEMOLDYN_GROUPB_LOG_H


#include <spdlog/spdlog.h>

namespace MolSim {
    class  Log {
        public:
            static void Init();

            inline static std::shared_ptr<spdlog::logger>& getConsoleLogger() {return consoleLogger;}
            inline static std::shared_ptr<spdlog::logger>& getFileLogger() {return fileLogger;}
        private:
            static std::shared_ptr<spdlog::logger> consoleLogger;
            static std::shared_ptr<spdlog::logger> fileLogger;

    };
}
// ConLog macros
#define LOGC_INFO(...)   ::MolSim::Log::getConsoleLogger()->info(__VA_ARGS__)
#define LOGC_TRACE(...)  ::MolSim::Log::getConsoleLogger()->trace(__VA_ARGS__)
#define LOGC_WARN(...)   ::MolSim::Log::getConsoleLogger()->warn(__VA_ARGS__)
#define LOGC_ERROR(...)  ::MolSim::Log::getConsoleLogger()->error(__VA_ARGS__)

// FileLog macros
#define LOGF_INFO(...)   ::MolSim::Log::getFileLogger()->info(__VA_ARGS__)
#define LOGF_TRACE(...)  ::MolSim::Log::getFileLogger()->trace(__VA_ARGS__)
#define LOGF_WARN(...)   ::MolSim::Log::getFileLogger()->warn(__VA_ARGS__)
#define LOGF_ERROR(...)  ::MolSim::Log::getFileLogger()->error(__VA_ARGS__)

#endif //PSEMOLDYN_GROUPB_LOG_H
