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

            inline static std::shared_ptr<spdlog::logger>& getLogger() {return s_Logger;}
        private:
            static std::shared_ptr<spdlog::logger> s_Logger;
    };
}
// Log macros
#define LOG_INFO(...)   ::MolSim::Log::getLogger()->info(__VA_ARGS__)
#define LOG_TRACE(...)  ::MolSim::Log::getLogger()->trace(__VA_ARGS__)
#define LOG_WARN(...)   ::MolSim::Log::getLogger()->warn(__VA_ARGS__)
#define LOG_ERROR(...)  ::MolSim::Log::getLogger()->error(__VA_ARGS__)

#endif //PSEMOLDYN_GROUPB_LOG_H
