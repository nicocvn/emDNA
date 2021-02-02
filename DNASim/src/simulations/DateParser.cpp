// DateParser class
// Nicolas Clauvelin


#include "simulations/DateParser.h"


namespace DNASim {


    // date function
    std::string DateParser::get_date() {
        return get_current_time();
    };


    // time getter function
    // time is formatted as YYYY-MM-DD HH::MM::SS
    std::string DateParser::get_current_time() {

        // time data
        time_t rawtime;
        struct tm* timeinfo;
        char buffer[80];

        // current time and conversion to local time
        time(&rawtime);
        timeinfo = localtime(&rawtime);

        // formatting
        strftime(buffer, 80, "%Y-%m-%d %H:%M:%S", timeinfo);

        // string
        return std::string(buffer);

    };


}
