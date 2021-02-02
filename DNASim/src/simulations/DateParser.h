// DateParser class
// Nicolas Clauvelin


// simple implementation of a date parser


#ifndef DNASim_DateParser_h
#define DNASim_DateParser_h


#include "DNASim_Includes.h"


namespace DNASim {


    class DateParser {


    public:

        // date function
        static std::string get_date();


    private:

        // private constructor and copy operator
        DateParser() = delete;
        ~DateParser() = delete;
        DateParser(const DateParser& date_parser) = delete;
        DateParser(DateParser&& date_parser) = delete;
        DateParser& operator=(const DateParser& date_parser) = delete;
        DateParser& operator=(DateParser&& date_parser) = delete;


        // time getter function
        // time is formatted as YYYY-MM-DD HH::MM::SS
        static std::string get_current_time();


    };


}


#endif  // DNASim_DateParser_h
