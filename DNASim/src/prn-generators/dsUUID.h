// UUID class
// Nicolas Clauvelin


// static class for random UUID generation


#ifndef DNASim_UUID_h
#define DNASim_UUID_h


#include <DNASim_Includes.h>


namespace DNASim {


	class dsUUID {


    public:

        // static UUID generation (as a string) method
        static std::string generate_new_UUID();


    private:

		dsUUID() = delete;
        ~dsUUID() = delete;
		dsUUID(const dsUUID& uuid) = delete;
        dsUUID(dsUUID&& uuid) = delete;
        dsUUID& operator=(const dsUUID& uuid) = delete;
        dsUUID& operator=(dsUUID&& uuid) = delete;


    };


}


#endif  // DNASim_UUID_h
