// UUID class
// Nicolas Clauvelin


#include "prn-generators/dsUUID.h"


#define ds_UUID_SIZE 16


namespace DNASim {


    // static UUID generation (as a string) method
    std::string dsUUID::generate_new_UUID() {

        // /dev/random opening
        FILE *dev_random;
        dev_random = fopen("/dev/random", "r");
        DS_ASSERT(dev_random != nullptr,
                  "cannot open /dev/random");

        // reading random bytes
        unsigned char uu[ds_UUID_SIZE];
        Size ret = fread(uu,
                         sizeof(unsigned char),
                         ds_UUID_SIZE,
                         dev_random);
        assert(ret != 0);

        // /dev/random closing
        fclose(dev_random);

        // convert as UUID
        char string_format[50];
        sprintf(string_format,
                "%02X%02X%02X%02X-%02X%02X-%02X%02X-"
                "%02X%02X-%02X%02X%02X%02X%02X%02X",
                uu[0], uu[1], uu[2], uu[3],
                uu[4], uu[5],
                uu[6], uu[7],
                uu[8], uu[9],
                uu[10], uu[11], uu[12], uu[13], uu[14], uu[15]);

        return std::string(string_format);

    };


}


#undef ds_UUID_SIZE
