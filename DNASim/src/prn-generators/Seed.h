// Seed class
// Nicolas Clauvelin


// static class for random seed initialization
// random initialization is based on /dev/random or /dev/urandom


#ifndef DNASim_Seed_h
#define DNASim_Seed_h


#include "DNASim_Includes.h"


namespace DNASim {
	
	
	class Seed {
		
		
	public:
		
		// seed initialization with /dev/random
		static UInteger init_seed_with_dev_random() {
			
            // file opening
			FILE *dev_random;
			dev_random = fopen("/dev/random", "r");
			DS_ASSERT(dev_random != nullptr,
					  "cannot open /dev/random");
			
			// reading
			UInteger seed;
			Size ret = fread(&seed, sizeof(seed), 1, dev_random);
            DS_ASSERT(ret != 0, "/dev/random is not readable");
            fclose(dev_random);
			
			// seed init
			srand(seed);
			
			return seed;

		};
		
		// seed initialization with /dev/random
		static UInteger init_seed_with_dev_urandom() {

			// file opening
			FILE *dev_urandom;
			dev_urandom = fopen("/dev/urandom", "r");
			DS_ASSERT(dev_urandom != nullptr,
					  "cannot open /dev/urandom");
			
			// reading
			UInteger seed;
			Size ret = fread(&seed, sizeof(seed), 1, dev_urandom);
            DS_ASSERT(ret != 0, "/dev/random is not readable");
            assert(ret != 0);
            fclose(dev_urandom);
			
			// seed init
			srand(seed);
			
			return seed;
			
		};
		
		
	private:
		
		Seed() = delete;
        ~Seed() = delete;
		Seed(const Seed& seed) = delete;
        Seed(Seed&& seed) = delete;
        Seed& operator=(const Seed& seed) = delete;
        Seed& operator=(Seed&& seed) = delete;
		
		
	};
	
	
}


#endif	// DNASim_Seed_h
