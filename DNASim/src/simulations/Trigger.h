// Trigger class
// Nicolas Clauvelin


// simple implementation of an incremental trigger


#ifndef DNASim_Trigger_h
#define DNASim_Trigger_h


#include <DNASim_Includes.h>


namespace DNASim {
	
	
	class Trigger {
		
		
	public:
		
		// constructors
		// default constructor creates a zero limit trigger
		Trigger() = default;
		Trigger(Size trigger_limit);
		Trigger(const Trigger& trigger) = default;
        Trigger(Trigger&& trigger) = default;
		~Trigger() = default;

        // copy and move operators
        Trigger& operator=(const Trigger& trigger) = default;
        Trigger& operator=(Trigger&& trigger) = default;
		
		// trigger limit accessor/modifer
		void set_trigger_limit(Size trigger_lim);
		Size trigger_limit() const;
		
		// management method
		bool click(); // true if the limit is realized
		void reset();
		
		
	private:
		
		Size m_trigger_limit;
		Size m_count;
		
		
	};
	
	
}


#endif	// DNASim_Trigger_h

