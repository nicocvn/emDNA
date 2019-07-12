// Trigger class
// Nicolas Clauvelin


#include <Trigger.h>


namespace DNASim {
	
	
	// class constructor with trigger limit initialization
	Trigger::Trigger(Size trigger_lim):
        m_trigger_limit(trigger_lim),
        m_count(0) {};
	
	
	// trigger limit accessor/modifer
	void Trigger::set_trigger_limit(Size trigger_lim) {
		m_trigger_limit = trigger_lim;
	};
	Size Trigger::trigger_limit() const {
		return m_trigger_limit;
	};
	
	
	// tick method
	bool Trigger::click() {
		
		// increment the counter
		++m_count;
		
		// return true if the limit is realized
		return m_count >= m_trigger_limit ? true : false;
		
	};
	
	
	// reset method
	void Trigger::reset() {
		m_count = 0;
	};
	
	
}
