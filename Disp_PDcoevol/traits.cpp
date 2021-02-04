#include"traits.h"

trait::trait()
{
	g_pref = 0;
	g_display = 0;
	g_emig = 0;

	p_display = 0;
	p_pref = 0;
	p_emig = 0;

	viability = 0.0;


}

trait::~trait()
{
	pref1.clear();
	pref2.clear();
	
	emig1.clear();
	emig2.clear();

	display1.clear();
	display2.clear();

}

// deconstructor not empty because we want to clear the vectors