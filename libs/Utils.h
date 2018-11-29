#ifndef _Utils_h_
#define _Utils_h_

namespace Utils
{
	template <typename T>
	inline T f_pow(T a, int b)
	{
		T r = 1.0;
	    while(b > 0){
	        r *= a;
	        --b;
	    }
	    return r;
	}		
}

#endif
