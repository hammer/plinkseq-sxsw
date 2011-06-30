#include "extra.h"

#include <iostream>
#include "pseq.h"

#define MATHLIB_STANDALONE
#include "Rmath.h"

extern GStore * GP;

void extra_code()
{

  for (int r = 0 ; r < 100 ; r++ )
    {
      for (int n = 0 ; n < 100 ; n++)
	for (int x = 0 ; x <= n ; x++ )
	  for (int pi = 0 ; pi <= 10; pi++)
	    {
	      double p = pi / 10.0;
	      std::cout << x << " " << n << " " << p << "\t" 
			<< dbinom( x , n , p , false ) << "\t"
			<< Statistics::dbinom( x , n , p ) << "\n";	      
	    }
    }

}
