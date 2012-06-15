#ifndef __PSEQ_DAVIES_H__
#define __PSEQ_DAVIES_H__

#include <vector>

extern "C" { 
void davies_qfc(double* lb1, 
		double* nc1, 
		int* n1, 
		int *r1, 
		double *sigma, 
		double *c1, 
		int *lim1, 
		double *acc, 
		double* trace, 
		int* ifault, 
		double *res);
}

class Davies { 

 public:

  Davies() 
    {

    }

  // for 'j' non-central chi-squared random variables:

  double pvalue( double q , 
		 std::vector<double> & lb , 		 
		 std::vector<int> & df , 
		 bool * okay , 
		 std::vector<bool> * ifaults = NULL )
  {
    
    //
    
    /*   lambdas=as.double(lambda),         lb[j]            coefficient of j-th chi-squared variable	  */
    /*   noncentral=as.double(delta),	     nc[j]            non-centrality parameter			  */
    /*   df=as.integer(h),		     n[j]             degrees of freedom     			  */
    /*     r=as.integer(r),		     lb.size() */
    /*   sigma=as.double(sigma),	     sigma            coefficient of standard normal variable	  */
    /*   q=as.double(q),		     c                point at which df is to be evaluated	  */
    /*   lim=as.integer(lim),		     lim              maximum number of terms in integration	  */
    /*   acc=as.double(acc),		     acc              maximum error				  */
    /*   trace=as.double(rep(0,7)),	     trace[7] */
    /*   ifault=as.integer(0),		     ifault[6]         */
    /*   res=as.double(0),                  res              double * result; */
    
    int r = lb.size();
    std::vector<double> nc(r,0); // central chi-sq under null
    double sigma = 0;     //
    double acc = 0.000001; // maximum error 
    int    lim   = 10000; // maximum number of terms in integration
    double trace[7];  
    int ifault[6];
    double res = 0; //result;
    
    davies_qfc( &(lb[0]) , 
		&(nc[0]) , 
		&(df[0]) , 
		& r , 
		& sigma , 
		& q , 
		& lim , 
		& acc ,        
		&(trace[0]) , 
		&(ifault[0]) , 
		&res );
    
    
    if ( ifault[0] || ifault[2] || ifault[3] || ifault[4] ) 
      {
	*okay = false;
	if ( ifaults ) 
	  {
	    ifaults->resize( 6 );
	    for (int i=0;i<6;i++)
	      (*ifaults)[i] = ifault[i]; 
	  }
      }
    else 
      *okay = true;

    // ifault = 1       required accuracy NOT achieved
    //          2       round-off error possibly significant
    //          3       invalid parameters
    //          4       unable to locate integration parameters
    //          5       out of memory
    
    // trace[0]         absolute sum
    // trace[1]         total number of integration terms
    // trace[2]         number of integrations
    // trace[3]         integration interval in final integration
    // trace[4]         truncation point in initial integration
    // trace[5]         s.d. of initial convergence factor
    // trace[6]         cycles to locate integration parameters 
    
    return 1 - res;    
    
    // Outputs
    
  }

 private:
  
};

#endif
