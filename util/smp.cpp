#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <stdio.h>
#include <stdlib.h>



#include "char_tok.h"

int main( int argc , char ** argv )
{
  
  if ( argc != 3 ) 
    {
      std::cerr << "expecting ./smp matrix.dat test.set \n" ;
      exit(1);
    }

  // simple utility to read in    
  //   a) a T string vector of test types (e.g. BURDEN UNIQ)

  //   b) a G by 1 + (R+1)*T , of gene-name, T original stats, R*T permuted stats (grouped by T1/REP1 T2/REP1 T1/REP2 T2/REP2 ... )
  //      i.e. output from --dump-null-matrix
  
  //   c) a .set file, that contain rows of set-T pairs
  
  // task evaluate sum-of-statistics in the original versus the null data

  // assume T is genic, and R not tremendously high (so should fit in memory okay)
  
  
  // read first row of (n+1)T original tests (where n is # of test types (i.e. BURDEN, UNIQ, etc)


  std::ifstream NMAT( argv[1] );
  if ( ! NMAT.good() ) 
    {
      std::cerr << "problem opening " << argv[1] << "\n";
      exit(1);
    }
  
  
  // Get test names from first row
  int ntest_types;
  std::vector<std::string> test_names;
  std::string line;
  std::getline( NMAT , line );
  if ( line == "" ) exit(1);
  char_tok tok( line , &ntest_types , '\t' );

  // first obs should be # perms
  --ntest_types;
  int nrep = atoi( tok(0) );
  
  std::cerr << "read " << nrep << " replicates\n";
  std::cerr << "read " << ntest_types << " test types\n";

  for (int tt=0;tt<ntest_types;tt++) 
    {
      test_names.push_back( tok(tt+1) );
    }
  
  //
  // Get all data 
  // G, T original vales, R*T null values
  //
  
  std::map<std::string,int> orig_elem;
  std::vector< std::vector<double> > orig_val( ntest_types );    // TEST * ELEM 
  std::vector< std::vector<std::vector<double> > > null_val( nrep );  // REP * TEST * ELEM 
  for (int r = 0; r < nrep; r++) null_val[r].resize( ntest_types );

  int nelems = 0;

  while ( ! NMAT.eof() )
    {

      std::string elem;

      NMAT >> elem;

      //      std::cout << "elem = " << elem << "\n";

      if ( elem == "" ) continue;
      
      // track 'slot' for this element 
      orig_elem[ elem ] = nelems;

      
      for (int j=0;j<ntest_types;j++)
	{
	  double x;
	  NMAT >> x;
	  //	  std::cout << "orig = " << x << "\n";
	  orig_val[j].push_back( x );
	}

      // now read all reps
      for (int r=0;r<nrep;r++)
	for (int j=0;j<ntest_types;j++)
	  {
	    double x;
	    NMAT >> x;
	    //	    std::cout << "rep, test " << r << " " << j << " = " << x << "\n";
	    null_val[r][j].push_back( x );
	  }
      
      ++nelems;

    }

  NMAT.close();
  
  std::cerr << "read " << nelems << " elements\n";


  //
  // read in sets to test
  //

  std::ifstream IN1( argv[2] );
  if ( ! IN1.good() ) 
    {
      std::cerr << "problem opening " << argv[2] << "\n";
      exit(1);
    }


  std::map<std::string,std::vector<std::string> > sets;
  
  while ( ! IN1.eof() ) 
    {
      // format GENE  SET  (i.e. same as INRICH)
      std::string s, t;
      IN1 >> t >> s;
      if ( s != "" ) sets[s].push_back(t);      
    }
  IN1.close();

  std::cerr << "read " << sets.size() << " sets\n";


  const int ntests = sets.size();

  std::map<std::string,std::vector<std::string> >::iterator ii = sets.begin();
  while ( ii != sets.end() )
    {

      //      std::cout << "considering " << ii->first << "\n";

      // for each test type
      for (int j=0;j<ntest_types;j++)
	{
	  
	  double s = 0;
	  std::vector<double> ns( nrep , 0 );

	  std::vector<std::string> & elems = ii->second;
	  
	  for (int e=0;e<elems.size();e++)
	    {

	      std::map<std::string,int>::iterator ee = orig_elem.find( elems[e] );
	      
	      if ( ee != orig_elem.end() ) 
		{
		  // original
		  s += orig_val[j][ ee->second ];
		  
		  // nulls 
		  for (int r=0;r<nrep;r++)
		    ns[r] += null_val[r][j][ ee->second ]; 
		  
		}
	    }
	  
	  // calculate empirical p-value
	  int pv = 1;
	  for (int r=0; r<nrep; r++ ) 
	    if ( ns[r] >= s ) ++pv;
	  
	  // output

	  std::cout << ii->first << "\t" 
		    << test_names[j] << "\t"
		    << (double)pv / (double)(nrep+1)   // both denom and numer inc. the +1 alreadt
		    << "\n";
	  
	  
	} // next test
      
      ++ii;
    } // next set


  exit(0);

}
