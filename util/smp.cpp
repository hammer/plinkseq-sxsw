#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cstdio>
#include <cstdlib>


#include "pseq.h"
#include "char_tok.h"

std::map<std::string,std::vector<std::string> > read_sets( const std::string & sets_filename , 
							   const std::map<std::string,int> & orig_elem ,
							   long int * elem_cnt );

int main( int argc , char ** argv )
{
  
  if ( argc != 4 ) 
    {
      std::cerr << "usage:   ./smp list   sets.list data.mat\n"
		<< "         ./smp matrix sets.mat  data.mat\n"
		<< "         ./smp make   sets.list\n";

      exit(1);
    }

  std::string run_mode = argv[1];
  std::string sets_filename = argv[2];
  std::string data_filename = argv[3];
    

  // simple utility to read in    
  //   a) a T string vector of test types (e.g. BURDEN UNIQ)

  //   b) a G by 1 + (R+1)*T , of gene-name, T original stats, R*T permuted stats (grouped by T1/REP1 T2/REP1 T1/REP2 T2/REP2 ... )
  //      i.e. output from --dump-null-matrix
  
  //   c) a .set file, that contain rows of set-T pairs
  
  // task evaluate sum-of-statistics in the original versus the null data

  // assume T is genic, and R not tremendously high (so should fit in memory okay)
  
  // read first row of (n+1)T original tests (where n is # of test types (i.e. BURDEN, UNIQ, etc)

  
  //  Notation: 
  
  //    E_i  is element (i.e. typically genic score element i)
  //    S_j  is additive set score = Sum_{'i' in 'j'} E_i 
  
  //    Evaluate empirical null distribution of S_j by application to R sets of E calculated under the null
  //    (phenotype shuffling in PSEQ, that preserves LD between tests)

  //    2) Obtain per-gene revised scores, given any set-based enrichment: score for element 'i' 
  //       (performed within test-type currently)
  //      
  //     E'_i =  E_i * Sum_{ all sets containing i } Sum_{ j=all other genes in set } (1 - r^2_ij ) E_j 

  //  Note: consider a version of this where (to alleviate concerns of over-lapping sets) we perform QR 
  //   decomposition on the original (implicit) S * G matrix (from 0s, 1s) to produce a new 'continuous' 
  //   set of sets that a) are orthogonal, b) span all genes. We then proceed with scoring these sets for the
  //   purpose of collapsing back for a per-gene score. 



  
  std::ifstream NMAT( data_filename.c_str() );
  if ( ! NMAT.good() ) 
    {
      std::cerr << "problem opening " << data_filename << "\n";
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

  long int elem_cnt = 0;
  
  std::map<std::string,std::vector<std::string> > sets = read_sets( sets_filename , orig_elem, &elem_cnt );

  const int ntests = sets.size();

  std::map<std::string,std::vector<std::string> >::iterator ii = sets.begin();
  while ( ii != sets.end() )
    {
      
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
	  
	  std::cout << test_names[j] << "\t"
		    << (double)pv / (double)(nrep+1) << "\t"  // both denom and numer inc. the +1 alreadt
		    << ii->second.size() << "\t"
		    << ii->first  		    
		    << "\n";
	  
	  
	} // next test
      
      ++ii;
    } // next set


  exit(0);

}




std::map<std::string,std::vector<std::string> > read_sets( const std::string & sets_filename , 
							   const std::map<std::string,int> & orig_elem , 						   
							   long int * elem_cnt )
{
  
  std::map<std::string,std::vector<std::string> > sets;

  std::ifstream IN1( sets_filename.c_str() );
  if ( ! IN1.good() ) 
    {
      std::cerr << "problem opening " << sets_filename << "\n";
      exit(1);
    }

  while ( ! IN1.eof() ) 
    {
      // format GENE /tab SET (blurb)  (i.e. same as INRICH)
      
      std::string line;
      std::getline( IN1 , line );
      if ( line == "" ) continue;
      int n;
      char_tok tok( line , &n , '\t' );
      if ( n < 2 ) { std::cerr << "problem with format of set file\n"; exit(1); }
      std::string t = tok(0);
      std::string s = tok(1);
      
      // only load in elements found in this result set
      std::map<std::string,int>::const_iterator ee = orig_elem.find( t );
      if ( ee != orig_elem.end() ) 
	{
	  sets[ s ].push_back( t );
	  ++(*elem_cnt);
	}
    }
  IN1.close();

  std::cerr << "read " << elem_cnt << " elements in " << sets.size() << " sets\n";
  
  return sets;
}

