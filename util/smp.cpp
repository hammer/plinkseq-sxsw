#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cstdio>
#include <cstdlib>

#include "pseq.h"
#include "char_tok.h"

double mean( const Data::Vector<double> & x ) ;
double variance( const Data::Vector<double> & , double ) ; 

Data::Vector<double> residualise( const Data::Vector<double> & y , const Data::Vector<double> & x , double vx , double mx );

void read_sets( const std::string & sets_filename , 
		const std::map<std::string,int> & orig_elem ,
		long int * elem_cnt , 
		std::vector<std::vector<int> > & , 
		std::vector<std::string> & , 
		std::vector<std::vector<int> > & ); 


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
  //     E'_i =  E_i * Sum_{ all genes in any set containing i (1 - r^2_ij ) E_j * min(Set_SIZE) 


  // TODO -- multiple testing correction for 
  //         for set-based tests?  
  
  //      -- and set-weighted gene-based tests.
  
  // 'Set-weighted gene-based tests' 
  
  


  
  std::ifstream NMAT( data_filename.c_str() );
  if ( ! NMAT.good() ) 
    {
      std::cerr << "problem opening " << data_filename << "\n";
      exit(1);
    }
  
  
  //
  // Get test names and number of replicates from first row
  //

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
  
  std::map<std::string,int> gene2slot;        // map gene-name to column in E matrix
  std::vector<std::string> slot2gene; 

  std::vector< Data::Matrix<double> > E( ntest_types );  // For each test type, E is rep X gene matrix
  
  //    ---->  Genes 
  //
  // |
  // |
  // V Reps
  
  
  //
  // Read in all data
  //

  int nelems = 0;

  while ( ! NMAT.eof() )
    {

      std::string elem;

      NMAT >> elem;

      if ( elem == "" ) continue;
      
      // track 'slot' for this element 
      gene2slot[ elem ] = nelems;
      slot2gene.push_back( elem );

      std::vector< Data::Vector<double> > dt( ntest_types );
      for (int j=0;j<ntest_types;j++)
	{
	  double x;
	  NMAT >> x;	  
	  dt[j].push_back( x );
	}

      // now read all reps
      for (int r=0;r<nrep;r++)
	for (int j=0;j<ntest_types;j++)
	  {
	    double x;
	    NMAT >> x;	    
	    dt[j].push_back( x );
	  }
    
      // add to main matrix      
      for (int j=0;j<ntest_types;j++)
	E[j].add_col( dt[j] );
      
      ++nelems;

    }

  NMAT.close();
  
  std::cerr << "read " << nelems << " elements\n";



  //
  // read in sets to test
  //

  long int elem_cnt = 0;

  std::vector<std::string> set_names;  
  std::vector<std::vector<int> > sets;      // for set i, all genes
  std::vector<std::vector<int> > gene2sets; // for gene i, all sets


  read_sets( sets_filename , gene2slot , &elem_cnt , sets , set_names , gene2sets );
  
  const int ntests = sets.size();

  //
  // Simple gene-based (element-wise) empirical p-values, and experiment-wide corrected
  //
  
  // max test value per replicate, for each test

  std::vector<std::vector<double> > mx( ntest_types );  // test x rep(+1)

  for (int j=0;j<ntest_types;j++)
    {
      Data::Matrix<double> & T = E[j];
      std::vector<double> & x = mx[j];
      x.resize( nrep+1 ); // includes original data in [0]
      for (int r = 0 ; r <= nrep ; r++ ) 
	{
	  x[r] = T(r,0);
	  for (int e = 1 ; e < nelems ; e++ ) 
	    if ( T(r,e) > x[r] ) x[r] = T(r,e);
	}
    }
  
  // Now for each gene, get pointwise and corrected p-values
  
  for (int e = 0 ; e < nelems ; e++ ) 
    {
      for (int j=0;j<ntest_types;j++)
	{
	  int p1 = 1;
	  int p2 = 1;
	  Data::Matrix<double> & T = E[j];
	  std::vector<double> & x = mx[j];

	  for (int r = 1 ; r <= nrep ; r++ ) 
	    {
	      if ( T(r,e) >= T(0,e) ) ++p1;
	      if ( x[r]   >= T(0,e) ) ++p2;
	    }
	  
	  std::cout << "ELEM\t"
		    << test_names[ j ] << "\t" 		    
		    << (double)p1 / (double)( nrep + 1 ) << "\t"
		    << (double)p2 / (double)( nrep + 1 ) << "\t"
		    << slot2gene[e] << "\n"; 
	}
    }



  //
  // Tests set enrichment
  //

  for (int i = 0 ; i < sets.size(); i++ ) 
    {
      
      // for each test type
      
      for (int j=0;j<ntest_types;j++)
	{
	  
	  std::vector<double> ns( nrep + 1 , 0 );
	  std::vector<int> & elems = sets[i];
	  
	  // Sum over set elements
	  
	  for (int e=0;e<elems.size();e++)
	    {
	      // orig, followed by nulls 
	      for (int r=0;r<=nrep;r++)
		ns[r] += E[j](r,elems[e]);
	    }
	  
	  
	  // Calculate empirical p-value
	  
	  int pv = 1;
	  const double & s = ns[0];
	  for (int r=1; r <= nrep; r++ ) 
	    if ( ns[r] >= s ) ++pv;
	  

	  // Output
	  
	  std::cout << "SET\t" 
		    << test_names[j] << "\t"
		    << (double)(pv+1) / (double)(nrep+1) << "\t" 
		    << elems.size() << "\t"
		    << set_names[ i ] 
		    << "\n";
	  
	  
	  
	} // next test
      

    } // next set
  

  



  //
  // ---------------------- Gene-based re-scoring ------------------------------
  //
  
  // Now for each gene, get pointwise and corrected p-values
  
   for (int e = 0 ; e < nelems ; e++ ) 
     {
       	 
       for (int j=0;j<ntest_types;j++)
 	{

	  // original score * Sum_S ( Sum_(-E) E_res ) 	      
	  // sum over sets this element belongs to 	  
	  
	  // get mean/var of gene 'e'
	  Data::Vector<double> & evec = E[j].col(e);
	  double meane = mean( evec );
	  double vare = variance( evec , meane );
	  
	  // consider each set to which this guy belongs
	  // and create a matrix of residuals (excluding 'e')
	  
	  // For now, simply sum all Sets  (only adjust for correl between index element 
	  // and other genes in a fellow set, so we can have just a single accumulator)
	  // outside the set loo
	  
	  // Weight by inverse of smallest set size of set for the index and the fellow
	  //  then again, does it not mean more for a set w/ >1 gene to be significant too?
	  //  i.e. you might not want set size == 2 to contibute a lot, if it is only a single 
	  //       other gene contributing...

	  std::set<int> fellows;
	  std::map<int,double> wgt;
	  
	  for (int s=0;s<gene2sets[e].size(); s++ ) 
	    {	      
	      std::vector<int> & inset = sets[ gene2sets[e][s] ];	      
	      for (int f = 0 ; f < inset.size(); f++ ) 
		{
		  if ( inset[f] == e ) continue; // skip if same gene
		  fellows.insert( inset[f] );

		  // minimum set-size weight
		  if ( wgt[ inset[f] ] < 1.0/(double)inset.size() ) 
		    {
		      wgt[ inset[f] ] = 1.0/(double)inset.size() ;
		    }
		}
	    }
	  
	  std::vector<double> accum( nrep+1 , 0 );
	  
	  std::set<int>::iterator fi = fellows.begin();
	  while ( fi != fellows.end() )
	    {
	      Data::Vector<double> resid = residualise( E[j].col( *fi ) , evec , vare , meane ) ;
	      for (int r=0;r<nrep;r++) accum[r] += wgt[ *fi ] * resid[r];
	      ++fi;
	    }
	  
	  

	  //
	  // Get empirical p-value
	  //
	  
	  int pv1 = 1;  // orig * fellow-score 
	  int pv2 = 1;  // fellow-score
	  
	  // has at least one fellow?
	  if ( fellows.size() > 0 ) 
	    {
	      for (int r=1;r<=nrep;r++)
		{
		  if ( accum[r] * E[j](r,e) >= accum[0] * E[j](0,e) ) ++pv1;
		  if ( accum[r] >= accum[0] ) ++pv2;	      
		}

	      // Output
	      
	      std::cout << "REVISED\t" 
			<< test_names[j] << "\t"
			<< gene2sets[e].size() << "\t"
			<< fellows.size() << "\t"			
			<< (double)(pv1) / (double)(nrep+1) << "\t" 
			<< (double)(pv2) / (double)(nrep+1) << "\t" 
			<< slot2gene[e] << "\n";
	      
	      
	    }
	  else
	    {
	      for (int r=1;r<=nrep;r++)
		{
		  if ( E[j](r,e) >= E[j](0,e) ) ++pv1;
		}
	      
	      // Output
	      
	      std::cout << "REVISED\t" 
			<< test_names[j] << "\t"
			<< gene2sets[e].size() << "\t"
			<< fellows.size() << "\t"			
			<< (double)(pv1) / (double)(nrep+1) << "\t" 
			<< "NA" << "\t"
			<< slot2gene[e] << "\n";	      
	      
	    }

	}
     }	  


  exit(0);

}




void read_sets( const std::string & sets_filename , 
		const std::map<std::string,int> & gene2slot ,  
		long int * elem_cnt , 
		std::vector<std::vector<int> > & sets , 
		std::vector<std::string> & set_names , 
		std::vector<std::vector<int> > & gene2sets ) 

{
  
  std::map<std::string,int> setmap;
  gene2sets.resize( gene2slot.size() );

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
      std::map<std::string,int>::const_iterator ee = gene2slot.find( t );
      if ( ee != gene2slot.end() ) 
	{
	  
	  std::map<std::string,int>::iterator si = setmap.find( s ) ;
	  int sidx;
	  if ( si == setmap.end() ) 
	    {
	      sidx = setmap.size();
	      setmap[ s ] = sidx;
	      std::vector<int> dummy;
	      sets.push_back( dummy );
	      set_names.push_back( s );
	    }
	  else 
	    sidx = setmap[ s ];
	  
	  sets[ sidx ].push_back( ee->second );
	  gene2sets[ ee->second ].push_back( sidx );
					    
	  ++(*elem_cnt);
	}
    }
  IN1.close();

  std::cerr << "read " << *elem_cnt << " elements in " << sets.size() << " sets\n";

  return;
}


double mean( const Data::Vector<double> & x ) 
{
  double m = 0;
  const int n = x.size();
  for (int i=0;i<n;i++) { m += x[i]; }
  return m / (double)n;
}

double variance( const Data::Vector<double> & x , double mx )
{
  double sx = 0;
  const int n = x.size();
  for (int i=0;i<n;i++) { sx += ( x[i] - mx ) * ( x[i] - mx ) ; } 
  return sx/((double)n-1.0);
}

Data::Vector<double> residualise( const Data::Vector<double> & y , const Data::Vector<double> & x , double vx , double mx )
{

  // b_yx = cov(y,x)/var(x)
  // return y - b_yx * x 

  double my = mean(y);
  double sxy = 0;
  const int n = y.size();
  for (int i=0;i<n;i++) sxy = ( y[i] - my ) * ( x[i] - mx );
  sxy /= (double)(n-1.0);
  double beta = sxy / vx;  
  Data::Vector<double> r(n);
  for (int i=0;i<n;i++) r[i] = y[i] - beta * x[i];
  return r;
}

