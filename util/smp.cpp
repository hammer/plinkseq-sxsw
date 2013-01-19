#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cstdio>
#include <cstdlib>

#include "plinkseq.h"

#include "char_tok.h"

double mean( const Data::Vector<double> & x ) ;
double variance( const Data::Vector<double> & , double ) ; 

Data::Vector<double> residualise( const Data::Vector<double> & y , const Data::Vector<double> & x , double vx , double mx );

void read_sets( const std::string & sets_filename , 
		const std::map<std::string,int> & orig_elem ,
		long int * elem_cnt , 
		std::vector<std::vector<int> > & , 
		std::vector<std::string> & , 
		std::vector<std::string> & , 
		std::vector<std::vector<int> > & ); 


struct elem_t 
{
  elem_t( double x , const std::string & n ) : name(n) , statistic(x) { } 
  std::string name;
  double statistic;
  bool operator< ( const elem_t & rhs ) const
  {
    // sort in reverse order
    if ( statistic > rhs.statistic ) return true;
    if ( statistic < rhs.statistic ) return false;
    return name < rhs.name;
  }
};

struct anon_elem_t 
{
  // sorted, but kept unique too
  anon_elem_t( double x , int p ) : pos(p) , statistic(x) { } 
  int pos;
  double statistic;
  bool operator< ( const anon_elem_t & rhs ) const
  {
    // sort in reverse order
    if ( statistic > rhs.statistic ) return true;
    if ( statistic < rhs.statistic ) return false;
    return pos < rhs.pos;
  }
};


int main( int argc , char ** argv )
{

  
  if ( argc != 5 && argc != 6 ) 
    {
      std::cerr << "usage:   ./smp {-m} gene.list equiv.list sets.list data.mat\n";
      exit(1);
    }

  int off = argc == 6 ? 1 : 0;
  bool calc_max = false;
  bool calc_gene = false;
  if ( off ) 
    {
      std::string t = argv[1];
      if ( t == "-m" ) calc_max = true;
      if ( t == "-g" ) calc_gene = true;
      if ( t == "-mg" || t == "-gm" ) calc_max = calc_gene = true;
    }

  std::string genes_filename = argv[1+off];
  std::string equiv_filename = argv[2+off];
  std::string sets_filename  = argv[3+off];
  std::string data_filename  = argv[4+off];

  

  // simple utility to read in    
  //   a) a T string vector of test types (e.g. BURDEN UNIQ)

  //   b) a G by 1 + (R+1)*T , of gene-name, T original stats, R*T permuted stats (grouped by T1/REP1 T2/REP1 T1/REP2 T2/REP2 ... )
  //      i.e. output from --dump-null-matrix
  
  //   c) a .set file, that contain rows of set-T pairs
  
  // task evaluate sum-of-statistics in the original versus the null data

  // assume T is genic, and R not tremendously high (so should fit in memory okay)
  
  // read first row of (n+1)T original tests (where n is # of test types (i.e. BURDEN, UNIQ, etc)

  // Also: for each set, get the ranked ordering and get p-value cut-offs therein per gene


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
  
  std::map<std::string,std::string> gene_annot;
  
  std::ifstream GENES( genes_filename.c_str() , std::ios::in );
  
  if ( GENES.good() ) 
    {
      while ( ! GENES.eof() ) 
	{
	  std::string line;
	  std::getline( GENES , line );
	  if ( line == "" ) continue;
	  int ncol;
	  char_tok tok( line , &ncol , '\t' );
	  if ( ncol != 2 ) 
	    { 
	      std::cerr << "skipping line in " << genes_filename 
			<< " that doesn't have 2 tab-delim fields\n"; 
	      continue; 
	    } 
	  gene_annot[ tok(0) ] = tok(1);
	}      
      std::cerr << "read " << gene_annot.size() << " gene/anotations\n";
    }
  else
    std::cerr << "skipping a gene-list\n";

  GENES.close();


 
  //
  // data
  //

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
  
  std::map<std::string,int> gene2slot;         // map gene-name to column in E matrix
  std::vector<std::string> slot2gene; 

  std::map<std::string,std::string> gene2result; 
  std::map<int,int> agenecnt;
  std::map<int,int> ugenecnt;

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
      
      // altA, altU, nA, nU
      // (ignore for now)

      double alta, altu;
      double na, nu;

      NMAT >> alta >> altu >> na >> nu;

      gene2result[ elem ] = Helper::int2str( int( alta * na + 0.5 ) ) + "/" + Helper::int2str( int (altu * nu + 0.5 ) );      
            
      // do we want this gene?
      if ( gene_annot.size() > 0 ) 
	{
	  if ( gene_annot.find( elem ) == gene_annot.end() ) 
	    {
	      for (int j=0;j<ntest_types;j++)
		{
		  double x;
		  NMAT >> x;	  		  
		}

	      for (int r=0;r<nrep;r++)
		for (int j=0;j<ntest_types;j++)
		  {
		    double x;
		    NMAT >> x;	    
		  }
	      continue;
	    }
	}
      
      // track 'slot' for this element 
      gene2slot[ elem ] = nelems;
      slot2gene.push_back( elem );      

      // track cnts
      agenecnt[ nelems ] = int( alta * na + 0.5 );
      ugenecnt[ nelems ] = int( altu * nu + 0.5 );


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
  // Given actually observed genes, read and create equivalence map
  //
  
  // we've already populated index 'i' for 'gene' name
  //gene2slot[ gene ] = i;
  //slot2gene[ i ] = gene;


  //
  // equivalence sets for genes
  //

  // for gene 'i', store all equiv. genes in vector 
  std::map<int,std::vector<int> > eq;

  if ( ! ( equiv_filename == "-" || equiv_filename == "." ) )
    {
      int eset = 0, fset = 0;
      std::ifstream EQ( equiv_filename.c_str() );
      while ( ! EQ.eof() )
	{
	  // each row should be a complete group
	  std::string line;
	  std::getline( EQ , line );
          if ( line == "" ) continue;
          int ncol;
          char_tok tok( line , &ncol , '\t' );
	  if ( ncol > 1 ) 
	    {
	      ++fset;
	      std::vector<int> x;
	      for (int i=0;i<ncol;i++) 
		{
		  std::map<std::string,int>::iterator ii = gene2slot.find( tok(i) );
		  if ( ii != gene2slot.end() ) x.push_back( ii->second );
		}
	      
	      if ( x.size() > 1 ) 
		{
		  for (int i=0;i<x.size();i++)
		    for (int j=0;j<x.size();j++)
		      if ( i != j ) eq[x[i]].push_back(x[j]);
		  ++eset;
		}
	    }
	}
      std::cerr << "read " << fset << " gene equivalence sets, " << eset << " of which map to observed genes\n";
      EQ.close();
    }
  


  //
  // read in sets to test
  //

  long int elem_cnt = 0;

  std::vector<std::string> set_names;  
  std::vector<std::string> set_descriptions;
  std::vector<std::vector<int> > sets;      // for set i, all genes
  std::vector<std::vector<int> > gene2sets; // for gene i, all sets
  

  read_sets( sets_filename , gene2slot , &elem_cnt , sets , set_names , set_descriptions, gene2sets );
  
  const int ntests = sets.size();
  
  
  //
  // Note which genes are in equivalence sets
  //
  
  std::set<int> ineq; 
  std::map<int,std::vector<int> >::iterator ee = eq.begin();
  while ( ee != eq.end() ) { ineq.insert( ee->first ); ++ee; }
  
  
  //
  // Tests set enrichment
  //

  for (int i = 0 ; i < sets.size(); i++ ) 
    {

      std::vector<int> & elems = sets[i];
      
      // make a quick store of the genes in this set
      std::set<int> inset;
      for (int e=0;e<elems.size();e++) inset.insert(elems[e]);
      
      //
      // for each test type
      //

      for (int j=0;j<ntest_types;j++)
	{

	  // basic statistic;
	  std::vector<double> ns( nrep + 1 , 0 );

	  // count number of actual independent elements that will be added:
	  
	  std::set<elem_t> indep_elem;
	  std::vector<std::set<anon_elem_t> > ordered ( nrep );
	  
	  //
	  // Sum over independent set elements
	  //
	  
	  int acnt = 0 , ucnt = 0;

	  for (int e=0;e<elems.size();e++)
	    {
	     
	      bool in_equiv_set = ineq.find( elems[e] ) != ineq.end() ;

	      if ( ! in_equiv_set )
		{
		  
		  // for independent elements, just add own value; 
		  // orig, followed by nulls 
		  
		  for (int r=0;r<=nrep;r++)
		    ns[r] += E[j](r,elems[e]);
		  
		  // store original 
		  acnt += agenecnt[ elems[e] ];
		  ucnt += ugenecnt[ elems[e] ]; 
		    
		  if ( calc_max )
		    {
		      // track original statistics for this set, this will order by original statistic also
		      indep_elem.insert( elem_t( E[j](0,elems[e]) , slot2gene[elems[e]] ) );
		      
		      // track
		      for (int r=1;r<=nrep;r++)
			{
			  int pos = ordered[r-1].size();
			  ordered[r-1].insert( anon_elem_t( E[j](r,elems[e]) , pos ) );
			}		  
		    }
	
		}	      
	      else
		{
	
	  // otherwise, if in an equivalence set, add only the max from the set;
		  // considering only elements that are actually in the set
		  
		  std::vector<int> & t = eq[ elems[e] ];
		  std::vector<int> x;
		  for (int ee=0;ee<t.size();ee++)
		    {
		      if ( inset.find( t[ee] ) != inset.end() ) 
			{
			  x.push_back( t[ee] );
			}		      
		    }

		  
		  for (int r=0;r<=nrep;r++)
		    {

		      // Current value; only add this if it is the max.
		      
		      double mx = E[j](r,elems[e]);
		      int    mxi = elems[e];
		      
		      // if ( r == 0 ) 
		      // 	std::cerr << "considering " << set_names[i] << " " << slot2gene[ mxi ] << " " << mx << "\n";
		      
		      // so, see if any other equiv, in-set element scores higher.
		      // for ties, always take the lower element number; that way
		      // we will avoid double-counting equivalently-scored equivalent elements
		      
		      bool add_this = true;

		      for (int f=0;f<x.size();f++)
			{
			  
			  // is score bigger (or tied, but lower-indexed element)?
			  
			  if ( E[j](r,x[f]) > mx ) 
			    { 
			      add_this = false; 
			      break; 
			    }

			  if ( E[j](r,x[f]) == mx && x[f] < mxi ) 
			    { 
			      add_this = false; 
			      break; 
			    }
			}
		      
		      if ( add_this )
			{
			  ns[ r ] += mx;
			  
			  // store original 
			  if ( r== 0 )
			    {
			      acnt += agenecnt[ mxi ];
			      ucnt += ugenecnt[ mxi ]; 
			    }
			  
			  // track null replicates only here
			  if ( calc_max )
			    {
			      if ( r ) 
				{
				  int pos = ordered[r-1].size();
				  ordered[r-1].insert( anon_elem_t( mx , pos ) );
				}
			      else indep_elem.insert( elem_t( mx , slot2gene[ mxi ] ) );  // and originals w/ labels
			    }
			}
		      
		    }	  
		  
		}
	    }

	  
	  //
	  // Calculate empirical p-value
	  //
	  
	  int pv = 1;
	  const double & setstat = ns[0];
	  for (int r=1; r <= nrep; r++ ) 
	    if ( ns[r] >= setstat ) ++pv;
	  
	  
	  //
	  // Output
	  //

	  std::cout << "SET\t" 
		    << test_names[j] << "\t"
		    << (double)(pv) / (double)(nrep+1) << "\t" 
		    << elems.size() << "\t"
		    << acnt << "/" << ucnt << "\t";

	  if ( ucnt > 0 ) 
	    std::cout << (double)acnt / (double)ucnt << "\t";
	  else 
	    std::cout << ".\t";
	  
	  std::cout << set_names[ i ] 
		    << "\n";
	  
	  
	  
	  //
	  // Set-test version 2: ordered ranked cumulative sums
	  //
	  
	  if ( calc_max )
	    {

	      // For original, determine order of results (Num Indep Elements)
	      
	      int nie = indep_elem.size();
	      std::vector<double> origmx( nie , 0 );
	      std::vector<double> origmxr( nie , 0 ); // reverse stat
	      std::set<elem_t>::iterator ii = indep_elem.begin();
	      std::set<elem_t>::reverse_iterator rr = indep_elem.rbegin();
	      int ix = 0;
	      while ( ii != indep_elem.end() ) 
		{
		  // make cumulative sum
		  if ( ix ) origmx[ ix ] = origmx[ ix-1 ] + ii->statistic ;
		  else origmx[ ix ] = ii->statistic ;
		  
		  // reverse statistic
		  if ( ix ) origmxr[ ix ] = origmxr[ ix-1 ] + rr->statistic ;
		  else origmxr[ ix ] = rr->statistic ;

		  ++ix;
		  ++ii;
		  ++rr;
		}
	      
	      
	      // make ordered[] sets cumulative
	      std::vector<int> pvalmx( nie , 0 );  // for best subset
	      std::vector<int> pvalmx2( nie , 0 ); // pvalue for remainder
	      for (int r=1;r<nrep;r++)
		{
		  double s = 0;  // best N
		  double rs = 0; // reverse statistic (taking worst N)

		  int ix = 0;
		  std::set<anon_elem_t>::iterator ii = ordered[r-1].begin();
		  std::set<anon_elem_t>::reverse_iterator rr = ordered[r-1].rbegin();
		  while ( ii != ordered[r-1].end() )
		    {
		      // 0..k  statistic
		      s += ii->statistic;
		      if ( s >= origmx[ ix ] ) pvalmx[ ix ]++;
		      
		      // k+1..n statistic (i.e. any signal in the remainder?)
		      rs += rr->statistic;
		      if ( rs >= origmxr[ ix ] ) pvalmx2[ ix ]++;

		      ++ix;
		      ++ii;
		      ++rr;
		    }
		}
	      
	      
	      // print output
	      ii = indep_elem.begin();	  
	      for (int m=1;m<=nie;m++)
		{
                  
		  //std::cout << "origmx = " << origmx[m-1] << " " << ns[0] << "\n";
		  
		  std::cout << "MX\t" 
			    << test_names[j] << "\t"
			    << (double)(pvalmx[m-1]+1) / (double)(nrep+1) << "\t" 
			    << (double)(pvalmx2[nie-m]+1) / (double)(nrep+1) << "\t" 
			    << m << "\t"
			    << ( ns[0] > 0 ? (double)origmx[m-1] / (double)ns[0] : 1.0 ) << "\t"
			    << ii->name << "\t"
			    << gene2result[ ii->name ] << "\t"
			    << set_names[ i ] 		
			    << "\n";
		  ++ii;
		}
	    }


	  
	} // next test
      
    } // next set
  

  //
  // unless explicitly requested, just do SET level tests
  //

  if ( ! calc_gene ) exit(0);


  //
  // Obtain max(statistic) distribution (per test, per replicate)
  //
  
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
  


  //
  // ---------------------- Gene-based re-scoring ------------------------------
  //
  


  // Now for each gene, get pointwise and corrected p-values
  
   for (int e = 0 ; e < nelems ; e++ ) 
     {

       std::cerr << "testing " << e << " of " << nelems << "        \r";


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
	  // outside the set loop
	  
	  // Weight each fellow inversely by its size of the smallest shared set
	  // (up tp N=1000 genes -- after this point, being in the same set does not count)
	  // Might want to explore other possibilites here, this is somewhat arbitrary...
	  
	  std::set<int> fellows;
	  std::map<int,double> wgt;
	  std::map<int,int> minset;
	  
	  // for each SET 's' that index ELEM 'e' belongs to
	  for (int s=0;s<gene2sets[e].size(); s++ ) 
	    {	      
	      
	      // consider all other ELEM 'f' in these SETS ('friends')

	      std::vector<int> & inset = sets[ gene2sets[e][s] ];	      

	      if ( inset.size() < 1000 ) 
		{
		  double w = 1.0/(double)inset.size();
		  
		  for (int f = 0 ; f < inset.size(); f++ ) 
		    {
		      
		      if ( inset[f] == e ) continue; // skip if same gene
		      
		      fellows.insert( inset[f] );
		      
		      // minimum set-size weight
		      
		      if ( wgt[ inset[f] ] < w )
			{
			  wgt[ inset[f] ] = w;
			  minset[ inset[f] ] =  gene2sets[e][s];
			}
		    }
		}
	    }

	  
	  //
	  // figure out which fellows are contributing significantly to 
	  // this set-weighted test
	  //
	  
	  std::vector<std::vector<double> > fellow_tests( nrep+1 ) ;
	  
	  std::vector<double> accum( nrep+1 , 0 );
	  
	  std::set<int>::iterator fi = fellows.begin();
	  while ( fi != fellows.end() )
	    {
	      Data::Vector<double> resid = residualise( E[j].col( *fi ) , evec , vare , meane ) ;
	      for (int r=0;r<=nrep;r++) 
		{
		  accum[r] += wgt[ *fi ] * resid[ r ];
		  fellow_tests[r].push_back( wgt[ *fi ] * resid[ r ] );
		}
	      
	      ++fi;
	    }
	  

	  //
	  // and then standardize fellow contributions to be a
	  // proportion of the total sum (which means the weighting
	  // should be implicitly counted now
	  //

	  for (int r=0;r<=nrep;r++) 
	    for (int f=0;f<fellows.size();f++)
	      fellow_tests[r][f] /= accum[r];
		
	  //
	  // Get empirical p-value
	  //
	  
	  int pv0 = 1;  // basic, pointwise EMP1
	  int pvm = 1;  // EMP2 for gene
	  int pv1 = 1;  // orig * fellow-score 
	  int pv2 = 1;  // fellow-score
	  	  
	  
	  // has at least one fellow?
	  if ( fellows.size() > 0 ) 
	    {
	      for (int r=1;r<=nrep;r++)
		{
		  if ( E[j](r,e) >= E[j](0,e) ) ++pv0;
		  if ( mx[j][r]  >= E[j](0,e) ) ++pvm;
		  if ( accum[r] * E[j](r,e) >= accum[0] * E[j](0,e) ) ++pv1;
		  if ( accum[r] >= accum[0] ) ++pv2;	      
		}
	      

	      // Output
	      
	      std::cout << "GENE\t" 
			<< test_names[j] << "\t"
			<< (double)(pv0) / (double)(nrep+1) << "\t" 
			<< (double)(pvm) / (double)(nrep+1) << "\t" 
			<< (double)(pv2) / (double)(nrep+1) << "\t" 
			<< (double)(pv1) / (double)(nrep+1) << "\t" 			
			<< gene2sets[e].size() << "\t"
			<< fellows.size() << "\t"			
			<< slot2gene[e] << "\t"
			<< gene2result[ slot2gene[e] ] << "\t";
	      
	      std::string ann = gene_annot[ slot2gene[e] ] ;
	      std::cout << ( ann == "" ? "." : ann ) ;

	      if ( gene2sets[e].size() == 0 ) std::cout << "\t.";	      
	      else for (int s=0;s<gene2sets[e].size();s++)		
		     std::cout << ( s > 0 ? "|" : "\t" ) << set_names[ gene2sets[e][s] ] ;
			  
	      std::cout << "\n";
	      
	      // if the fellow-test is nominally significant, give
	      // some extra output as to what is driving the signal.

	      if ( (double)(pv1) / (double)(nrep+1)  < 0.05 && (pv2) / (double)(nrep+1) < 0.05 ) 
		{

		  std::vector<int> pv( fellows.size() , 1 );
		  for (int r=1;r<=nrep;r++)
		    for (int f=0;f<fellows.size();f++)
		      if ( fellow_tests[r][f] >= fellow_tests[0][f] ) ++pv[f];

		  double totwgt = 0;
		  std::set<int>::iterator fi = fellows.begin();
		  while ( fi != fellows.end() )
		    { totwgt += wgt[ *fi ]; ++fi; } 
		  
		  int f = 0;
		  fi = fellows.begin();
		  while ( fi != fellows.end() )
		    {
		      double pval = (double)pv[f] / (double)(nrep+1);
		  
		      if ( pval < 0.05 ) 
			{

			  std::cout << "FELLOW" << "\t"
				    << slot2gene[e] << "\t"
				    << slot2gene[ *fi ] << "\t"
				    << pval << "\t"
				    << sets[ minset[ *fi ] ].size() << "\t"
				    << wgt[ *fi ] / totwgt << "\t"
				    << set_names[ minset[ *fi ] ] << "\t"
				    << set_descriptions[ minset[ *fi ] ] << "\n";
			  
			}
		      ++f;
		      ++fi;
		    }
		  
		}
	    }
	  else // if no fellows...
	    {
	      for (int r=1;r<=nrep;r++)
		{
		  if ( E[j](r,e) >= E[j](0,e) ) ++pv0;
		  if ( mx[j][r]  >= E[j](0,e) ) ++pvm;
		}
	      
	      // Output
	      
	      std::cout << "GENE\t" 
			<< test_names[j] << "\t"
			<< (double)(pv0) / (double)(nrep+1) << "\t" 
			<< (double)(pvm) / (double)(nrep+1) << "\t" 
			<< "NA" << "\t"
			<< (double)(pv0) / (double)(nrep+1) << "\t" 
			<< gene2sets[e].size() << "\t"
			<< fellows.size() << "\t"			
			<< slot2gene[e] << "\t"
			<< gene2result[ slot2gene[e] ] << "\t";
				      
	      std::string ann = gene_annot[ slot2gene[e] ] ;
	      std::cout << ( ann == "" ? "." : ann ) ;
	      
	      if ( gene2sets[e].size() == 0 ) std::cout << "\t.";
	      else for (int s=0;s<gene2sets[e].size();s++)		
		     std::cout << ( s > 0 ? "|" : "\t" ) << set_names[ gene2sets[e][s] ] ;
	      
	      std::cout << "\n";
	      
	      
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
		std::vector<std::string> & set_descriptions , 
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
	      set_descriptions.push_back( n == 3 ? tok(2) : "." );
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

