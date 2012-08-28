#include "pops.h"
#include "func.h"
#include "plinkseq.h"
#include <cmath>

extern GStore g;

double Pseq::IPop::min_af = 1e-4;

// helper function to sum log-probabilities
double sumLogProbs( const std::vector<double> & logprobs)
{
  double max = 0;
  unsigned int i;
  for (i = 0; i<logprobs.size(); i++) {
    if (i==0 || logprobs[i]>max)
      max = logprobs[i];
  }
  if ( std::isinf(max) ) // the largest probability is 0 (log prob= -inf)
    return max;   // return log 0
  double p = 0;
  for (i = 0; i<logprobs.size(); i++) {
    p += exp(logprobs[i]-max);
  }
  return max + log(p);
}


Pseq::IPop::IPop( const std::string & filename , SeqDBase * seqdb ) : seqdb(seqdb) 
{

  Helper::checkFileExists( filename );
  
  // expecting format lable 
  // row 1 : N pop labels
  // rows 2 +  :   chr:start..end   a1  a2   { N a1 allele frequencies }
  // note -- a1, a2 assumed on +ve strand; not assumed to be major/minor or ref/alt alleles
  //  i.e. p can be >0.5

  InFile F( filename );

  // header
  
  std::string l = F.readLine();
  
  if ( l.size() == 0 ) Helper::halt( "no data to read" );
  
  if ( l[0] == '#' ) l = l.substr(1); // skip comment char
  int npop = 0;
  Helper::char_tok tok( l , &npop , '\t' );
  
  for (int i=0;i<npop;i++)
    {
      pop_t p;
      p.label = tok(i);
      p.prior = 1.0/(double)npop;
      populations.push_back( p );
    }
  
  // variant/freq data

  int nv = 0 ;
  std::set<int> bad_sites;

  while ( ! F.eof() )
    {
      // expecting chr1:100..101  A1  A2  freqs
      
      // where frequencies reported are expected to be for A1

      // Internally, we always want to use ALT freqs; thus we need to check that A2 == REF
      // If not, flip (after checking that A1 matches REF )
      
      // If not SEQDB is attached, print a warning and assume that A1 == ALT and A2 == REF 
      //  (i.e. we are still giving ALT allele frequencies)


      std::string l = F.readLine();
      if ( l == "" ) continue;
      int n2 =0;

      Helper::char_tok tok( l , &n2 , '\t' );
      if ( n2 != npop + 3 ) 
	Helper::halt( "invalid # of fields:\n[" + l + "]" ); 
      
      bool okay = false;
      Region pos( std::string( tok(0) ) , okay );
      
      std::string a1 = tok(1);
      std::string a2 = tok(2);
      
      //
      // Check which is reference : should always be 'A2'
      //

      bool flip = false;
      bool bad_site = false;

      if ( seqdb )
	{
	  std::string sref = seqdb->lookup( pos );

	  Helper::str2upper( sref );

	  if ( sref != a1 && sref != a2 ) 
	    {
	      plog.warn( "mismatch with SEQDB " , std::string(tok(0)) + ", " +
			 a1+"/"+a2+" in panel,  versus " + sref + " in SEQDB" );
	      bad_site = true;
	    }
	  else
	    {
	      if ( a1 == sref ) 
		{
		  std::string tmp = a2;
		  a2 = a1;
		  a1 = tmp;
		  flip = true;
		}
	    }
	}

      // variant positions (in reverse-map)
      positions[ pos ] = nv;

      // allele labels
      alleles.push_back( var_t( a1 , a2 ) );
      
      // A1 frequencies per population
      for (int i=0; i<npop; i++) 
	{
	  double x = 0;
	  if ( ! Helper::str2dbl( tok(i+3) , x ) ) 
	    Helper::halt( "invalid numeric entry in allele freq file" );

	  if ( flip ) x = 1 - x;
	  
	  // ensure we do not have absolute 0 as a MAF
	  if ( x < min_af ) x = min_af;
	  else if ( x > ( 1 - min_af ) ) x = 1 - min_af ; 
	  
	  populations[i].p.push_back( x );
	}

      // filter that may be applied downstream; value should be false
      bad.push_back( bad_site );

      // track number of variants
      ++nv;

    }

  plog << "\nRead frequency data on " << npop << " populations for " << nv << " variants\n";

  // accumulate genotype data ( indiv / variant )
  const int ni = g.indmap.size();
 
  genotypes.resize( ni );
  for (int i=0;i<ni;i++) genotypes[i].resize(nv);
  
  bad.resize( nv , false );

  //
  // Priors -- add option to read in variable priors, but for now
  // 1/NPOP was set above
  //

  F.close();

}


void f_ipop( Variant & v , void * p )
{

  // TODO/check -- do not allow multiallelic markers; assume direct
  // match of reference and alternate sequences

  
  Pseq::IPop * aux = (Pseq::IPop*)p;

  // get genotype per individual, and populate aux->genotypes
  
  // does this variant match one of interest (it should, add we would have added a mask
  // above

  Region reg( v.chromosome() , v.position() , v.stop() );

//   std::map<Region,int>::iterator ii = aux->positions.begin();
//   while ( ii != aux->positions.end() )
//     {
//       std::cout << ii->first << "\t" << ii->second  << "\n";
//       ++ii;
//     }
  
  // must match a position in the list exactly.
  if ( ! aux->has( reg ) ) 
    {      
      return;
    }

  const int j = aux->positions[ reg ];
  
  const int n = v.size();

  // sanity check on VARDB alleles too (as well as SEQDB earlier)
  int ref = 0;
  
  if ( aux->alleles[j].a2 != v.reference() ) 
    {
      aux->bad[j] = true;
      plog.warn( "mismatch in VARDB" , v.displaycore() + " : " + v.reference() + " as REF, expecting " + aux->alleles[j].a2 );
      return;
    }
  else if ( aux->alleles[j].a1 != v.alternate() )
    {
      aux->bad[j] = true;
      plog.warn( "mismatch in VARDB" , v.displaycore() + " : " + v.alternate() + " as ALT, expecting " + aux->alleles[j].a1 );
      return;
    }


  //
  // Get genotypes
  //

  for (int i = 0 ; i < n ; i++ ) 
    {

      if ( v(i).null() ) aux->genotypes[i][j] = -1; // missing
      else 
	{
	  int ac = v(i).allele_count( );	      
	  if ( ac < 0 || ac > 2 ) ac = -1;
	  aux->genotypes[i][j] = ac;	  
	}      
    }
  
}


void Pseq::IPop::calculate()
{

  // Now we have the table of allele freqs, all genotype matrices for
  // all samples.  Calculate posteriors
  
  const int ni = genotypes.size();
  const int nv = bad.size();
  const int np = populations.size();
  
  Out & pout = Out::stream( "ipop" );
  
  pout << std::setprecision(4) ;

  pout << "ID" << "\t"
       << "NVAR" << "\t"
       << "BEST" << "\t"
       << "POP" << "\t"
       << "PROB" << "\n";

  for (int i=0;i<ni;i++) 
    {
      std::vector<double> pp( np );
      const std::vector<int8_t> & geno = genotypes[i];

      int good_snps = 0;

      // each pop
      for (int p = 0 ; p < np ; p++ ) 
	{

	  double x = 0;  // log(P(G|P))
	  
	  // each marker, P(G|P) = \prod P(G_v|P)
	  
	  for (int v = 0 ; v < nv ; v++ )
	    {

	      // just skip missing or bad genotypes

	      if ( bad[v] ) continue;
	      if ( geno[v] == -1 ) continue;
	      
	      const double & frq = populations[p].p[v];
	      
	      // by default, if we did not see the variant, all sites will be reference for the individual
	      // for now, assume this is a reasonable assumption
	      
	      // might also want some diagnostic in the multi-sample
	      // case about the probability of not seeing at least one
	      // non-ref. allele, etc 
	      
	      // use SEQDB to figure out whether A1 or A2 is in fact the reference.

 	      if ( geno[v] == 0 ) 
 		x += log( (1-frq) * (1-frq) );
 	      else if ( geno[v] == 2 ) 
 		x += log( frq * frq );
 	      else if ( geno[v] == 1 ) 
 		x += log( 2 * frq * (1-frq) );
	      
	      if ( p == 0 )  // only need to count this once (@ first population)
		++good_snps; 

	    }
	  
	  // P(G|P) P(P) 
	  
	  x += log( populations[p].prior );
	  	  
	  pp[p] = x;
	  
	} // next population
      
      // Denominator
      double s = sumLogProbs( pp ) ;
      for (int p=0;p<np;p++) { pp[p] -= s; pp[p] = exp( pp[p] ); }
      
      // Get max
      int mx = -1;
      double mp = 0;
      for (int p=0;p<np;p++) 
	if ( pp[p] > mp ) 
	  { mp = pp[p]; mx = p; }
           

      // Now display

      for (int p=0;p<np;p++)
	{
	  pout << g.indmap(i)->id() << "\t"
	       << good_snps << "\t"
	       << ( p == mx ? "B" : "." ) << "\t"
	       << std::fixed << pp[p] << std::scientific << "\t"
	       << populations[p].label << "\n";
	}
    
    } // next individual

}



