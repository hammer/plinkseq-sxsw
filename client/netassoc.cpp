
#include "netassoc.h"
#include "plinkseq.h"
#include "assoc.h"
#include "util.h"

extern GStore g;
extern Pseq::Util::Options args;


// Notes
//   1. uses a DOMINANT coding of rare variants (despite gcounts files having '1' or '2' place-holder)
//   2. we screen for when self-self interactions listed -- in this context, do not want to add to set
//   3. option to exclude any variants present in the seed gene from the test statistic

bool Pseq::NetDB::lookup( const std::string & db , const std::string & gene , const std::string & grp , const int steps )
{

  if ( ! g.locdb.attached() ) Helper::halt( "no attacged LOCDB" ); 

  NetDBase netdb;
  netdb.attach( db );
  if ( ! netdb.attached() ) Helper::halt( "no attached NETDB" );
  netdb.set_locdb( &g.locdb , grp );
  
  // depth=1, threshold=0.0; add-self = T
  std::set<std::string> regs = netdb.connections( gene , 1 , 0 , true );
  
  Out & pout = Out::stream( "partners" );

  std::set<std::string>::iterator ii = regs.begin();
  while ( ii != regs.end() )
    {
      pout << gene << " <---> " << *ii << "\n";
      ++ii;
    }
  return true;

}


bool Pseq::Assoc::net_assoc_test( Mask & m , const Pseq::Util::Options & args )
{

 
  //
  // Have the g-scores been pre-calculated? In this case, read from file
  //
    
  if ( args.has( "file" ) )
    {
      std::map<std::string,Pseq::Assoc::NetDB::Aux_netdet> d = Pseq::Assoc::NetDB::read_scores( args.as_string( "file" ) );
      Pseq::Assoc::NetDB::driver( d , args , m );
      return true;
    }

  
  //
  // Otherwise, iterate over genes, collecting for each gene a list of individuals containing a rare variant 
  // and their 'score' (which could be 0/1, or weighted based on function, frequency, # of variants, etc)
  //
  
  std::map<std::string,Pseq::Assoc::NetDB::Aux_netdet> gscore;
  
  g.vardb.iterate( g_net_assoc_collector , &gscore , m );
  

  //
  // And then, optionally, write these as a self-contained dataset
  //
  
  if ( args.has( "out" ) ) 
    {
      Out output( "netdet" , "intermediate output from net-assoc" );      
      write_gscores( gscore );
    }

  //
  // Call actual network test, need only gscore
  // 

  if ( args.has( "netdb" ) )
    Pseq::Assoc::NetDB::driver( gscore , args , m );
  
  return true;
}



void Pseq::Assoc::NetDB::write_gscores( const std::map<std::string,Pseq::Assoc::NetDB::Aux_netdet> & gscore )
{
  
  Out & FO = Out::stream( "netdet" );

  // header of IDs
  for (int i = 0 ; i < g.indmap.size(); i++ ) 
    FO << g.indmap(i)->id() << "\t";
  FO << "\n";
  
  std::map<std::string,Aux_netdet>::const_iterator i = gscore.begin();
  while ( i != gscore.end() )
    {
      FO << i->first << "\t"
	 << i->second.nvar << "\t"
	 << i->second.ind.size();
      
      std::map<int,double>::const_iterator j = i->second.ind.begin();
      while ( j != i->second.ind.end() )
	{
	  FO << "\t" << j->first 
	     << "\t" << j->second;
	  ++j;
	}
      FO << "\n";
      ++i;
    }  

}


std::map<std::string, Pseq::Assoc::NetDB::Aux_netdet> Pseq::Assoc::NetDB::read_scores( const std::string & filename )
{
  
  Helper::checkFileExists( filename );
  
  InFile F1( filename );

  // format: 
  //  header: individual IDs 
  //          gene-name nvar { id1 sc1 }  { id2 sc2 } ... \n
  //          gene-name nvar { id1 sc1 }  { id2 sc2 } ... \n

  // read individuals
  std::vector<std::string> ids = F1.tokenizeLine();
  std::vector<int> idmap( ids.size() , -1 );
  for (int i=0; i<ids.size(); i++)
    idmap[ i ] = g.indmap.ind_n( ids[i] );
    
  std::map<std::string,Pseq::Assoc::NetDB::Aux_netdet> d;
  
  // read for each gene, counts/scores
  while ( !F1.eof() )
    {     

      std::string gname;
      Aux_netdet a;
      F1 >> gname;
      F1 >> a.nvar;
      if ( gname == "" ) continue;
      int nind;
      F1 >> nind;

      for (int i=0; i<nind; i++)
	{
	  int id;
	  double sc;
	  F1 >> id >> sc;
	  if ( id >= 0 && id < idmap.size() && idmap[id] != -1 )
	    a.ind[ idmap[id] ] = sc > 0 ;	// note -- using a dominant coding here, 0 vs. 1 or 2    
	}
      
      d[ gname ] = a;

    }
  F1.close();

  return d;
}


bool Pseq::Assoc::NetDB::driver( const std::map<std::string,Aux_netdet> & gscore_str , 
				 const Pseq::Util::Options & args , Mask & m )
{
  
  //
  // need a NETDB
  //
  
  std::string netdb_filename;
  if ( args.has( "netdb" ) ) netdb_filename = args.as_string( "netdb" );
  else Helper::halt( "no NETDB specified" );
  NetDBase netdb;
  netdb.attach( netdb_filename );

  

  //
  // Parameters for test
  //

  // get total # of non-missing cases/controls
  int na = 0, nu = 0;
  for (int i=0;i<g.indmap.size();i++)
    {
      if ( ! g.indmap(i)->missing() ) 
	{
	  if ( g.indmap(i)->affected() == CASE ) ++na;
	  else if ( g.indmap(i)->affected() == CONTROL ) ++nu;
	}
    }
  
  if ( na == 0 || nu == 0 ) Helper::halt( "net-assoc requires 1 or more cases and controls" );

  Aux_net_param param( na , // total number of cases (note -- ignores ploidy issues)		       
		       nu , // total number of controls
		       ! args.has( "exclude-seed" )   // do not include SEED in association statistic (default to include)
		       );
  
  //
  // Encode gene-names numerical for faster lookups
  //
  
  std::map<std::string,int> genemap;
  int cnt = 0;
  std::map<std::string,Aux_netdet>::const_iterator i = gscore_str.begin();
  std::map<int,Aux_netdet> gscore;
  while ( i != gscore_str.end() )
    {
      genemap[ i->first ] = cnt;
      gscore[ cnt ] = i->second;
      gscore[ cnt ].name = i->first;
      ++cnt;
      ++i;
    }

  // can free the original storage
  //  gscore_str.clear();

  //
  // Do we want to test all genes?
  //

  std::set<int> grp = m.included_loc();
  std::set<std::string> testset;
  if ( grp.size() ) testset = m.subset_loc( *grp.begin() );
  
  //
  // Find associated network cliques
  //

  long int rseed = time(0);
  int ntests = 1;
  int nrep = -1;
  if ( args.has( "perm" ) ) nrep = args.as_int( "perm" );
  g.perm.initiate( nrep , ntests  );

  std::map<int,Pseq::Assoc::NetDB::Aux_netdet>::const_iterator i1 = gscore.begin();
  while ( i1 != gscore.end() )
    {
      Pseq::Assoc::NetDB::net_test( i1->first , param , rseed, netdb , testset, genemap, gscore );
      ++i1;
    }

  
  //
  // Questions: what if we want to restrict which genes connect in NETDB? (i.e. a Mask for NETDB?) 
  // e.g. 'only brain genes' for example? 
  //

 
}


void g_net_assoc_collector( VariantGroup & vars , void * p )
{

  std::map<std::string,Pseq::Assoc::NetDB::Aux_netdet> * aux = (std::map<std::string,Pseq::Assoc::NetDB::Aux_netdet>*)p;
  
  const std::string name = vars.name();
  
    (*aux)[ name ].nvar = vars.size();
  
  // Current weighting is simply 0 or 1 for at least one variant in gene
  // No weights

  for ( int v = 0 ; v < vars.size(); v++ )
    {      
      int c     = 0; // minor allele
      int c_tot = 0; // total counts	  
      bool altmin = vars(v).n_minor_allele( &c , &c_tot );      
      const int n = vars.n_individuals();
      for (int i = 0; i < n; i++)
	{	  	  
	  if ( ! vars.geno(v,i).null() )
	    {
	      // note, also using a dominent coding here
	      if ( vars.geno(v,i).minor_allele_count( altmin ) )  { (*aux)[ name ].ind[i] = 1; }
	    }
	}
    }
}


void Pseq::Assoc::NetDB::net_test( const int seed , 
				   const Aux_net_param & param ,
				   long int rseed , 				   
				   NetDBase & netdb , 
				   const std::set<std::string> & testset , 
				   const std::map<std::string,int> & genemap , 
				   const std::map<int,Aux_netdet> & gscore )
{
  
  // make a name(string) --> id(int) table to handle genes
  
  // Get all network partners of seed

  std::map<int,Aux_netdet>::const_iterator k = gscore.find( seed );
  std::string gname = k->second.name;

  // do we want to test this gene? 
  if ( testset.size() > 0 ) 
    {
      if ( testset.find( gname ) == testset.end()) return; 
    }


  std::set<int> connections = netdb.connections( gname , genemap );
  
  // Set up permutations
  g.perm.seed( rseed );
  g.perm.reset();

  double ca = 0 , cu = 0;
  int nvar = 0;
  std::set<Aux_connection> endset;
  double original = net_statistic( seed , param , &netdb, connections, &endset, gscore , genemap, &nvar, &ca , &cu );
  g.perm.score( original );

  // permute
  
  const int R = g.perm.replicates();  

  for (int p = 1 ; p < R ; p++ )
    {
      double stat = net_statistic( seed , param , &netdb, connections, NULL , gscore , genemap ); 
      if ( ! g.perm.score( stat ) ) break; 
    }

  Out & pout = Out::stream( "net.assoc" );

  // report : first seed 

  pout << "SEED" << "\t"
       << gname << "\t" 
       << connections.size() << "\t" 
       << endset.size() << "\t"
       << nvar << "\t" ;
  
  if ( ca + cu == 0 )  // can happen when --exclude-seed but no PPI partners
    pout  << 0 << ":" << 0 << "\t" 
	  << 1 << "\t";
  else
    pout  << ca << ":" << cu << "\t" 
	  << g.perm.pvalue(0) << "\t";
  
  pout << g.locdb.alias( gname , false ) << "\n";
  
  // partners

  std::set<Aux_connection>::iterator i = endset.begin();
  while ( i != endset.end() ) 
    {
      pout << "PARTNER" << "\t"
	   << gname << "\t";
      
      std::string n = gscore.find( i->extension )->second.name ; 
      
      if ( n != gname && i->parent != seed && i->parent != -1 ) 
	pout << gscore.find( i->parent )->second.name << " --> ";

      pout << n << "\t"
	   << i->acnt << "\t"
	   << i->ucnt << "\t";
      pout << g.locdb.alias( gscore.find( i->extension )->second.name , false ) << "\n";
      ++i;
    }
}



double Pseq::Assoc::NetDB::net_statistic( const int seed, 
					  const Aux_net_param & param , 
					  NetDBase * netdb , 
					  std::set<int> & connections , 
					  std::set<Aux_connection> * endset , 
					  const std::map<int,Aux_netdet > & gscore , 
					  const std::map<std::string,int> & genemap , 
					  int * pnvar , double * pa , double * pu )
{
  
  // calculate case/control enrichment   

  // optimize over all possible gene parnters, using greedy algorithm.
  
  // start point is always the 1st degree connections that are listed in endset
    
  std::set<next_node_t> next_node;

  std::set<int> inset;

  //
  // Original score is single gene test (inset is empty)
  //
  
  double stat0 =  stat_adder( seed, param , inset , gscore );
    
  //
  // initially, populate with 1st degree connections
  //
  
  std::set<int>::iterator i = connections.begin();
  while( i != connections.end() )
    {
      next_node.insert( next_node_t( *i ) ); 
      ++i;
    }
  
  //
  // Optimisation loop 
  //

  while ( 1 ) 
    {

      bool found_better = false;
      
      // test each potential extension and select the best

      // retain if it is better than original test (by a margin?)

      next_node_t * best = NULL;
      
      std::set<next_node_t>::iterator i = next_node.begin();
      
      int cnt = 0;

      while ( i != next_node.end() )
	{
	  
	  // skip genes we've already seen and made the score worse
	  if ( ! i->retain )
	    {
	      ++i;
	      ++cnt;
	      continue;
	    }

	  std::set<int> testset = inset;
	  testset.insert( i->extension );
	  
	  // New statistic
	  double stat1 = stat_adder( seed , param , testset , gscore );
	  
	  if ( stat1 > stat0 ) 
	    {
	      best = (next_node_t*)&(*i);
	      stat0 = stat1;	      
	    }
	  else 
	    const_cast<bool&>(i->retain) = false;

	  ++i;
	  ++cnt;
	}
      
      // we didn't do any better
      if ( ! best ) 
	{
	  break;
	}

      // otherwise, add that gene to the list, and all of it's neighbours (up to a point)
      
      inset.insert( best->extension );
      
      if ( best->depth < 2 )
	{
	  std::map<int,Aux_netdet>::const_iterator k = gscore.find( best->extension );
	  std::string gname = k->second.name;

	  std::set<int> ext = netdb->connections( gname , genemap );	  
	  std::set<int>::iterator ii = ext.begin();
	  while ( ii != ext.end() )
	    {	      	      
	      if ( *ii != seed ) // ensure do not add seed to the list 
		next_node.insert( next_node_t( *ii , best->extension , best->depth + 1 ) );
	      ++ii;
	    }
	}
      
      
      // now go back to try to add another gene
    }


  //  
  // for original, copy over the best set, with some other info
  //

  if ( endset ) 
    {

      std::set<next_node_t>::iterator i = next_node.begin();

      while ( i != next_node.end() )
	{

	  // only consider genes in set 

	  if ( inset.find( i->extension ) == inset.end() ) 	    
	    { ++i; continue; }
	    	      
	  Aux_connection c;
	  c.extension = i->extension;
	  c.parent = i->parent;
	  	  
	  // get gene specific A/U counts
	  int d_pnvar;
	  double d_a;
	  double d_u; 
	  std::set<int> dummy;
	  double dummy_stat = stat_adder( - i->extension , param , dummy , gscore , &d_pnvar , &d_a , &d_u );
	  
	  c.acnt = (int)d_a;
	  c.ucnt = (int)d_u;
	  endset->insert( c );

	  ++i;
	}
      

      // and insert seed node in endset too
      
      if ( param.add_seed ) 
	{
	  int d_pnvar;
	  Aux_connection c;
	  double d_a;
	  double d_u; 
	  std::set<int> dummy;
	  double dummy_stat = stat_adder( -seed , param , dummy , gscore , &d_pnvar , &d_a , &d_u );
	  c.extension = seed;
	  c.parent = -1;
	  c.acnt = (int)d_a;
	  c.ucnt = (int)d_u;
	  endset->insert( c );
	}
 
    }

  //
  // run for a last time with best model, to populate some stats
  //

  return stat_adder( seed, param , inset , gscore , pnvar, pa , pu );
  
}


double Pseq::Assoc::NetDB::stat_adder( const int seed , 
				       const Aux_net_param & param , 
				       std::set<int> & inset , 				       
				       const std::map<int, Aux_netdet > & s , 
				       int * pnvar , double * pa , double * pu )
  
{
  
  double ca = 0 , cu = 0;

  // add 'seed' (optionally)

  // note, if seed < 0 this means we want to force-add the seed

  // this is because we use this function to calculate C/C counts for the
  // 'endset' also, and so this is a cludge to ensure they are added
  // in this case, 'inset' above will always be empty too.
  
  bool eval_seed_only_kludge  = seed < 0 ;
  
  if ( param.add_seed || eval_seed_only_kludge ) 
    {
      
      std::map<int, Aux_netdet>::const_iterator i = s.find( eval_seed_only_kludge ? -seed : seed );
      if ( i == s.end() ) return 0;
      
      if ( pnvar ) *pnvar += i->second.nvar;
      std::map<int,double>::const_iterator ii = i->second.ind.begin();  
      while ( ii != i->second.ind.end() )
	{
	  affType aff = g.indmap( g.perm.pos( ii->first ) )->affected();
	  if ( aff == CASE ) ca += ii->second;
	  else if ( aff == CONTROL ) cu += ii->second;
	  ++ii;
	}  
    }

  
  // all partners in 'inset'
  std::set<int>::iterator k = inset.begin();
  while ( k != inset.end() )
    {
      std::map<int, Aux_netdet>::const_iterator ii = s.find( *k );
      if ( ii == s.end() ) 
	{
	  ++k;
	  continue;
	}
      
      if ( pnvar ) *pnvar += ii->second.nvar;
      std::map<int,double>::const_iterator iii = ii->second.ind.begin();
      while ( iii != ii->second.ind.end() )
	{
	  affType aff = g.indmap( g.perm.pos( iii->first ) )->affected();
	  if ( aff == CASE ) ca += iii->second;
	  else if ( aff == CONTROL ) cu += iii->second;
	  ++iii;
	}
      
      ++k;
    }
  
  if ( pa ) *pa = ca;
  if ( pu ) *pu = cu;
  

  // based on dominant coding, determine chi-sq
  
  const int da = param.tot_case - ca;
  const int du = param.tot_control - cu;
  
  if ( da < 0 || du < 0 ) 
    Helper::halt( "internal error in net-assoc: da < 0 || du < 0" ); 
  
  bool zero_count = ca == 0 || cu == 0 || da == 0 || du == 0; 
  
  double odds = 0;
  
  if ( zero_count ) 
    odds = ( ( ca + 0.5 ) * ( du + 0.5 ) ) / ( ( cu + 0.5 ) * ( da + 0.5 ) ) ;   
  else
    odds = (double)( ca * du ) / (double)( cu * da );

  // 1-sided test -- looking only for increase in risk
  if ( odds <= 1 ) return 0;
  
  // return chi-sq are test statistic
  return Helper::chi2x2( ca , cu , da , du ) ; 
	  
  // Original statistic
  //return ( ca / (ca+cu) ) * ( ca-cu ) ;
  
}


