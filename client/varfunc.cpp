
#include "func.h"
#include "pseq.h"
#include "assoc.h"
#include "genic.h"

extern GStore g;
extern Pseq::Util::Options options;

struct AuxConcordanceStats 
{  
  std::map<int,int>  once;
  std::map<int2,int> conc;
  std::map<int2,int> disc;
  std::map<int2,int> disc_ref_nonref;
};


struct AuxConcordanceCheck
{

  AuxConcordanceCheck() 
  {
    report_all = false;
  }

  bool report_all;

  // per-individual statistics
  std::map< std::string , AuxConcordanceStats > istats;

  // total dataset
  AuxConcordanceStats tstats;
  
  
  void evaluate( const Variant & v )
  {
    

    //
    // Only need to consider variants with >1 sample
    //

    if ( ! v.multi_sample() ) return;


    //
    // Track output for this one variant
    //

    // # of indiv with 1+ non-missing copy of this genotype
    int n_multiple_obs = 0;

    // # of indiv with only 1 non-missing genotype | 2+ observed
    int n_concordant_obs = 0;
    
    // list of discordant individuals
    std::set<std::string> disc_indiv;

    // track ref/alt, alt/ref, alt/alt status for discordant calls

    int cnt_ra = 0 , cnt_ar = 0 , cnt_aa = 0;
    
    //
    // Look at each individual
    //

    for (int i=0; i<v.size(); i++)
      {
	
	// Get all genotypes
	std::map<int, const Genotype *> gm = v.all_genotype(i);

	// Only consider if >1 genotype seen
	if ( gm.size() < 2 ) continue;
	
	// Consider all unique pairs of files
	std::map<int, const Genotype *>::iterator j1 = gm.begin();
	while ( j1 != gm.end() )
	  {
	    std::map<int, const Genotype *>::iterator j2 = j1;
	    ++j2;
	    while ( j2 != gm.end() )
	      {
		const Genotype * g1 = j1->second;
		const Genotype * g2 = j2->second;
		
		// both genotypes observed?
		// this should handle multi-allelic markers? (check)

		if ( g1->notnull() && g2->notnull() )
		  record( v.ind(i)->id() , 
			  j1->first , j2->first, 
			  g1, g2, 
			  cnt_ra , cnt_ar , cnt_aa , 
			  v );
		
		++j2;
	      }
	    ++j1;
	  }
	

	//
	// track seen_once stats; also accumulate per-variant stats
	//
	
	int nobs = 0;
	std::set<std::string> uniq_obs;

       	j1 = gm.begin();
	while ( j1 != gm.end() )
	  {

	    if ( j1->second->notnull() ) 
	      {
		record( v.ind(i)->id() , j1->first );
		++nobs;
		
		// plog << v.ind(i)->id() << "\t" 
		// << j1->first << "\t" 
		// << v.psample( j1->first )->label( *(j1->second) , true )  << "\n";
		
		uniq_obs.insert( v.psample( j1->first )->label( *(j1->second) , true ) );
	      }

	    ++j1;
	  }
	
	//
	// Accumulate per-variant statistics
	//

	if ( nobs > 1 ) 
	  {
	    ++n_multiple_obs;
	    if ( uniq_obs.size() == 1 ) 
	      ++n_concordant_obs;
	    else
	      {
		disc_indiv.insert( v.ind(i)->id() );
		
	      }
	  }
	
      }
    
    
    //
    // Report per-variant stats
    //
    
    if ( n_multiple_obs > 0 || report_all )
      {
	plog << "_VARIANT" << "\t"
		  << v << "\t" 
		  << v.n_samples() << "\t"
		  << n_multiple_obs << "\t"
		  << n_concordant_obs << "\t";

	if ( n_multiple_obs > 1 )
	  plog << (double)(n_concordant_obs)/(double)n_multiple_obs << "\t";
	else
	  plog << "NA\t";

	if ( v.n_samples() == 2 ) 
	  plog << cnt_ra << "\t" 
		    << cnt_ar << "\t"
		    << cnt_aa << "\t";
 	else
	  plog << "NA" << "\t" 
		    << "NA" << "\t"
		    << "NA" << "\t";

	if ( disc_indiv.size() == 0 ) 
	  plog << ".";
	else
	  {
	    std::set<std::string>::iterator i = disc_indiv.begin();
	    while ( i != disc_indiv.end() )
	      {
		if ( i != disc_indiv.begin() ) plog << ",";
		plog << *i;
		++i;
	      }
	  }
	plog << "\n";
      }
  }


  void report_helper( AuxConcordanceStats & stats , const std::string & label )
  {
    
    std::map<int,int>::iterator i = stats.once.begin();
    while ( i != stats.once.end() )
      {
	std::map<int,int>::iterator j = i;
	while ( j != stats.once.end() )
	  {
	    
	    if ( i == j ) { ++j; continue; }
	    
	    plog << label << "\t"
		      << g.vardb.file_tag( i->first ) << "\t"
		      << g.vardb.file_tag( j->first ) << "\t"
		      << i->second << "\t"
		      << j->second << "\t";
	    
	    int conc = i->first < j->first ? 
	      stats.conc[int2(  i->first ,  j->first ) ] :
	      stats.conc[int2(  j->first ,  i->first ) ] ;
	    
	    int disc = i->first < j->first ? 
	      stats.disc[int2(  i->first ,  j->first ) ] :
	      stats.disc[int2(  j->first ,  i->first ) ] ;
	    
	    int disc_ref_nonref = stats.disc_ref_nonref[ int2(  i->first ,  j->first ) ] ;
	    int disc_nonref_ref = stats.disc_ref_nonref[ int2(  j->first ,  i->first ) ] ;
	    
	    plog << disc << "\t" 
		      << conc << "\t"
		      << ( conc+disc>0 ? Helper::flt2str((double)conc/(double)(conc+disc)) : "NA" ) << "\t"
		      << disc_ref_nonref << "\t" 
		      << disc_nonref_ref << "\t"
		      << ( disc - disc_ref_nonref - disc_nonref_ref ) << "\n";

	    ++j;
	  }

	++i;
    }

  }


  //
  // Final report after all variants considered
  //

  void report()
  {

    //
    // For each individual
    //

    std::map< std::string , AuxConcordanceStats >::iterator ii = istats.begin();
    while ( ii != istats.end() )
      {
	report_helper( ii->second , "_INDIV\t" + ii->first );
	++ii;
      }
	
    	
    //
    // for all files seen, give total variant count
    //

    report_helper( tstats , "_TOTAL" );    
    
  }


  //
  // Record pairs concordance statistics
  //

  void record( const std::string & id , 
	       const int f1 , const int f2 , 
	       const Genotype * g1 , const Genotype * g2 , 
	       int & ra, int & ar , int & aa ,  
	       const Variant & v )
  {
    
    bool conc = v.concordant( f1, g1, f2, g2 );
    
    // genotypes will be non-missing here, so:
    bool ref1 = conc ? false : g1->reference();
    bool alt1 = ! ref1;

    bool ref2 = conc ? false : g2->reference();
    bool alt2 = ! ref2;
    
    if ( !conc )
      {
	if ( ref1 && alt2 ) ++ra;
	else if ( alt1 && ref2 ) ++ar;
	else ++aa;
      }
    
    if ( conc ) 
      tstats.conc[ int2(f1,f2) ]++;    
    else
      {
	tstats.disc[ int2(f1,f2) ]++;    
	if ( ref1 && alt2 ) tstats.disc_ref_nonref[ int2(f1,f2) ]++;
	else if ( alt1 && ref2 ) tstats.disc_ref_nonref[ int2(f2,f1) ]++;	
      }
    
    if ( conc ) 
      istats[id].conc[ int2(f1,f2) ]++;
    else
      {
	istats[id].disc[ int2(f1,f2) ]++;
	if ( ref1 && alt2 ) istats[id].disc_ref_nonref[ int2(f1,f2) ]++;
	else if ( alt1 && ref2 ) istats[id].disc_ref_nonref[ int2(f2,f1) ]++;
      }

    if ( ! conc )
      plog << "_GENO" << "\t"
		<< id << "\t"
		<< f1 << "\t"
		<< f2 << "\t"
		<< v << "\t"
		<< v.psample(f1)->label( *g1 ) << "\t"
		<< v.psample(f2)->label( *g2 ) << "\t"
                << g1->nonreference() << "\t"
		<< g2->nonreference() << "\n";

  }

  //
  // Record that non-missing call seen in file
  //

  void record( const std::string & id , int f1 )
  {
    tstats.once[f1]++;
    istats[id].once[f1]++;
  }
  

};


void f_conc( Variant & v , void * p )
{
  ((AuxConcordanceCheck*)p)->evaluate(v);
}

bool Pseq::VarDB::check_concordance( Mask & m )
{
  AuxConcordanceCheck aux;

  if ( options.key("report-all") )
    aux.report_all = true;

  plog << "_VARIANT" << "\t"
	    << "VAR" << "\t" 
	    << "N_SAMPLES" << "\t"
	    << "N_MULT_OBS" << "\t"
	    << "N_CONC_OBS" << "\t"
	    << "CONC_RATE" << "\t"
	    << "DISC_INDIV" << "\n";

  g.vardb.iterate( f_conc , &aux , m );
  aux.report();  
}



bool Pseq::VarDB::lookup_list( const std::string & filename , Mask & mask )
{

  Helper::checkFileExists( filename );

  // From LOCDB, take set of gene-groups
  // From REFDB, take set of ref-variants
  
  std::set<std::string> locs = options.get_set( "loc" );
  std::set<std::string> refs = options.get_set( "ref" );
  
  bool append_phe = g.vardb.attached() && g.inddb.attached() && g.phmap.type() == PHE_DICHOT;    
  bool append_loc = g.locdb.attached() && locs.size() > 0;
  bool append_ref = g.refdb.attached() && refs.size() > 0;
  bool append_seq = g.seqdb.attached();

  bool vardb = g.vardb.attached();
  bool append_annot = g.seqdb.attached() && ( options.key("annot") || options.key( "annotate" ) );

  if ( append_annot ) 
    {
      std::string annot_transcripts = PLINKSeq::DEFAULT_LOC_GROUP() ;
      if ( options.key( "annot" ) ) annot_transcripts = options.as<std::string>( "annot" );
      else if ( options.key( "annotate" ) ) annot_transcripts = options.as<std::string>( "annotate" );      
      Annotate::load_transcripts( LOCDB, annot_transcripts );
    }
  
  InFile F1( filename );
  
  if ( vardb )
    {
      if ( append_phe ) plog << "##casecon,1,String,\"Case/control minor allele counts\"\n";
      else plog << "##count,1,String,\"Minor allele counts\"\n";
    }
  
  if ( append_loc )
    {
      std::set<std::string>::iterator i = locs.begin();
      while ( i != locs.end() )
	{
	  plog << "##" << *i << ",String,\"LOCDB group\"\n";
	  ++i;
	}
    }
  
  if ( append_ref ) 
    {
      std::set<std::string>::iterator i = refs.begin();
      while ( i != refs.end() )
	{
	  plog << "##" << *i << ",String,\"REFDB group\"\n";
	  ++i;
	}
    }
  
  if ( append_annot )
    {
      plog << "##func,1,String,\"Genomic annotation\"\n";
      plog << "##transcript,1,String,\"Transcript ID\"\n";
      plog << "##genomic,1,String,\"Genomic DNA change\"\n";
      plog << "##codon,1,String,\"Codon change\"\n";
      plog << "##protein,1,String,\"Any nonsynon amino acid change\"\n";
      plog << "##worst,1,String,\"Worst annotation\"\n";
      plog << "##class,1,String,\"Summary of all annotations\"\n";
    }

  if ( append_seq )
    {
      plog << "##seqdb_ref,1,String,\"SEQDB reference sequence\"\n";
    }

  if ( append_ref ) 
    {

      std::set<std::string>::iterator i = refs.begin();
      while ( i != refs.end() )
	{
	  plog << "##ref_" << *i << ",.,String,\"REFDB annotation for group " << *i << "\"\n";
	  ++i;
	}
    }

  if ( ! ( vardb || append_loc || append_ref || append_seq || append_annot ) ) 
    Helper::halt("no information to append");
  
  // Assume format is either a) a single site (1 col) OR b) single, REF and ALT (3 col)
  while ( !F1.eof() )
    {
      
      bool okay = true;
      
      std::vector<std::string> line = F1.tokenizeLine( " \t" );
      
      if ( ! ( line.size() == 1 || line.size() == 3 ) ) continue;
      
      const std::string s = line[0];

      Region region( s , okay);
      
      std::string ref_allele = line.size() == 3 ? line[1] : ".";
      std::string alt_allele = line.size() == 3 ? line[2] : ".";
      
      if ( append_annot && okay ) 
	{
	  
	  Variant var;
	  var.chromosome( region.start.chromosome() );
	  var.position( region.start.position() );
	  var.stop( var.position() + ref_allele.size() - 1 );
	  var.consensus.reference( ref_allele );
	  var.consensus.alternate( alt_allele );
	  
	  bool exonic = Annotate::annotate( var );
	  
	  std::string annot = var.meta.get1_string( PLINKSeq::ANNOT() );
	  
	  // detailed annotation vector, primary annotation
	  
	  plog << var.coordinate() << "\t"
	       << "func" << "\t"
	       << var.meta.as_string( PLINKSeq::ANNOT_TYPE() , "," ) << "\n";
	  
	  plog << var.coordinate() << "\t"
	       << "transcript" << "\t"
	       << var.meta.as_string( PLINKSeq::ANNOT_GENE() , "," ) << "\n";

	  plog << var.coordinate() << "\t"
	       << "genomic" << "\t"
	       << var.meta.as_string( PLINKSeq::ANNOT_CHANGE() , "," ) << "\n";

	  plog << var.coordinate() << "\t"
	       << "codon" << "\t"
	       << var.meta.as_string( PLINKSeq::ANNOT_CODON() , "," ) << "\n";
				      
	  plog << var.coordinate() << "\t"
	       << "protein" << "\t"
	       << var.meta.as_string( PLINKSeq::ANNOT_PROTEIN() , "," ) << "\n";

	  // worst-case consensus annotation
	  
	  plog << var.coordinate() << "\t"
	       << "worst" << "\t"
	       << annot << "\n";
          
	  // summary annotation
	  
	  plog << var.coordinate() << "\t"
	       << "class" << "\t"
	       << var.meta.get1_string( PLINKSeq::ANNOT_SUMMARY()) << "\n";
	  
	}
      

      if ( ! okay ) 
	{
	  plog.warn( "not a valid region" , s );
	  continue;
	}
      

      // Fetch from SEQDB

      if ( append_seq ) 
	{
	  plog << s << "\t" 
	       << "seqdb_ref" << "\t"
	       << g.seqdb.lookup( region ) << "\n";	  
	}

      
      // Fetch from VARDB 
      
      std::set<Variant> vars = g.vardb.fetch( region );
      
      if ( vars.size() == 0  )
	{
	  if ( vardb ) 
	    plog << "_VAR" << "\t"
		 << s << "\t"
		 << 0 << "\t"
		 << "_NO_VARIANTS" << "\t"
		 << ( append_phe ? "NA\tNA\t" : "NA\t" )
		 << "\n";
	}
      
      int cnt = 0;
      std::set<Variant>::iterator v = vars.begin();
      while ( v != vars.end() )
	{
	  
	  plog << "_VAR" << "\t"
	       << s << "\t"
	       << ++cnt << "\t"
	       << *v << "\t";		    


	  // 
	  // Either allele counts; or stratify by case/control
	  //

	  int case_n = 0 , control_n = 0;
	  
	  for (int j=0; j < v->size() ; j++)
	    {
	      
	      if ( (*v)(j).nonreference() )
		{
		  if ( append_phe ) 
		    {
		      if ( v->ind(j)->affected() == CASE )
			++case_n;
		      else
			++control_n;
		    }
		  else
		    ++control_n;
		}
	    }
	   
   
	  if ( append_phe )
	    plog << case_n << "\t" ;
	  plog << control_n << "\t";
	  
	  ++v;
	  
	  plog << "\n";
	}

      	  
      //
      // Locus DB ? 
      //
      
      if ( append_loc ) 
	{
	  std::set<std::string>::iterator i = locs.begin();
	  while ( i != locs.end() )
	    {
	      std::set<Region> rregs = g.locdb.get_regions( *i , region );
	      if ( rregs.size() == 0 ) 
		{
		  plog << "_REFLOC" << "\t"
		       << s << "\t"
		       << *i << "\t"
		       << "_NO_LOCI" << "\n";
		}
	      
	      std::set<Region>::iterator j = rregs.begin();
	      while ( j != rregs.end() )
		{
		  plog << "_REFLOC" << "\t"
		       << s << "\t"
		       << *i << "\t"
		       << j->coordinate() << ":"
		       << j->name << "\n";
		  ++j;
		}
	      ++i;
	    }
	}
      
      
      //
      // Reference variants? 
      //
      
      if ( append_ref ) 
	{
	  std::set<std::string>::iterator i = refs.begin();
	  while ( i != refs.end() )
	    {

	      std::set<RefVariant> rvars = g.refdb.lookup( region , *i );
	      
	      if ( rvars.size() == 0 ) 
		{
		  plog << s << "\t"
		       << "ref_" << *i << "\t"
		       << "." << "\n";
		}
	      
	      std::set<RefVariant>::iterator j = rvars.begin();
	      while ( j != rvars.end() )
		{
		  plog << s << "\t"
		       << "ref_" << *i << "\t"
		       << *j << "\n";
		  ++j;
		}
	      ++i;
	    }
	}


    } // next input region
  
  F1.close();
  return true;

}

struct Aux_annotate_loc {
  Aux_annotate_loc( LocDBase * locdb , const std::string & grp , int nm, bool s ) 
    : locdb(locdb) , grp(grp) , nm(nm), show_sub(s) { }
  LocDBase * locdb;
  std::string grp;
  bool show_sub;
  int nm; 
};

void f_annotate_loc( Variant & v , void * p )
{
  Aux_annotate_loc * aux = (Aux_annotate_loc*)p;

  std::set<Region> loc = aux->locdb->get_regions( aux->grp , v );

  if ( loc.size() == 0 ) 
    {
      plog << v << "\t"
	   << 0 << "\t"
	   << "." << "\t"
	   << "." << "\t"
	   << ( aux->show_sub ? "." : "" ) << "\t";
      for (int j=0;j<aux->nm;j++) plog << (j!=0 ? "\t" : "" ) << ".";
      plog << "\n";
    }
  else
    {
      int cnt = 1;
      std::set<Region>::iterator i = loc.begin();
      while ( i != loc.end() ) 
	{
	  plog << v << "\t"
	       << cnt++ << "\t"
	       << i->name << "\t"
	       << aux->locdb->alias( i->name , false ) << "\t"
	       << ( aux->show_sub ? "0\t" : "" ) 
	       << i->meta.display_row(".") << "\n";	  
	  ++i;
	}
    }
}

bool Pseq::VarDB::annotate_loc( const std::string & grp, Mask & m )
{  
  
  plog << "VAR" << "\t"
       << "N" << "\t"
       << "LOC" << "\t"
       << "ALIAS" << "\t"
       << MetaInformation<LocMeta>::display_header() << "\n";
  

  int nm = MetaInformation<LocMeta>::n_keys();
  
  // display headers for LOCDB group 'grp'
  Aux_annotate_loc aux( &g.locdb , 
			grp , 
			nm , 
			options.key("show-subregions"));
  g.vardb.iterate( f_annotate_loc , &aux , m );

  return true;
}



bool Pseq::VarDB::simple_sim()
{

  // a very basic wrapper purely to generate some data to be fed to 
  // to various genic tests, primarily in a sanity-check/debug mode
  // although could also be used for power calculations, etc

  Pseq::Assoc::AuxGenic a;
  a.g     = &g;
  a.rseed = time(0);

  plog.data_reset();  
  plog.data_group_header( "LOCUS" );
  plog.data_header( "POS" );
  plog.data_header( "ALIAS" );
  plog.data_header( "NVAR" );
  plog.data_header( "TEST" );
  plog.data_header( "P" );
  plog.data_header( "I" );
  plog.data_header( "DESC" );  
  plog.data_header_done();
  
  a.vanilla = a.burden = a.uniq = true;
  a.mhit = options.key("mhit" );
  a.vt = options.key("vt");
  a.fw = options.key("fw");
  a.calpha = options.key("calpha");
  a.cancor = options.key("cancor");

  const int ntests = a.n_tests();
  const int nrep = -1;
  
  // Format:

  // # cases
  // # controls
  // # genes  { # of times to repeat the blockset below }

  // { a block } 
  // 10   # num of instances
  // 2    # num of sites in haplotype (i.e. 1 == SNP)
  // 2    # num of alternate alleles
  // 01   # first alt haplotype  
  // 11   # second
  // 2    # cases, first alternate
  // 0    #        second
  // 1    # controls, first alternate
  // 1    #           second
  
  // { next block }

  // 0 # implies end of block
  


  int ncase, ncontrol, ngene;
  std::cin >> ncase >> ncontrol >> ngene;

  Mask dummy;
  
  //
  // construct our own phenotype map
  //

  std::vector<std::string> id_list;
  for (int i=0;i<ncase;i++) id_list.push_back( "A" + Helper::int2str(i) );
  for (int i=0;i<ncontrol;i++) id_list.push_back( "U" + Helper::int2str(i) );
  int n = g.indmap.populate( id_list );
  for (int i=0; i< n; i++ )
    {
      Individual * person = g.indmap.ind(i);
      if ( i < ncase ) person->affected( CASE );
      else person->affected( CONTROL );
    }

  
  // Atfer creating the IndMap, we can initiate the permutation class
  // (which requires a constructed IndMap, to specify the # of individuals)

  g.perm.initiate(  nrep , ntests );
  a.fix_null_genotypes = options.key("fix-null");


  // read block-level data
  
  std::vector<int> blk_nvar;
  std::vector<int> blk_nsite;
  std::vector<int> blk_alt;
  std::vector< std::vector<std::string> > blk_allele;
  std::vector< std::vector<int> > blk_case;
  std::vector< std::vector<int> > blk_control;

  while ( ! std::cin.eof()  )     
    {
      int nvar, nsite, alt;
      std::cin >> nvar;
      if ( std::cin.eof() || nvar == 0 ) break;
      std::cin >> nsite >> alt;

      blk_nvar.push_back( nvar );      
      blk_nsite.push_back( nsite );
      blk_alt.push_back( alt );
      std::vector<std::string> tmp;
      std::vector<int> tmpi;
      blk_allele.push_back( tmp );
      blk_case.push_back( tmpi );
      blk_control.push_back( tmpi );

      for (int i=0;i<alt;i++)
	{
	  std::string allele;
	  int na, nu;
	  std::cin >> allele >> na >> nu;
	  blk_allele[ blk_allele.size()-1 ].push_back( allele );
	  blk_case[ blk_case.size()-1 ].push_back( na );
	  blk_control[ blk_control.size()-1 ].push_back( nu );
	}
      
    }



  // build gene and run test, ngene times
  // assumes biallelic allele coding (0,1) in specification text
  
  int seed = CRandom::rand( 10000000 );
  
  for (int i=0; i<ngene; i++)
    {
      // if we do not do this, the genic tests reset the RNG after each 
      // gene: TODO< allow access to a separate, "safe" RNG

      CRandom::srand( seed );

      VariantGroup vars( dummy );
      vars.name( "GENE" + Helper::int2str( i ) );

      int totv = 0;

      for (int b=0; b<blk_nvar.size(); b++) // block
	{

	  for (int v=0; v<blk_nvar[b]; v++) // variant
	    {
	      
	      // select random cases/controls chromosomes
	      // -1 means ref; alt-alleles 0,..,na-1
	      std::vector<int> acasepat( ncase ,-1);
	      std::vector<int> acasemat( ncase ,-1);
	      std::vector<int> acontrolpat( ncontrol ,-1);
	      std::vector<int> acontrolmat( ncontrol ,-1);
	      
	      int na = blk_allele[b].size();
	      for (int a=0; a < na ; a++)
		{

		  // cases
		  for (int p=0; p<ncase; p++) 
		    {
		      bool a1 = CRandom::rand() <= (blk_case[b][a] / (double)(2*ncase));
		      bool a2 = CRandom::rand() <= (blk_case[b][a] / (double)(2*ncase));
		      if ( a1 && acasepat[p] == -1 ) acasepat[p] = a;
		      if ( a2 && acasemat[p] == -1 ) acasemat[p] = a;
		    }
		  
		  // controls
		  for (int p=0; p<ncontrol; p++) 
		    {
		      bool a1 = CRandom::rand() <= (blk_control[b][a]/ (double)(2*ncontrol));
		      bool a2 = CRandom::rand() <= (blk_control[b][a] / (double)(2*ncontrol));
		      if ( a1 && acontrolpat[p] == -1 ) acontrolpat[p] = a;
		      if ( a2 && acontrolmat[p] == -1 ) acontrolmat[p] = a;
		    }		
		}
	    

	      // now we've assigned all individuals to haplotypes; figure out the
	      // implied per-variant allele coding and store, for 'ns' SNPs in the gene
	      
	      int ns = blk_nsite[b];
	      for (int w=0; w<ns; w++) // sites at variant
		{
		  
		  Variant var;
		  var.consensus.reference("A");
		  var.consensus.alternate("B");
		  var.chromosome(1);
		  var.position( ++totv );
		  var.attach( &g.indmap );
		  var.consensus.calls.size( ncase + ncontrol );
		  int cnt = 0;
		  for ( int i=0; i < ncase; i++ )
		    {
		      int geno = 0;
		      if ( acasepat[i] >= 0 && blk_allele[b][acasepat[i]].substr(w,1) == "1" ) ++geno;
		      if ( acasemat[i] >= 0 && blk_allele[b][acasemat[i]].substr(w,1) == "1" ) ++geno;		      
		     		      
		      var(i).set_alternate_allele_count( geno );
		      ++cnt;
		    }
		  
		  for ( int i=0; i < ncontrol; i++ )
		    {
		      int geno = 0;
		      if ( acontrolpat[i] >= 0 && blk_allele[b][acontrolpat[i]].substr(w,1) == "1" ) ++geno;
		      if ( acontrolmat[i] >= 0 && blk_allele[b][acontrolmat[i]].substr(w,1) == "1" ) ++geno;
		      var( cnt ).set_alternate_allele_count( geno );
		      ++cnt;
		    }
		  
		  vars.force_add( var );

		}

	    } // next variant
	} // next block
      

      // view
      //      std::cout << vars.dump() << "\n";
      
      seed = CRandom::rand( 1000000 );

      // apply genic-association tests
      g_set_association( vars , &a );

      
      
    }
        

  return true; 

}



//
// Proximity Scan
//

struct Aux_proximity{ 

  int dist;
  
  Variant lastvar;
  
  std::set<Variant> store;
  
};


void f_proxscan( Variant & var , void * p )
{

  Aux_proximity * aux = (Aux_proximity*)p;
  
  // Is this close enough to the last variant seen?

  if ( var.chromosome() == aux->lastvar.chromosome() && 
       var.position() - aux->lastvar.stop() <= aux->dist ) 
    {

      aux->store.insert( aux->lastvar );

      const int n = var.size();
      
      std::set<Variant>::iterator i = aux->store.begin();
      while ( i != aux->store.end() )
	{

	  plog.data_group( i->displaycore() + "," + var.displaycore() );
	  
	  // Counts

	  int n1, m1, n2, m2;
	  bool altmin1 = var.n_minor_allele( n1, m1 );
	  bool altmin2 = aux->lastvar.n_minor_allele( n2, m2 );
	  
	  // DIST
	  plog.data( var.position() - aux->lastvar.stop() );

	  // Calculate correlation coefficient between pair of SNPs
	  
	  Data::Vector<double> alt1 = VarFunc::alternate_allele_count( var );	  
	  Data::Vector<double> alt2 = VarFunc::alternate_allele_count( *i  );

	  Data::Matrix<double> d( n , 0 );
	  d.add_col( alt1 );
	  d.add_col( alt2 );
	  
	  int mis = 0;
	  for (int ik=0; ik<d.dim1(); ik++)
	    if ( d.masked(ik) ) ++mis; 

	  d = d.purge_rows();

	  Data::Matrix<double> cov = Statistics::covariance_matrix( d );
	  
	  double r = cov(0,1) / sqrt( cov(0,0) * cov(1,1) ) ; 
	  
	  const int na = d.dim1();

	  // Either/both

	  int g1 = 0 , g2 = 0;
	  for (int ii=0; ii<na; ii++)
	    {
	      if ( d(ii,0) ) ++g1;
	      if ( d(ii,1) ) ++g2;
	    }
	  
	  // Either/both (not pairwise missing)
	  int either = 0, both = 0;
	  for (int ii=0; ii<n; ii++) 
	    {
	      if ( ! ( alt1.masked(ii) || alt2.masked(ii) ) )
		{
		  either += alt1(ii) || alt2(ii);
		  both += alt1(ii) && alt2(ii);
		}
	    }
	 	  
	  // Non-ref genotype counts
	  
	  // FREQ1, OBS1
	  plog.data( n1 ) ; 
	  plog.data( n2 ) ; 

	  plog.data( g1 );
	  plog.data( g2 );
	  plog.data( na );	  

	  // Both/Either

	  plog.data( both  );
	  plog.data( either  );
	  if ( either ) 
	    plog.data( both/(double)either );
	  else 
	    plog.data( "NA" );

	  if ( Helper::realnum( r ) ) plog.data( r ) ; 
	  else plog.data( "NA" ) ; 	  
	  
	  plog.print_data_group();    

	  ++i;
	}
    }
  else // if not, we can purge
    {
      aux->store.clear();
    }
  
  aux->lastvar = var;

}

bool Pseq::VarDB::proximity_scan( Mask & mask )
{
  
  Aux_proximity aux;
  aux.dist = options.key( "dist" ) ? options.as<int>( "dist" ) : 2 ;
  
  plog.data_reset();	

  plog.data_group_header( "PAIR" );

  plog.data_header( "DIST" );

  plog.data_header( "FRQ1" );
  plog.data_header( "FRQ2" );

  plog.data_header( "GENO1" );
  plog.data_header( "GENO2" );

  plog.data_header( "NP" );
  plog.data_header( "BOTH" );
  plog.data_header( "EITHER" );
  plog.data_header( "CONC" );

  plog.data_header( "R" );

  plog.data_header_done();
  
  IterationReport report = g.vardb.iterate( f_proxscan , &aux , mask );
  
  return true;
}

