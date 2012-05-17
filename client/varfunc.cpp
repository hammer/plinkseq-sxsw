
#include "util.h"
#include "pseq.h"
#include "assoc.h"
#include "genic.h"

extern GStore g;
extern Pseq::Util::Options args;

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

		if ( ! ( g1->null() || g2->null() ) )
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

	    if ( ! j1->second->null() ) 
	      {
		record( v.ind(i)->id() , j1->first );
		++nobs;
		
		//		uniq_obs.insert( v.psample( j1->first )->label( *(j1->second) , true ) );
		
		
		// Insert unphased genotype label here
		uniq_obs.insert( v.geno_label( j1->first , *(j1->second) ) );

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
	   << v.geno_label( f1 , *g1 ) << "\t"
	   << v.geno_label( f2 , *g2 ) << "\t"
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

  if ( args.has("report-all") )
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




struct AuxLookup 
{ 
  bool append_phe;
  bool append_loc;
  bool append_aliases;
  bool append_ref;
  bool append_seq;
  bool vardb;
  bool append_annot;    
  std::set<std::string> locs;
  std::set<std::string> refs;
  std::set<std::string> aliases;
};


void f_lookup_annotator( Variant & var , void * p )
{
  
  std::cout << "in here..\n";

  AuxLookup * aux = (AuxLookup*)p;

  Region region( var );
  
  if ( aux->append_annot ) 
    {

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

  std::cout << "s2\n";

  std::string s = var.coordinate();

  std::cout << "s3\n";

  // Fetch from SEQDB
  
  if ( aux->append_seq ) 
    {
      if ( region.length() <= 10 )
	plog << s << "\t" 
	     << "seqdb_ref" << "\t"
	     << g.seqdb.lookup( region ) << "\n";	        
      else
	plog << s << "\t" 
	     << "seqdb_ref" << "\t"
	     << "." << "\n";	        

    }
  
  // Fetch from VARDB 

  std::set<Variant> vars = g.vardb.fetch( region );

  std::cout << "sX\n";

  if ( vars.size() == 0  )
    {
      if ( aux->vardb ) 
	{
	  
	  plog << s << "\t"
	       << "var" << "\t"
	       << "NA" << "\n";
	  
	  if ( aux->append_phe ) 
	    plog << s << "\t"
		 << "case" << "\t"
		 << "NA" << "\n"
		 << s << "\t"
		 << "con" << "\t"
		 << "NA" << "\n";
	}
    }

  plog << s << "\t"
       << "nvar" << "\t"
       << vars.size() << "\n";		    

  int cnt = 0;
  std::set<Variant>::iterator v = vars.begin();

  while ( v != vars.end() )
    {
      
      ++cnt;

      if ( vars.size() > 1 ) 
	plog << s << "\t"
	     << "var_" << cnt << "\t"
	     << *v << "\n";		    
      else
	plog << s << "\t"
	     << "var" << "\t"
	     << *v << "\n";		    	
      
      // 
      // Either allele counts; or stratify by case/control
      //
      
      int case_n = 0 , control_n = 0;
      
      for (int j=0; j < v->size() ; j++)
	{
	  
	  if ( (*v)(j).nonreference() )
	    {
	      if ( aux->append_phe ) 
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

      if ( aux->append_phe )
	{
	  if ( vars.size() > 1 ) 
	    plog << s << "\t"
		 << "case_" << cnt << "\t" 
		 << case_n << "\n" 
		 << s << "\t"
		 << "con_" << cnt << "\t" 
		 << control_n << "\n";
	  else
	    plog << s << "\t"
		 << "case" << "\t" 
		 << case_n << "\n" 
		 << s << "\t"
		 << "con" << "\t" 
		 << control_n << "\n";
	}
      else
	{
	  if ( vars.size() > 1 )
	    plog << s << "\t"
		 << "cnt_" << cnt << "\t" 
		 << control_n << "\n";
	  else
	    plog << s << "\t"
		 << "cnt" << "\t" 
		 << control_n << "\n";
	}

      ++v;
      
    }
  

  //
  // Locus DB ? 
  //
  
  if ( aux->append_loc ) 
    {
      std::set<std::string>::iterator i = aux->locs.begin();
      while ( i != aux->locs.end() )
	{
	  std::set<Region> rregs = g.locdb.get_regions( *i , region );
	  if ( rregs.size() == 0 ) 
	    {
	      plog << s << "\t"
		   << "loc_" << *i << "\t"
		   << "." << "\n";
	    }
	  
	  std::set<Region>::iterator j = rregs.begin();
	  while ( j != rregs.end() )
	    {
	      plog << s << "\t"
		   << "loc_" << *i << "\t"
		   << j->coordinate() << ":"
		   << j->name << "\n";

	      if ( aux->append_aliases )
		{
		  std::set<std::string>::iterator jj = aux->aliases.begin();
		  while ( jj != aux->aliases.end() )
		    {
		      
		      std::set<std::string> a = g.locdb.targetted_lookup_alias( j->name , *i , *jj );
		      std::set<std::string>::iterator ia = a.begin();
		      while ( ia != a.end() )
			{
			  plog << s << "\t"
			       << "loc_" << *jj << "\t"			       
			       << *ia << "\n";
			  ++ia;
			}
		      
		      ++jj;
		    }

		}
	      
	      ++j;
	    }
	  ++i;
	}
    }
  
      
  //
  // Reference variants? 
  //

  
  if ( aux->append_ref ) 
    {
      std::set<std::string>::iterator i = aux->refs.begin();
      while ( i != aux->refs.end() )
	{
	  
	  std::set<RefVariant> rvars = g.refdb.lookup( region , *i );
	  
	  if ( rvars.size() == 0 ) 
	    {
	      plog << s << "\t"
		   << "ref_" << *i << "\t"
		   << "." << "\t"
		   << "." << "\n";
	    }
	  
	  std::set<RefVariant>::iterator j = rvars.begin();
	  while ( j != rvars.end() )
	    {
	      plog << s << "\t"
		   << "ref_" << *i << "\t"
		   << *j << "\t"
		   << j->value() << "\n";
	      ++j;    
	    }
	  ++i;
	}
    }
  
}


bool Pseq::VarDB::lookup_list( const std::string & filename , 
			       Mask & mask , 
			       const std::vector<Region> * regs )
{
 
  AuxLookup aux;

  // From LOCDB, take set of gene-groups
  // From REFDB, take set of ref-variants

  aux.locs = args.get_set( "loc" );
  aux.refs = args.get_set( "ref" );
  aux.aliases = args.get_set( "alias" );

  aux.append_phe = g.vardb.attached() && g.inddb.attached() && g.phmap.type() == PHE_DICHOT;    
  aux.append_loc = g.locdb.attached() && aux.locs.size() > 0;
  aux.append_aliases = aux.append_loc && aux.aliases.size() > 0;
  aux.append_ref = g.refdb.attached() && aux.refs.size() > 0;
  aux.append_seq = g.seqdb.attached();
  aux.vardb = g.vardb.attached();
  aux.append_annot = g.seqdb.attached() && args.has( "annotate" ) ;

  if ( ! ( aux.vardb || aux.append_loc || aux.append_ref || aux.append_seq || aux.append_annot ) ) 
    Helper::halt("no information to append");
  
  if ( aux.append_annot ) 
    {
      std::string annot_transcripts = PLINKSeq::DEFAULT_LOC_GROUP() ;      
      if ( args.has( "annotate" ) ) annot_transcripts = args.as_string( "annotate" );      
      Annotate::load_transcripts( LOCDB, annot_transcripts );
    }
  

  // Headers

  if ( aux.vardb )
    {
      plog << "##nvar,1,Integer,\"Number of overlapping variants\"\n";
      plog << "##var,1,String,\"Variant ID\"\n";
      if ( aux.append_phe ) {
	plog << "##case,1,Integer,\"Case minor allele counts\"\n";
	plog << "##con,1,Integer,\"Control minor allele counts\"\n";
      }
      else plog << "##cnt,1,Integer,\"Minor allele counts\"\n";
    }
  
  if ( aux.append_loc )
    {
      std::set<std::string>::iterator i = aux.locs.begin();
      while ( i != aux.locs.end() )
	{
	  plog << "##" << "loc_" << *i << ",.,String,\"LOCDB group\"\n";
	  ++i;
	}
    }
  
  if ( aux.append_ref ) 
    {
      std::set<std::string>::iterator i = aux.refs.begin();
      while ( i != aux.refs.end() )
	{
	  plog << "##" << "ref_" << *i << ",.,String,\"REFDB group\"\n";
	  ++i;
	}
    }
  
  if ( aux.append_annot )
    {
      plog << "##func,1,String,\"Genomic annotation\"\n";
      plog << "##transcript,1,String,\"Transcript ID\"\n";
      plog << "##genomic,1,String,\"Genomic DNA change\"\n";
      plog << "##codon,1,String,\"Codon change\"\n";
      plog << "##protein,1,String,\"Any nonsynon amino acid change\"\n";
      plog << "##worst,1,String,\"Worst annotation\"\n";
      plog << "##class,1,String,\"Summary of all annotations\"\n";
    }
  
  if ( aux.append_seq )
    {
      plog << "##seqdb_ref,1,String,\"SEQDB reference sequence\"\n";
    }

    
  //
  // Annotate variants internally
  //
  
  if ( filename == "." && ! regs ) 
    {
      g.vardb.iterate( f_lookup_annotator , &aux , mask );
      return true;
    }


  //
  // Else, annotate variant list from file, or from a list of regions from the command line
  //

  if ( filename != "." )
    {
      Helper::checkFileExists( filename );

      InFile F1( filename );
      
      // Assume format is either a) a single site (1 col) OR b) single, REF and ALT (3 col)
      
      while ( !F1.eof() )
	{
	  
	  bool okay = true;
	  std::vector<std::string> line = F1.tokenizeLine( " \t" );
	  
	  if ( line.size() == 0 ) continue;
	  
	  if ( ! ( line.size() == 1 || line.size() == 3 ) ) 
	    {
	      plog.warn( "not a valid line of input" , line );
	      continue;
	    }
	  
	  // Region REF ALT
	  
	  Region region( line[0] , okay );     
	  std::string ref_allele = line.size() == 3 ? line[1] : ".";
	  std::string alt_allele = line.size() == 3 ? line[2] : ".";
	  
	  if ( okay ) 
	    {	  
	      Variant var;
	      var.chromosome( region.start.chromosome() );
	      var.position( region.start.position() );
	      
	      if ( ref_allele == "." )
		var.stop( region.stop.position() ); 
	      else
		var.stop( var.position() + ref_allele.size() - 1 );
	      
	      var.consensus.reference( ref_allele );
	      var.consensus.alternate( alt_allele );	  

	      f_lookup_annotator( var , &aux );

	    }
	  else
	    plog.warn( "not a valid region" , line[0] );
	  
	} // next input region

      F1.close();      
    }

  
  // And/or Regions from the command line (won't have allelic encoding)

  if ( regs ) 
    {
      for (int r = 0 ; r < regs->size(); r++)
	{
	  Variant var;
	  var.chromosome( (*regs)[r].start.chromosome() );
	  var.position(   (*regs)[r].start.position()   );
	  var.stop(       (*regs)[r].stop.position()    );
	  f_lookup_annotator( var , &aux );
	}      
    }
  

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
			args.has("show-subregions"));
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
  a.mhit = args.has( "assoc" , "mhit" );
  a.vt = args.has( "vt" );
  a.fw = args.has( "fw" );
  a.calpha = args.has( "calpha" );
  a.cancor = args.has( "cancor" );

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
  a.fix_null_genotypes = args.has( "fix-null" );


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
		  var.resize( ncase + ncontrol );
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
	  bool altmin1 = var.n_minor_allele( &n1, &m1 );
	  bool altmin2 = aux->lastvar.n_minor_allele( &n2, &m2 );
	  
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
  aux.dist = args.has( "distance" ) ? args.as_int( "distance" ) : 2 ;
  
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


//
//
//

struct aux_transmission_summary { 

  bool has_parents;

  int missing; // 1 or more missing genotypes in trio
  // PAT x MAT --> OFF 

  int aa_aa_aa;    
  int aa_aa_ab; // **    
  int aa_aa_bb; // **

  int aa_ab_aa;    
  int aa_ab_ab;    
  int aa_ab_bb; // **  
 
  int aa_bb_aa; // **
  int aa_bb_ab; 
  int aa_bb_bb; // **


  int ab_aa_aa;    
  int ab_aa_ab; 
  int ab_aa_bb; // **

  int ab_ab_aa;    
  int ab_ab_ab; 
  int ab_ab_bb; 
  
  int ab_bb_aa; // **    
  int ab_bb_ab; 
  int ab_bb_bb; 


  int bb_aa_aa; // **   
  int bb_aa_ab;     
  int bb_aa_bb; // **

  int bb_ab_aa; // **     
  int bb_ab_ab;    
  int bb_ab_bb; 
 
  int bb_bb_aa; // **
  int bb_bb_ab; // **
  int bb_bb_bb; 
  
  int dcount;

  int complete_transmissions() const
  {
    return aa_aa_aa +
      aa_aa_ab +
      aa_aa_bb +
      aa_ab_aa +
      aa_ab_ab +
      aa_ab_bb +
      aa_bb_aa +
      aa_bb_ab +
      aa_bb_bb +
      ab_aa_aa +
      ab_aa_ab +
      ab_aa_bb +
      ab_ab_aa +
      ab_ab_ab +
      ab_ab_bb +
      ab_bb_aa +   
      ab_bb_ab +
      ab_bb_bb +
      bb_aa_aa +  
      bb_aa_ab +
      bb_aa_bb +
      bb_ab_aa +    
      bb_ab_ab +
      bb_ab_bb +
      bb_bb_aa +
      bb_bb_ab +
      bb_bb_bb ;
  }

  int nonmendelian() const
  {
    return 
      aa_aa_ab + 
      aa_aa_bb +
      aa_ab_bb +
      aa_bb_aa + 
      aa_bb_bb +
      ab_aa_bb +
      ab_bb_aa +   
      bb_aa_aa +  
      bb_aa_bb +
      bb_ab_aa +    
      bb_bb_aa +
      bb_bb_ab ;
  }

  int trans_ref_from_het() const
  {
    return 
      aa_ab_aa +
      ab_aa_aa + 
      ab_ab_aa + 
      ab_ab_ab + 
      ab_bb_ab + 
      bb_ab_ab ;
  }


  int trans_alt_from_het() const
  {

    return 
      aa_ab_ab +
      ab_aa_ab +
      ab_ab_ab + 
      ab_ab_bb + 
      ab_bb_bb + 
      bb_ab_bb ;
  }

  int potential_denovo() const
  {
    return aa_aa_ab ;
  }

  int passing_denovo() const
  {
    return dcount;
  }

  aux_transmission_summary() 
  {
      aa_aa_aa =
      aa_aa_ab =
      aa_aa_bb =
      aa_ab_aa =
      aa_ab_ab =
      aa_ab_bb =
      aa_bb_aa =
      aa_bb_ab =
      aa_bb_bb = 0;
    
    ab_aa_aa =
      ab_aa_ab =
      ab_aa_bb =
      ab_ab_aa =
      ab_ab_ab =
      ab_ab_bb = 
      ab_bb_aa =
      ab_bb_ab =
      ab_bb_bb = 0;
    
    bb_aa_aa =
      bb_aa_ab =
      bb_aa_bb =
      bb_ab_aa = 
      bb_ab_ab =
      bb_ab_bb = 
      bb_bb_aa =
      bb_bb_ab =
      bb_bb_bb = 0;

    has_parents = false;
    missing = 0;
    dcount = 0;
  }

};


struct Aux_transmission_summary
{

  Aux_transmission_summary(const int n) 
    { 
      res.resize(n); 
      dp_kid = dp_par = 0;
      ab_kid_min = 0;
      ab_kid_max = 1 ;
      ab_par = 1;
      pl_kid = pl_par = 0;
    } 
  
  aux_transmission_summary * indiv(const int i) { return &res[i]; }
  
  std::vector<aux_transmission_summary> res;
  
  double dp_kid;
  double dp_par;
  
  double ab_kid_min, ab_kid_max;
  double ab_par;
  
  double pl_kid;
  double pl_par;
   
};

void f_denovo_scan( Variant & v , void * p )
{
  
  Aux_transmission_summary * aux = (Aux_transmission_summary*)p;
  
  const int n = v.size();

  for (int i=0; i<n; i++)
    {
      
      Individual * pat = v.ind(i)->pat();      
      Individual * mat = v.ind(i)->mat();
      
      // only consider individuals where both parents are present"
      if ( pat == NULL || mat == NULL ) continue;

      // there will be a better way to get the actual individual slot 
      // but use for now...
      
      int patn = g.indmap.ind_n( pat->id() );
      int matn = g.indmap.ind_n( mat->id() );
      Genotype & go = v(i);
      Genotype & gp = v(patn);
      Genotype & gm = v(matn);
      
      bool denovo = false;

      aux_transmission_summary * p = aux->indiv(i);
      
      if ( go.null() || gp.null() || gm.null() ) p->missing++;
      else
 	{
	  
 	  int ao = go.reference() ? 0 : go.heterozygote() ? 1 : go.alternate_homozygote() ? 2 : 0 ; 
 	  int ap = gp.reference() ? 0 : gp.heterozygote() ? 1 : gp.alternate_homozygote() ? 2 : 0 ; 
 	  int am = gm.reference() ? 0 : gm.heterozygote() ? 1 : gm.alternate_homozygote() ? 2 : 0 ; 
	  
 	  if ( ap == 0 ) 
 	    {
 	      if ( am == 0 ) 
 		{
 		  if      ( ao == 0 ) ++p->aa_aa_aa;
 		  else if ( ao == 1 ) { ++p->aa_aa_ab; denovo = true; }
 		  else if ( ao == 2 ) ++p->aa_aa_bb;
 		}
 	      else if ( am == 1 ) 
 		{
 		  if      ( ao == 0 ) ++p->aa_ab_aa;
 		  else if ( ao == 1 ) ++p->aa_ab_ab;
 		  else if ( ao == 2 ) ++p->aa_ab_bb;
 		}
 	      else // am == 2 
 		{
 		  if      ( ao == 0 ) ++p->aa_bb_aa;
 		  else if ( ao == 1 ) ++p->aa_bb_ab;
 		  else if ( ao == 2 ) ++p->aa_bb_bb;
 		}
 	    }
 	  else if ( ap == 1 ) 
 	    {
 	      if ( am == 0 ) 
 		{
 		  if      ( ao == 0 ) ++p->ab_aa_aa;
 		  else if ( ao == 1 ) ++p->ab_aa_ab;
 		  else if ( ao == 2 ) ++p->ab_aa_bb;
 		}
 	      else if ( am == 1 ) 
 		{
 		  if      ( ao == 0 ) ++p->ab_ab_aa;
 		  else if ( ao == 1 ) ++p->ab_ab_ab;
 		  else if ( ao == 2 ) ++p->ab_ab_bb;
 		}
 	      else // am == 2 
 		{
 		  if      ( ao == 0 ) ++p->ab_bb_aa;
 		  else if ( ao == 1 ) ++p->ab_bb_ab;
 		  else if ( ao == 2 ) ++p->ab_bb_bb;
 		}
 	    }
 	  else
 	    {
 	      if ( am == 0 ) 
 		{
 		  if      ( ao == 0 ) ++p->bb_aa_aa;
 		  else if ( ao == 1 ) ++p->bb_aa_ab;
 		  else if ( ao == 2 ) ++p->bb_aa_bb;
 		}
 	      else if ( am == 1 ) 
 		{
 		  if      ( ao == 0 ) ++p->bb_ab_aa;
 		  else if ( ao == 1 ) ++p->bb_ab_ab;
 		  else if ( ao == 2 ) ++p->bb_ab_bb;
 		}
 	      else // am == 2 
 		{
 		  if      ( ao == 0 ) ++p->bb_bb_aa;
 		  else if ( ao == 1 ) ++p->bb_bb_ab;
 		  else if ( ao == 2 ) ++p->bb_bb_bb;
 		}
 	    }
	  
 	}
      
      
      // would this putative de novo pass special denovo filters?
      if ( denovo ) 
	{

	  // Assume individual metrics: DP, PL, AD

	  if ( aux->dp_kid > 0 ) 
	    {
	      if ( ! go.meta.has_field( "DP" ) ) denovo = false; 
	      else if ( go.meta.get1_int( "DP" ) < aux->dp_kid ) denovo = false;
	    }

	  if ( aux->dp_par > 0 ) 
	    {
	      if ( ! gp.meta.has_field( "DP" ) ) denovo = false; 
	      else if ( gp.meta.get1_int( "DP" ) < aux->dp_par ) denovo = false;
	      
	      if ( ! gm.meta.has_field( "DP" ) ) denovo = false; 
	      else if ( gm.meta.get1_int( "DP" ) < aux->dp_par ) denovo = false;
	    }
	  
	  
	  // PLs
	  
	  if ( aux->pl_kid > 0 ) 
	    {
	      if ( ! go.meta.has_field( "PL" ) ) denovo = false; 
	      else 
		{
		  std::vector<int> pl = go.meta.get_int( "PL" ) ;
		  if ( pl.size() != 3 ) denovo = false;
		  else if ( pl[0] < aux->pl_kid || pl[2] < aux->pl_kid ) denovo = false;
		}
	    }
	  
	  if ( aux->pl_par > 0 ) 
	    {
	      
	      if ( ! gp.meta.has_field( "PL" ) ) denovo = false; 
	      else 
		{
		  std::vector<int> pl = gp.meta.get_int( "PL" ) ;
		  if ( pl.size() != 3 ) denovo = false;
		  else if ( pl[1] < aux->pl_par || pl[2] < aux->pl_par ) denovo = false;
		}

	      if ( ! gm.meta.has_field( "PL" ) ) denovo = false; 
	      else 
		{
		  std::vector<int> pl = gm.meta.get_int( "PL" ) ;
		  if ( pl.size() != 3 ) denovo = false;
		  else if ( pl[1] < aux->pl_par || pl[2] < aux->pl_par ) denovo = false;
		}
	    
	    }


	  // ABs
	  
	  if ( aux->ab_kid_min > 0 || aux->ab_kid_max < 1 ) 
	    {
	      if ( ! go.meta.has_field( "AD" ) ) denovo = false;
	      else 
		{
		  std::vector<int> ad = go.meta.get_int( "AD" );
		  if ( ad.size() != 2 ) denovo = false;
		  else
		    {
		      // prop. of NR reads
		      double ab = ad[1] / (double)( ad[0] + ad[1] );
		      if ( ab < aux->ab_kid_min || ab > aux->ab_kid_max ) denovo = false;
		    }
		}
	    }

	  // AB in parents
	  if ( aux->ab_par < 1 )
	    {
	      if ( ! gp.meta.has_field( "AD" ) ) denovo = false;
	      else 
		{
		  std::vector<int> ad = gp.meta.get_int( "AD" );
		  if ( ad.size() != 2 ) denovo = false;
		  else
		    {
		      double ab = ad[1] / (double)( ad[0] + ad[1] );
		      if ( ab > aux->ab_par ) denovo = false;
		    }
		}

	      if ( ! gm.meta.has_field( "AD" ) ) denovo = false;
	      else 
		{
		  std::vector<int> ad = gm.meta.get_int( "AD" );
		  if ( ad.size() != 2 ) denovo = false;
		  else
		    {
		      double ab = ad[1] / (double)( ad[0] + ad[1] );
		      if ( ab > aux->ab_par ) denovo = false;
		    }
		}
	    }
	}
      
      
      // directly output possible de novo events (only REF x REF --> HET)
      // that also passed any above, de-novo specific filters
      
      if ( denovo ) 
	{
	  
	  // track # of actual 'passing' de novo calls
	  p->dcount++;
	  
	  // get allele frequencies
	  int c = 0 , c_tot = 0;
	  v.n_minor_allele( &c , &c_tot ); 
	  
	  plog << "Variant\tRefxRef->Het\t" 
	       << v << "\t" 
	       << c << "\t"
	       << c_tot << "\t"
	       << v.ind(i)->id() << "\t"
	       << v.label( patn , "," ) << " x "
	       << v.label( matn , "," ) << " -> "
	       << v.label( i , "," ) << "\t"
	       << "[" << v.gmeta_label(patn) << "]" << "\t"
	       << "[" << v.gmeta_label(matn) << "]" << "\t"
	       << "[" << v.gmeta_label(i) << "]" << "\n";

	}      
    }
}


bool Pseq::VarDB::denovo_scan( Mask & mask )
{
  
  // did we have some special values 
  
  const int n = g.indmap.size();

  // store summary transmission data 
  Aux_transmission_summary aux(n);

  if ( args.has("param") )
    {
      std::vector<double> p = args.as_float_vector( "param" );
      if ( p.size() != 7 ) 
	Helper::halt( "expect --param DP(kid) DP(par) PL(kid) PL(par) AB(kid,lwr) AB(kid,upr) AB(par,upr)" );

      aux.dp_kid = p[0];
      aux.dp_par = p[1];
      aux.pl_kid = p[2];
      aux.pl_par = p[3];
  
      aux.ab_kid_min = p[4];
      aux.ab_kid_max = p[5];
      aux.ab_par = p[6];  

    }

  
  // Attach parents
  for (int i=0;i<n;i++)
    {      
      Individual * person = g.indmap(i);
      g.inddb.fetch( person );
      Individual * p = g.indmap.ind( person->father() );
      Individual * m = g.indmap.ind( person->mother() );
      if ( p ) person->pat( p );
      if ( m ) person->mat( m );
      if ( p && m ) aux.indiv(i)->has_parents = true;      
    }
  
  g.vardb.iterate( f_denovo_scan , &aux , mask );
  
  // display summaries
  
  for (int i=0;i<n;i++)
    {
      aux_transmission_summary * p = aux.indiv(i);
      if ( p->has_parents ) 
	plog << "Individual\t" 
	     << g.indmap(i)->id() << "\t" 
	     << p->complete_transmissions() << "\t"
	     << p->missing << "\t"
	     << p->trans_ref_from_het() << "\t"
	     << p->trans_alt_from_het() << "\t"
	     << p->trans_alt_from_het() / (double)( p->trans_ref_from_het() + p->trans_alt_from_het() ) << "\t"
	     << p->nonmendelian() << "\t"
	     << p->potential_denovo() << "\t"
	     << p->passing_denovo() << "\n";
    }
  
  return true;
}




bool Pseq::VarDB::write_lookup_matrix(Mask & m , const std::string & filename, const std::vector<std::string> & regs )
{

  // replace this functionality with a mask option instead ( reg.force=... )

//   // ultimately, might be better to adapt the reg mask function to
//   // have an option to force a position even when not at all present
//   // in the dataset?  For now, this will still be useful.
  
//   std::set<Region> lookups;
//   if ( filename != "." && Helper::fileExists( filename ) )
//     {
//       InFile IN1( filename );
//       while ( ! IN1.eof() )
// 	{
// 	  std::string s;
// 	  IN1 >> s;
// 	  bool okay = false;
// 	  Region r( s , okay );
// 	  if ( okay ) lookups.insert(r);
// 	}
//       IN1.close();
//     }
  
//   for (int i=0;i<regs.size();i++)
//     {
//       bool okay = false;
//       Region r( regs[i] , okay );
//       if ( okay ) lookups.insert(r);
//     }
  
//   std::cout << "found " << lookups.size() << " regions to lookup\n";
 

return false;
 }
