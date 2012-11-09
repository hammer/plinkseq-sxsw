
#include "util.h"
#include "plinkseq.h"
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
		
		// uniq_obs.insert( v.psample( j1->first )->label( *(j1->second) , true ) );
				
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
	
	Out & pvar = Out::stream( "concord.vars" ); 
	
	pvar << v << "\t" 
	     << v.n_samples() << "\t"
	     << n_multiple_obs << "\t"
	     << n_concordant_obs << "\t";
	
	if ( n_multiple_obs > 1 )
	  pvar << (double)(n_concordant_obs)/(double)n_multiple_obs << "\t";
	else
	  pvar << "NA\t";

	if ( v.n_samples() == 2 ) 
	  pvar << cnt_ra << "\t" 
	       << cnt_ar << "\t"
	       << cnt_aa << "\t";
 	else
	  pvar << "NA" << "\t" 
	       << "NA" << "\t"
	       << "NA" << "\t";
	
	if ( disc_indiv.size() == 0 ) 
	  pvar << ".";
	else
	  {
	    std::set<std::string>::iterator i = disc_indiv.begin();
	    while ( i != disc_indiv.end() )
	      {
		if ( i != disc_indiv.begin() ) pvar << ",";
		pvar << *i;
		++i;
	      }
	  }
	pvar << "\n";
      }
  }


  void report_helper( AuxConcordanceStats & stats , const std::string & label )
  {

    // if label == "" implies total sample
    // otherwise an individual ID
    
    bool whole_sample = label == "";

    Out & pout = Out::stream(  whole_sample ? "concord" : "concord.indiv" );
    
    std::map<int,int>::iterator i = stats.once.begin();
    while ( i != stats.once.end() )
      {
	std::map<int,int>::iterator j = i;
	while ( j != stats.once.end() )
	  {
	    
	    if ( i == j ) { ++j; continue; }
	    
	    pout << (whole_sample ? "TOTAL" : label ) << "\t"
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
	    
	    pout << disc << "\t" 
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
	report_helper( ii->second , ii->first );
	++ii;
      }
	
    	
    //
    // for all files seen, give total variant count
    //

    report_helper( tstats , "" );
    
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
      {
	Out & pgeno = Out::stream( "concord.geno" );

	pgeno << id << "\t"
	      << f1 << "\t"
	      << f2 << "\t"
	      << v << "\t"
	      << v.geno_label( f1 , *g1 ) << "\t"
	      << v.geno_label( f2 , *g2 ) << "\t"
	      << g1->nonreference() << "\t"
	      << g2->nonreference() << "\n";

      }

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

  Out & pvar = Out::stream( "concord.vars" );

  pvar << "VAR" << "\t" 
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
  bool append_prot;
  bool append_aliases;
  bool append_ref;
  bool append_ref_allelic;
  bool append_seq;
  bool vardb;
  bool append_annot;
  bool append_titv;

  ProtDBase protdb;
  std::set<std::string> locs;
  std::set<std::string> refs;
  std::set<std::string> refs_allelic;
  std::set<std::string> aliases;

};


void f_lookup_annotator( Variant & var , void * p )
{
  

  AuxLookup * aux = (AuxLookup*)p;

  Region region( var );

  Out & pout = Out::stream( "meta" );



  std::string s = var.coordinate();
  
  pout << s << "\t"
       << "allele_ref\t"
       << var.reference() << "\n"
       << s << "\t"
       << "allele_alt\t"
       << var.alternate() << "\n";

  //
  // Fetch from SEQDB
  //

  if ( aux->append_seq ) 
    {
      if ( region.length() <= 10 )
	pout << s << "\t" 
	     << "seqdb_ref" << "\t"
	     << g.seqdb.lookup( region ) << "\n";	        
      else
	pout << s << "\t" 
	     << "seqdb_ref" << "\t"
	     << "." << "\n";	        

    }


  
  if ( aux->append_annot ) 
    {

      bool exonic = Annotate::annotate( var );
  

      std::string annot = var.meta.get1_string( PLINKSeq::ANNOT() );
	  
      // detailed annotation vector, primary annotation
  
      pout << var.coordinate() << "\t"
	   << "func" << "\t"
	   << var.meta.as_string( PLINKSeq::ANNOT_TYPE() , "," ) << "\n";
      
      pout << var.coordinate() << "\t"
	   << "transcript" << "\t"
	   << var.meta.as_string( PLINKSeq::ANNOT_GENE() , "," ) << "\n";
      
      pout << var.coordinate() << "\t"
	   << "genomic" << "\t"
	   << var.meta.as_string( PLINKSeq::ANNOT_CHANGE() , "," ) << "\n";
      
      pout << var.coordinate() << "\t"
	   << "codon" << "\t"
	   << var.meta.as_string( PLINKSeq::ANNOT_CODON() , "," ) << "\n";
      
      pout << var.coordinate() << "\t"
	   << "protein" << "\t"
	   << var.meta.as_string( PLINKSeq::ANNOT_PROTEIN() , "," ) << "\n";
      
      // worst-case consensus annotation
      
      pout << var.coordinate() << "\t"
	   << "worst" << "\t"
	   << annot << "\n";
      
      // summary annotation
      
      pout << var.coordinate() << "\t"
	   << "class" << "\t"
	   << var.meta.get1_string( PLINKSeq::ANNOT_SUMMARY()) << "\n";
      
    }


  
  // Fetch from VARDB 
  
  std::set<Variant> vars = g.vardb.fetch( region );
  
  if ( vars.size() == 0  )
    {
      
      if ( aux->vardb ) 
	{
	  
	  pout << s << "\t"
	       << "var" << "\t"
	       << "NA" << "\n";
	  
	  if ( aux->append_phe ) 
	    pout << s << "\t"
		 << "case" << "\t"
		 << "NA" << "\n"
		 << s << "\t"
		 << "con" << "\t"
		 << "NA" << "\n";
	}
    }

  
  pout << s << "\t"
       << "nvar" << "\t"
       << vars.size() << "\n";		    
  

  int cnt = 0;
  std::set<Variant>::iterator v = vars.begin();
  
  std::stringstream ss_var, ss_case, ss_con;
  bool any_var = false, any_cc = false;
  
  while ( v != vars.end() )
    {
      
      ++cnt;
      
      any_var = true;
      
      if ( cnt > 1 )
	ss_var << "," << *v << ":" << v->reference() << "/" << v->alternate() ;
      else 
	ss_var << *v << ":" << v->reference() << "/" << v->alternate() ;
      
      
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
		    case_n += (*v)(j).minor_allele_count( true );  // true implies AAC, not MAC
		  else
		    control_n += (*v)(j).minor_allele_count( true );
		}
	      else
		control_n += (*v)(j).minor_allele_count( true );
	    }
	}
      
      if ( aux->append_phe )
	{

	  any_cc = true;
	  
	  if ( cnt > 1 )
	    {
	      ss_case << "," << case_n;
	      ss_con << "," << control_n;
	    }
	  else
	    {
	      ss_case << case_n;
	      ss_con << control_n;
	    }

	}
      else
	{
	  if ( cnt > 1 )
	    ss_con << "," << control_n;
	  else
	    ss_con << control_n;	 
	}

      ++v;
      
    }

  
  if ( any_var )    
    pout << s << "\t"
	 << "var" << "\t"
	 << ss_var.str() << "\n";
  
  if ( any_cc )
    {

      if ( aux->append_phe )
	pout << s << "\t"
	     << "case" << "\t"
	     << ss_case.str() << "\n"      
	  
	     << s << "\t"
	     << "con" << "\t"
	     << ss_con.str() << "\n";
      else
	pout << s << "\t"
	     << "cnt" << "\t"
	     << ss_con.str() << "\n";
    }


  if (aux->append_titv) {
	  std::string titv = ".";
	  if ( var.transition() )
		  titv = "transition";
	  else if (var.transversion())
		  titv = "transversion";

	  pout << s << "\t"
	       << "titv" << "\t"
	       << titv << "\n";
  }

  if ( aux->append_prot ) 
    {
      ProtFeatureSet fs = aux->protdb.lookup( var );

      // consider each transcript

      std::map<std::string,std::set<Feature> >::iterator ii = fs.feat.begin();
      while ( ii != fs.feat.end() )
	{
	  
	  // given this particular transcript, find the AA position
	  
	  int aapos = 0;
	  int bp = var.position();

	  bool negative_strand = false;
	  bool positive_strand = false; 
	  bool unknown_strand = true;
	  int aasize = 0;

	  // any features? if so, calcualte AA position for variant
	  
	  Region * region = Annotate::pointer_to_region( ii->first );
	  if ( region == NULL )
	  {
		  aapos = 0;
	  }
	  else
	  {

		  // determine the genomic position for the start/stop pairs
		  int ns = region->subregion.size();

		  // transcript on positive or negative strand?

		  for (int e = 0; e < ns ; e++ )
		  {
			  if ( region->subregion[e].CDS() )
			  {
				  aasize += region->subregion[e].stop.position() - region->subregion[e].start.position() + 1;
				  int s = region->subregion[e].meta.get1_int( PLINKSeq::TRANSCRIPT_STRAND() );
				  if ( s == 0 ) continue;
				  negative_strand = s < 0 ;
				  positive_strand = s > 0 ;
			  }

		  }

		  bool problem = false;
		  if ( negative_strand != positive_strand )
		  {
			  unknown_strand = false;
			  if ( aasize % 3 != 0 )
			  {
				  std::cerr << "problem, " << aasize << " " << aasize % 3 << "\n";
				  problem = true;
			  }
		  }


		  //
		  // this procedure also implicitly makes sure that the variant is in the CDS
		  // (i.e. not intronic for variants that span more than one exon
		  //

		  bool set_pos = false;

		  if (negative_strand)
		  {
			  for (int e = ns-1; e >= 0 ; e-- )
			  {
				  // only consider actual CDS regions
				  if ( region->subregion[e].CDS() )
				  {
					  // is the position in this region?
					  if ( bp >= region->subregion[e].start.position() && bp <= region->subregion[e].stop.position() )
					  {
						  aapos += region->subregion[e].stop.position() - bp;  // 0-based
						  set_pos = true;
						  break;
					  }
					  aapos += region->subregion[e].stop.position() - region->subregion[e].start.position() + 1;
				  }
			  }
		  }
		  else {
			  for (int e = 0; e < ns ; e++ )
			  {
				  // only consider actual CDS regions
				  if ( region->subregion[e].CDS() )
				  {
					  // is the position in this region?
					  if ( bp >= region->subregion[e].start.position() && bp <= region->subregion[e].stop.position() )
					  {
						  aapos += bp - region->subregion[e].start.position() ;  // 0-based
						  set_pos = true;
						  break;
					  }
					  aapos += region->subregion[e].stop.position() - region->subregion[e].start.position() + 1;
				  }
			  }
		  }


		  if ( problem || ! set_pos )
			  aapos = 0;
		  else
		  {
			  aapos = ( (int)aapos / (int)3 ) + 1; // 1-based AA number
		  }
	  }
	
      
      // consider all features, if maps to CDS
      
      if ( aapos ) 
	{
	  std::set<Feature>::iterator jj = ii->second.begin();	  
	  while ( jj != ii->second.end() )
	    {
	      pout << s << "\t"
		   << "prot_" << jj->source_id << "\t"
		   << jj->feature_id << "|" << jj->feature_name << "|" << jj->protein_id << "|";
	      pout << aapos << "/";   
	      pout << jj->pstart << ".." << jj->pstop << "|" << jj->mstr << "|" << ii->first << ( unknown_strand ? "(?strand)" : negative_strand ? "(-ve)" : "(+ve)" ) << "\n";
	      ++jj;
	    }
	}
      
      // next transcript
      ++ii;
    }
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
	      pout << s << "\t"
		   << "loc_" << *i << "\t"
		   << "." << "\n";
	    }
	  
	  std::set<Region>::iterator j = rregs.begin();
	  while ( j != rregs.end() )
	    {
	      pout << s << "\t"
		   << "loc_" << *i << "\t"
		//<< j->coordinate() << ":"
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
			  pout << s << "\t"
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
	      pout << s << "\t"
		   << "ref_" << *i << "\t"
		   << "." << "\t"
		   << "." << "\n";
	    }
	  
	  std::set<RefVariant>::iterator j = rvars.begin();
	  while ( j != rvars.end() )
	    {
	      pout << s << "\t"
		   << "ref_" << *i << "\t"
		   << *j << "\t"
		   << j->value() << "\n";
	      ++j;    
	    }
	  ++i;
	}
    }
  
  if (aux->append_ref_allelic) {
      std::set<std::string>::iterator i = aux->refs_allelic.begin();
      while ( i != aux->refs_allelic.end() )
	{

	  std::set<RefVariant> rvars = g.refdb.lookup( var , *i , true );

	  if ( rvars.size() == 0 )
	    {
	      pout << s << "\t"
		   << "ref_allelic_" << *i << "\t"
		   << "." << "\t"
		   << "." << "\n";
	    }

	  std::set<RefVariant>::iterator j = rvars.begin();
	  while ( j != rvars.end() )
	    {
	      pout << s << "\t"
		   << "ref_allelic_" << *i << "\t"
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

  Out & pout = Out::stream( "meta" );

  // From LOCDB, take set of gene-groups
  // From REFDB, take set of ref-variants

  aux.locs = args.get_set( "loc" );
  aux.refs = args.get_set( "ref" );
  aux.refs_allelic = args.get_set( "ref_allelic" );
  aux.aliases = args.get_set( "alias" );

  aux.append_phe = g.vardb.attached() && g.inddb.attached() && g.phmap.type() == PHE_DICHOT;    
  aux.append_loc = g.locdb.attached() && aux.locs.size() > 0;
  aux.append_prot = args.has( "protdb" );

  if ( aux.append_prot && g.locdb.attached() && ! aux.append_annot )
    {      
      // NOTE: assumes 'refseq' is named group in both PROTDB and LOCDB for now      
      plog << "assuming group 'refseq' exists in PROTDB and LOCDB\n";
      Annotate::setDB( LOCDB );
      if ( ! Annotate::set_transcript_group( PLINKSeq::DEFAULT_LOC_GROUP() ) ) Helper::halt( "trouble attaching 'refseq' group from LOCDB" );    
    }
  
  aux.append_aliases = aux.append_loc && aux.aliases.size() > 0;
  aux.append_ref = g.refdb.attached() && aux.refs.size() > 0;
  aux.append_ref_allelic = g.refdb.attached() && aux.refs_allelic.size() > 0;
  aux.append_seq = g.seqdb.attached();
  aux.vardb = g.vardb.attached();
  aux.append_annot = g.seqdb.attached() && args.has( "annotate" ) ;
  aux.append_titv = args.has( "titv" );
  
  if ( ! ( aux.vardb || aux.append_loc || aux.append_prot || aux.append_ref || aux.append_ref_allelic || aux.append_seq || aux.append_annot ) )
    Helper::halt("no information to append");
  
  
  if ( aux.append_annot ) 
    {
      std::string annot_transcripts = PLINKSeq::DEFAULT_LOC_GROUP() ;      
      if ( args.has( "annotate" ) ) annot_transcripts = args.as_string( "annotate" );      
      Annotate::setDB( LOCDB );
      Annotate::set_transcript_group( annot_transcripts );
    }
  

  // Headers
  
  if ( aux.vardb )
    {
      pout << "##nvar,1,Integer,\"Number of overlapping variants\"\n";
      pout << "##var,1,String,\"Variant ID\"\n";
      if ( aux.append_phe ) {
	pout << "##case,1,Integer,\"Case alternate allele counts\"\n";
	pout << "##con,1,Integer,\"Control alternate allele counts\"\n";
      }
      else pout << "##cnt,1,Integer,\"Alternate allele counts\"\n";
    }

  

   
  if ( aux.append_prot )
    {
      
      aux.protdb.attach( args.as_string( "protdb" ) );
      
      if ( ! aux.protdb.attached() ) Helper::halt( "could not attach PROTDB" );
      
      std::set<std::string> sources = aux.protdb.get_sources();
      std::set<std::string>::iterator ii = sources.begin();
      while ( ii != sources.end() )
	{
	  pout << "##prot_" << *ii << ",1,String,\"PROTDB " << *ii << " annotation\"\n";
	  ++ii;
	}
    }


  if ( aux.append_loc )
    {
      std::set<std::string>::iterator i = aux.locs.begin();
      while ( i != aux.locs.end() )
	{
	  pout << "##" << "loc_" << *i << ",.,String,\"LOCDB group\"\n";
	  ++i;
	}
    }
  

  if ( aux.append_aliases )
    {

      std::set<std::string>::iterator jj = aux.aliases.begin();
      while ( jj != aux.aliases.end() )
	{	  
	  pout << "##" << "loc_" << *jj << ",.,String,\"LOCDB alias group\"\n";
	  ++jj;
	}
    }


  if ( aux.append_ref ) 
    {
      std::set<std::string>::iterator i = aux.refs.begin();
      while ( i != aux.refs.end() )
	{
	  pout << "##" << "ref_" << *i << ",.,String,\"REFDB group\"\n";
	  ++i;
	}
    }

  if ( aux.append_ref_allelic )
    {
      std::set<std::string>::iterator i = aux.refs_allelic.begin();
      while ( i != aux.refs_allelic.end() )
	{
	  pout << "##" << "ref_allelic_" << *i << ",.,String,\"REFDB group\"\n";
	  ++i;
	}
    }

  if ( aux.append_annot )
    {
      pout << "##func,.,String,\"Genomic annotation\"\n";
      pout << "##transcript,.,String,\"Transcript ID\"\n";
      pout << "##genomic,.,String,\"Genomic DNA change\"\n";
      pout << "##codon,.,String,\"Codon change\"\n";
      pout << "##protein,.,String,\"Any nonsynon amino acid change\"\n";
      pout << "##worst,1,String,\"Worst annotation\"\n";
      pout << "##class,.,String,\"Summary of all annotations\"\n";
    }
  
  pout << "##allele_ref,1,String,\"VCF reference sequence\"\n"
       << "##allele_alt,.,String,\"VCF alternate sequence(s)\"\n";

  if ( aux.append_seq )
    {
      pout << "##seqdb_ref,1,String,\"SEQDB reference sequence\"\n";
    }

  if (aux.append_titv) {
      pout << "##titv,1,String,\"Transition, Transversion, or neither?\"\n";
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
	  
	  if ( line.size() == 0 || line[0].substr(0,1) == "#") continue;
	  
	  if ( ! ( line.size() == 1 || line.size() >= 3 ) )
	    {
	      plog.warn( "not a valid line of input" , line );
	      continue;
	    }
	  
	  // Region REF ALT
	  
	  Region region( line[0] , okay );     
	  std::string ref_allele = line.size() >= 3 ? line[1] : ".";
	  std::string alt_allele = line.size() >= 3 ? line[2] : ".";
	  
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

	      // Call make_consensus() so that the alleles vector is parsed out into individual alleles:
	      IndividualMap a;
	      //a.populate(std::vector<std::string>());
	      var.make_consensus(&a);

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
  

  //
  // clean-up
  //

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
  
  Out & pout = Out::stream( "loci" );
  
  if ( loc.size() == 0 ) 
    {
      pout << v << "\t"
	   << 0 << "\t"
	   << "." << "\t"
	   << "." << "\t"
	   << ( aux->show_sub ? "." : "" ) << "\t";
      for (int j=0;j<aux->nm;j++) pout << (j!=0 ? "\t" : "" ) << ".";
      pout << "\n";
    }
  else
    {
      int cnt = 1;
      std::set<Region>::iterator i = loc.begin();
      while ( i != loc.end() ) 
	{
	  pout << v << "\t"
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
  
  Out & pout = Out::stream( "loci" );

  pout << "VAR" << "\t"
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

  Out & pout = Out::stream( "assoc.sim" );

  // a very basic wrapper purely to generate some data to be fed to 
  // to various genic tests, primarily in a sanity-check/debug mode
  // although could also be used for power calculations, etc

  Pseq::Assoc::AuxGenic a;
  a.g     = &g;
  a.rseed = time(0);

  pout.data_reset();  
  pout.data_group_header( "LOCUS" );
  pout.data_header( "POS" );
  pout.data_header( "ALIAS" );
  pout.data_header( "NVAR" );
  pout.data_header( "TEST" );
  pout.data_header( "P" );
  pout.data_header( "I" );
  pout.data_header( "DESC" );  
  pout.data_header_done();
  
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

  Out & pout = Out::stream( "proximity" );

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

	  pout.data_group( i->displaycore() + "," + var.displaycore() );
	  
	  // Counts

	  int n1, m1, n2, m2;
	  bool altmin1 = var.n_minor_allele( &n1, &m1 );
	  bool altmin2 = aux->lastvar.n_minor_allele( &n2, &m2 );
	  
	  // DIST
	  pout.data( var.position() - aux->lastvar.stop() );

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
	  pout.data( n1 ) ; 
	  pout.data( n2 ) ; 

	  pout.data( g1 );
	  pout.data( g2 );
	  pout.data( na );	  

	  // Both/Either

	  pout.data( both  );
	  pout.data( either  );
	  if ( either ) 
	    pout.data( both/(double)either );
	  else 
	    pout.data( "NA" );

	  if ( Helper::realnum( r ) ) pout.data( r ) ; 
	  else pout.data( "NA" ) ; 	  
	  
	  pout.print_data_group();    

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

  Out & pout = Out::stream( "proximity" );

  Aux_proximity aux;

  aux.dist = args.has( "distance" ) ? args.as_int( "distance" ) : 2 ;
  
  pout.data_reset();	

  pout.data_group_header( "PAIR" );

  pout.data_header( "DIST" );

  pout.data_header( "FRQ1" );
  pout.data_header( "FRQ2" );

  pout.data_header( "GENO1" );
  pout.data_header( "GENO2" );

  pout.data_header( "NP" );
  pout.data_header( "BOTH" );
  pout.data_header( "EITHER" );
  pout.data_header( "CONC" );

  pout.data_header( "R" );

  pout.data_header_done();
  
  IterationReport report = g.vardb.iterate( f_proxscan , &aux , mask );
  
  return true;
}


//
//
//

struct aux_trio_transmission_summary { 

  bool has_parents;

  int _missing; // 1 or more missing (or non-diploid) genotypes in trio
  int _complete_transmissions;
  int _nonmendelian;
  int _trans_ref_from_het;
  int _trans_alt_from_het;
  int _potential_denovo;
  int _dcount;

  int missing() const
  {
	return _missing;
  }

  int complete_transmissions() const
  {
    return _complete_transmissions;
  }

  int nonmendelian() const
  {
    return _nonmendelian;
  }

  int trans_ref_from_het() const
  {
    return _trans_ref_from_het;
  }


  int trans_alt_from_het() const
  {

    return _trans_alt_from_het;
  }

  int potential_denovo() const
  {
    return _potential_denovo;
  }

  int passing_denovo() const
  {
    return _dcount;
  }

  aux_trio_transmission_summary() 
  {
	  _missing =
	  _complete_transmissions =
	  _nonmendelian =
	  _trans_ref_from_het =
	  _trans_alt_from_het =
	  _potential_denovo =
	  _dcount = 0;

    has_parents = false;
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

      printTransmission = false;
    } 
  
  aux_trio_transmission_summary * indiv(const int i) { return &res[i]; }
  
  std::vector<aux_trio_transmission_summary> res;
  
  double dp_kid;
  double dp_par;
  
  double ab_kid_min, ab_kid_max;
  double ab_par;
  
  double pl_kid;
  double pl_par;
   
  bool printTransmission;
};

#define MULTI_ALLELIC_DELIM "#"

#define VAR_DATA \
		v << "\t" \
		<< v.label( patn , MULTI_ALLELIC_DELIM ) << "\t" \
		<< v.label( matn , MULTI_ALLELIC_DELIM ) << "\t" \
		<< v.label( i , MULTI_ALLELIC_DELIM ) << "\t" \
		<< v.ind(patn)->id() << "\t" \
		<< v.ind(matn)->id() << "\t" \
		<< v.ind(i)->id() << "\t" \
		<< c << "\t" \
		<< c_tot << "\t" \
		<< "[" << v.gmeta_label(patn, MULTI_ALLELIC_DELIM) << "]" << "\t" \
		<< "[" << v.gmeta_label(matn, MULTI_ALLELIC_DELIM) << "]" << "\t" \
		<< "[" << v.gmeta_label(i, MULTI_ALLELIC_DELIM) << "]" << "\t" \
		<< v.alternate() << "\t" \
		<< v.alternate_label( patn , MULTI_ALLELIC_DELIM ) << "\t" \
		<< v.alternate_label( matn , MULTI_ALLELIC_DELIM ) << "\t" \
		<< v.alternate_label( i , MULTI_ALLELIC_DELIM )

void f_denovo_scan( Variant & v , void * p )
{

  Out & pDeNovos = Out::stream( "denovo.vars" );

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

      // PAT x MAT --> OFFSPRING:
      Genotype & go = v(i);
      Genotype & gp = v(patn);
      Genotype & gm = v(matn);

      aux_trio_transmission_summary* summary = aux->indiv(i);
      aux_trio_transmission_summary prevSummary = *summary;

      if ( go.null() || gp.null() || gm.null() || !go.diploid() || !gp.diploid() || !gm.diploid() ) {
    	  summary->_missing++;
      }
      else
 	{
    	  summary->_complete_transmissions++;

    	  int a1 = go.acode1();
    	  int a2 = go.acode2();

    	  std::set<int> pAll;
    	  pAll.insert(gp.acode1());
    	  pAll.insert(gp.acode2());

    	  std::set<int> mAll;
    	  mAll.insert(gm.acode1());
    	  mAll.insert(gm.acode2());

    	  bool a1Pat = pAll.find(a1) != pAll.end();
    	  bool a2Pat = pAll.find(a2) != pAll.end();

    	  bool a1Mat = mAll.find(a1) != mAll.end();
    	  bool a2Mat = mAll.find(a2) != mAll.end();

    	  bool a1Pat_a2Mat = a1Pat && a2Mat;
    	  bool a2Pat_a1Mat = a2Pat && a1Mat;
    	  bool bothAllelesInherited = a1Pat_a2Mat || a2Pat_a1Mat;

    	  if (bothAllelesInherited) {
    		  // Will count this trio's transmission at most once for Het->Ref and at most once for Het->Alt:
    		  bool transRefFromHet = false;
    		  bool transAltFromHet = false;

    		  /*
    		   * NOTE: Only looking at parents that are HETs with REF/ALT and not ALT1/ALT2, since want to assess transmission
    		   * for standard case (where genotypes are more "believable"):
    		   */
    		  if (gp.heterozygote() && (gp.a1IsReference() || gp.a2IsReference())) { // Father is a reference/non-reference het
    			  if (a1Pat_a2Mat) {
    				  if (go.a1IsReference())
    					  transRefFromHet = true;
    				  else
    					  transAltFromHet = true;
    			  }
    			  if (a2Pat_a1Mat) {
    				  if (go.a2IsReference())
    					  transRefFromHet = true;
    				  else
    					  transAltFromHet = true;
    			  }
    		  }

    		  if (gm.heterozygote() && (gm.a1IsReference() || gm.a2IsReference())) { // Mother is a reference/non-reference het
    			  if (a1Pat_a2Mat) {
    				  if (go.a2IsReference())
    					  transRefFromHet = true;
    				  else
    					  transAltFromHet = true;
    			  }
    			  if (a2Pat_a1Mat) {
    				  if (go.a1IsReference())
    					  transRefFromHet = true;
    				  else
    					  transAltFromHet = true;
    			  }
    		  }

    		  if (transRefFromHet)
    			  summary->_trans_ref_from_het++;
    		  if (transAltFromHet)
    			  summary->_trans_alt_from_het++;
    	  }
    	  else {
    		  summary->_nonmendelian++;

    		  /* Both parents are homozygous reference, and one of child alleles is an inherited *reference* allele.
    		   * So, only 1 of 2 child alleles are de novo, which is the case we want to consider (most "believable"):
    		   */
    		  bool bothParentsHomRef = gp.reference() && gm.reference();
    		  bool inheritedOneRef = ((a1Pat || a1Mat) && go.a1IsReference()) || ((a2Pat || a2Mat) && go.a2IsReference());

    		  if (bothParentsHomRef && inheritedOneRef)
    			  summary->_potential_denovo++;
    	  }


      bool denovo = (summary->potential_denovo() > prevSummary.potential_denovo());

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

      bool trans_ref_from_het = (summary->trans_ref_from_het() > prevSummary.trans_ref_from_het());
      bool trans_alt_from_het = (summary->trans_alt_from_het() > prevSummary.trans_alt_from_het());

	  // get allele frequencies
	  int c = 0 , c_tot = 0;
	  if (denovo || (aux->printTransmission && (trans_ref_from_het || trans_alt_from_het)))
		  /* Use instead n_alt_allele() of n_minor_allele(),
		   * since for de novo status and transmission we're specifically tracking:
		   *
		   * a. denovo: If the parents are both hom ref, then does the child have one alt allele
		   * (and we want to know how many total alt alleles *of any kind* are present at this site)
		   *
		   * b. transmission: If a parent is a het (ref+alt), then did the child receive the ref or the alt
		   * (and we want to know how many total alt alleles *of any kind* are present at this site)
		   */
		  v.n_alt_allele( &c , &c_tot );

      // directly output possible de novo events (only REF x REF --> HET)
      // that also passed any above, de-novo specific filters

      if ( denovo )
      {
    	  // track # of actual 'passing' de novo calls
    	  summary->_dcount++;

    	  pDeNovos << "Variant" << "\t"
    			  << "RefxRef->Het[Ref+Alt]" << "\t"
    			  << VAR_DATA << "\n";
      }

      if (aux->printTransmission) {
    	  Out & pTrans = Out::stream( "parent_transmission.vars" );

    	  if (trans_ref_from_het)
    		  pTrans << "Variant" << "\t"
    		  << "Het[Ref+Alt]->Ref_allele" << "\t"
    		  << VAR_DATA << "\n";

    	  if (trans_alt_from_het)
    		  pTrans << "Variant" << "\t"
    		  << "Het[Ref+Alt]->Alt_allele" << "\t"
    		  << VAR_DATA << "\n";
      }
    }
}
}

#define VAR_HEADER \
		"#DATA" << "\t" \
		<< "TRIO_STATUS" << "\t" \
		<< "LOCUS" << "\t" \
		<< "PAT_GT" << "\t" \
		<< "MAT_GT" << "\t" \
		<< "CHILD_GT" << "\t" \
		<< "PAT" << "\t" \
		<< "MAT" << "\t" \
		<< "CHILD" << "\t" \
		<< "AAC" << "\t" \
		<< "AN" << "\t" \
		<< "PAT_META" << "\t" \
		<< "MAT_META" << "\t" \
		<< "CHILD_META" << "\t" \
		<< "CONSENSUS_ALT_ALLELES" << "\t" \
		<< "PAT_ALT_ALLELES" << "\t" \
		<< "MAT_ALT_ALLELES" << "\t" \
		<< "CHILD_ALT_ALLELES"
  
bool Pseq::VarDB::denovo_scan( Mask & mask )
{
	Out & pDeNovos = Out::stream( "denovo.vars" );
	pDeNovos << VAR_HEADER << "\n";

  // did we have some special values

  const int n = g.indmap.size();

  // store summary transmission data
  Aux_transmission_summary aux(n);


  if ( args.has("param") )
    {
      std::vector<double> p = args.as_float_vector( "param" );
      if ( p.size() < 7 )
	Helper::halt( "expect --param DP(kid) DP(par) PL(kid) PL(par) AB(kid,lwr) AB(kid,upr) AB(par,upr) [printTransmission?]" );

      aux.dp_kid = p[0];
      aux.dp_par = p[1];
      aux.pl_kid = p[2];
      aux.pl_par = p[3];

      aux.ab_kid_min = p[4];
      aux.ab_kid_max = p[5];
      aux.ab_par = p[6];

      if (p.size() >= 8)
    	  aux.printTransmission = static_cast<bool>(p[7]);
    }

  Out* outputTrans = NULL;
  if (aux.printTransmission) {
	outputTrans = new Out( "parent_transmission.vars" , "for each parent het, output whether the ref or alt was transmitted (or one from each parent)" );

	Out & pTrans = Out::stream( "parent_transmission.vars" );
	pTrans << VAR_HEADER << "\n";
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

  if (aux.printTransmission)
	  delete outputTrans;

  // display summaries
  Out & pindiv = Out::stream( "denovo.indiv" );
  pindiv << "#DATA\t"
		  << "CHILD" << "\t"
		  << "COMPLETE_TRANSMISSIONS" << "\t"
		  << "MISSING" << "\t"
		  << "TRANS_REF_FROM_HET" << "\t"
		  << "TRANS_ALT_FROM_HET" << "\t"
		  << "ALT_TRANS_RATE" << "\t"
		  << "NON_MENDELIAN" << "\t"
		  << "POTENTIAL_DENOVO" << "\t"
		  << "PASSING_DENOVO" << "\n";

  for (int i=0;i<n;i++)
    {
      aux_trio_transmission_summary * p = aux.indiv(i);
      if ( p->has_parents )
	pindiv << "Individual\t"
	       << g.indmap(i)->id() << "\t"
	       << p->complete_transmissions() << "\t"
	       << p->missing() << "\t"
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
