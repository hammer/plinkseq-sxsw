#include "views.h"
#include "plinkseq/protdb.h"
#include "util.h"

#include <map>
#include <set>
#include <vector>

extern GStore g;
extern Pseq::Util::Options args;

void geneseq_define_Rfunc( Out * );

void g_geneseq( VariantGroup & vars , void * p )
{

  Out & pout = Out::stream( "gsview" );
  
  Opt_geneseq * aux = (Opt_geneseq*)p;
  
  Out * rout = aux->R_plot ? &Out::stream( "gsview.R" ) : NULL ;   
  
  if ( aux->R_plot ) geneseq_define_Rfunc( rout );

  Region region = g.locdb.get_region( PLINKSeq::DEFAULT_LOC_GROUP() , vars.name() ) ;  

  if ( region.subregion.size() == 0 ) return;

  int s = region.subregion[0].meta.get1_int( PLINKSeq::TRANSCRIPT_STRAND() ) ;
  bool positive_strand = s > 0 ;

  if ( s == 0 ) Helper::halt( "no strand info" );

  //
  // get list of all 'events' (variants, intron/exon boundaries, reference variants)
  //
  
  std::set<int> events;  
  std::map<int,const Variant*> evars;
  std::map<int,int> elocstart;
  std::map<int,int> elocstop;
  std::map<int,int> cds_num;
  std::map<int,const RefVariant*> erefvars;
  
  for (int v=0; v<vars.size(); v++)
    {
      int pos = positive_strand ? vars(v).position() : -vars(v).position() ;
      events.insert( vars(v).position() ) ;
      evars[ vars(v).position() ] = &vars(v);
    }
  
  
  //
  // enumerate CDS exons first
  //

  int num_cds_exons = 0;
  int cds_bp = 0;
  int genomic_start = 0;
  int genomic_stop = 0;

  for (int s = 0 ; s < region.subregion.size(); s++) 
    {
      // only consider CDS exons      

      if ( ! region.subregion[s].CDS() ) continue;      

      cds_num[s] = num_cds_exons;
      
      int pos1 = region.subregion[s].start.position();
      int pos2 = region.subregion[s].stop.position();
      
      if ( genomic_start == 0 ) genomic_start = pos1;
      if ( pos2 > genomic_stop ) genomic_stop = pos2;

      // track start & stop CDS (and record exon #)
      elocstart[ pos1 ] = num_cds_exons;
      elocstop[ pos2 ] = num_cds_exons;

      ++num_cds_exons;  
      
      // total CDS extent in base-pairs
      cds_bp += region.subregion[s].stop.position() - region.subregion[s].start.position()  + 1;
    }
  

  //
  // reference variant
  //

  std::set<RefVariant> rvars = g.refdb.lookup( region , aux->ref );
  std::set<RefVariant>::iterator i = rvars.begin();
  while ( i != rvars.end() ) 
    {

      // for now, just place 'start' of reference-variant events
      //int pos = positive_strand ? i->start() : - i->start() ; 
      int pos = i->start();
      events.insert( pos );
      erefvars[ pos ] = &(*i);
      ++i;

      //      std::cout << "found " << pos << "\n";

    }
  
  pout << vars.name() << " | "
       << g.locdb.alias( vars.name() , false ) << " | "
       << num_cds_exons << " CDS exons | " 
       << ( positive_strand ? "+ve strand | " : "-ve strand | " )
       << vars.size() << " variants | ";
  
  if ( rvars.size() > 0 ) 
    pout << rvars.size() << " refvars | ";      

  pout << region.coordinate() << " | " 
       << ( region.stop.position() - region.start.position() + 1 )/1000.0 << " kb | "      
       << cds_bp << " coding bases\n";
  
  //
  // attach LOCDB
  //

  LocDBase * locdb = g.resolve_locgroup( PLINKSeq::DEFAULT_LOC_GROUP() ) ;
  if ( ! locdb ) return;
  
  int gid = locdb->lookup_group_id( PLINKSeq::DEFAULT_LOC_GROUP() );  
  if ( gid == 0 ) return;


  //
  // Translate to AA sequence
  //

  Annotate::setDB( locdb , &g.seqdb );
  
  std::string aa = Annotate::translate_reference( region , false );
  

  //
  // Protein feature/domain annotations
  //
  
  std::map<int,std::string> pdm; //protein domain map

  if ( aux->protdb )
    {
      
      std::set<Feature> features = aux->protdb->fetch( vars.name() );

      bool all_prot = aux->protdom.find( "*" ) != aux->protdom.end() 
	|| aux->protdom.find( "ALL" ) != aux->protdom.end() 
	|| aux->protdom.find( "all" ) != aux->protdom.end(); 

      std::set<Feature>::iterator ii = features.begin();
      while ( ii != features.end() )
	{
	  if ( all_prot || aux->protdom.find( ii->source_id ) != aux->protdom.end() )
	    {
		  // Subtract 1 to have 0-based aa coordinates:
	      for (int aa = ii->pstart - 1; aa <= ii->pstop - 1; aa++)
		{		  
		  if ( pdm[aa] != "" ) pdm[aa] += " ";
		  pdm[ aa ] += ii->source_id + "::" + ii->feature_id + ":" + ii->feature_name + " ";
		}
	    }
          ++ii;	  
	}
    }

      
  
  //
  // Display full AA sequence, with 1-codon padding either side
  //

  int pmin = genomic_start ;
  int pmax = genomic_stop  ;
  
  if ( ! positive_strand )
    {
      int t = pmin;
      pmin = pmax;
      pmax = t;
    }
  

  pmin += positive_strand ? -9 : +9 ;
  pmax += positive_strand ? +9 : -9 ;
  
  int step = positive_strand ? +1 : -1;
  
  
  //
  // codon position, cycle 0,1,2
  // transcript CDS should always start at 0
  //

  int cpos = 0;  
  int gpos = 0;

  // will always start just before CDS, but in 'printable' region, thus:::
  
  bool cds = false;
  bool printing = true;
  
  int exon = 0;
  int exon0 = 0;
  int last_exon = 0;
  bool split_codon = false;

  std::string codon = "";
  std::string prt_codon = "";
  std::string prt_intronic = "";
  int chr = region.start.chromosome();
  
  int apos = 0;
  std::string refannot = "";
  std::string varannot = "";
  bool has_evar = false;
  int stop_here = positive_strand ? pmax+1 : pmax-1;
  
  std::stringstream ss;

  for ( int bp = pmin ; bp != pmax+step ; bp += step )
    {
      
      const int searchbp = bp;
      
      //
      // Are we entering or leaving a printable region?
      // (cpos-adjusted start, to make sure we are always in-sync)
      //
      
      if ( positive_strand )
	{	  
	  if      ( elocstart.find( searchbp + 9 + cpos ) != elocstart.end() ) 
	    {
	      printing = true;

	      // space a new one
	      if ( ! aux->only_variant_sites ) pout << "\n";

	    }
	  else if ( elocstop.find( searchbp - ( 7 + (3 - cpos)  )  ) != elocstop.end() ) printing = false; 	  	  
	}
      else
	{
	  if ( elocstop.find( searchbp - ( 9 + cpos ) ) != elocstop.end() ) 
	    {
	      printing = true;
	      // space a new one
	      if ( ! aux->only_variant_sites ) pout << "\n";
	    }
	  else if ( elocstart.find( searchbp + 7 + (3-cpos) ) != elocstart.end() ) printing = false;
	}
      

      //
      // hitting first position in a new CDS?
      //
      
      if ( ( positive_strand && elocstart.find( searchbp ) != elocstart.end() ) ||
	   ( (!positive_strand) && elocstop.find( searchbp ) != elocstop.end() ) )
	{	  
	  
	  cds = true;	  
	  
	  // track exon number

	  exon = positive_strand ? elocstart[ searchbp ] : elocstop[ searchbp ];
	  exon0 = exon;
	  
	  if ( ! positive_strand ) exon = num_cds_exons - (exon+1) ;	  
	  
	}


      
      //
      // Have we gone one past the end of a current exon?
      //
      
      if ( cds ) 
	if ( ( positive_strand && elocstop.find( searchbp - step ) != elocstop.end() ) 
	     ||
	     ( (!positive_strand) && elocstart.find( searchbp - step ) != elocstart.end() ) 
	     ) 
	  {
	    
	    cds = false;
	    
	    // track last exon code, for use printing introns
	    last_exon = exon+1;

	    if ( positive_strand )
	      {
		if ( elocstop[ searchbp - step ] != exon ) plog.warn("internal inconsistency");
	      }
	    else
	      {
		if ( elocstart[ searchbp - step ] != exon0 ) plog.warn("internal inconsistency");
	      }

	    //
	    // are we ending a CDS mid-frame? 
	    //
	    
	    
	    
	    if ( cpos == 1 ) 
	      {
		split_codon = true;
		//prt_codon = ".";
	      }
	    else if ( cpos == 2 )   
	      {
		split_codon = true;
		//prt_codon = "..";
	      }
	    
	  }
      
      if ( bp == stop_here ) { continue; } 

   
      //
      // append variant information
      //

      if ( printing && evars.find( searchbp ) !=  evars.end() )
	{
	  
	  const Variant * pvar = evars.find( searchbp )->second;

	  has_evar = true;
	  
	  std::stringstream ss2;
	  
	  ss2 << "[" ;
	  
	  int case_count = 0 , case_tot = 0; 
	  int control_count = 0 , control_tot = 0;	  
	  bool ma, control_ma;
	  
	  if ( aux->pheno )       
	    {
	      ma = pvar->n_minor_allele( &case_count , &case_tot , NULL , CASE );
	      control_ma = pvar->n_minor_allele( &control_count , &control_tot , NULL , CONTROL );
	    }
	  else    
	    ma = pvar->n_minor_allele( &case_count , &case_tot );
	  
	  // always report non-ref
	  if ( ! ma ) case_count = case_tot - case_count;
	  if ( ! control_ma ) control_count = control_tot - control_count;
	  
	  // print minor/major allele(s)
	  
	  if ( cpos == 0 ) 
	    ss2 << "" << pvar->alternate() << "..";
	  else if ( cpos == 1 ) 
	    ss2 << "." << pvar->alternate() << ".";
	  else
	    ss2 << ".." << pvar->alternate() ;
	  
	  if ( aux->pheno )
	    ss2 << "|A/U=" << case_count << ":" << control_count ;
	  else
	    ss2 << "|MAC=" << case_count ;
	  
	  ss2 << "|";
	  

	  // always use PSEQ engine for annotation
	  // 	  if ( pvar->meta.has_field( PLINKSeq::META_ANNOT() ) )
	  // 	    ss2 << pvar->meta.get1_string( PLINKSeq::META_ANNOT() );
	  // 	    {
	  // 	  else 
	  
	  bool exonic = Annotate::annotate( (Variant&)*pvar , region );
	      
	  // pull out correct transcript annotation
	  std::vector<std::string> annot = pvar->meta.get_string( PLINKSeq::ANNOT_TYPE() );
	  std::vector<std::string> annot_trans = pvar->meta.get_string( PLINKSeq::ANNOT_GENE() );
	  std::vector<std::string> annot_prot = pvar->meta.get_string( PLINKSeq::ANNOT_PROTEIN() );
	  
	  std::string a = ".";
	  std::string ap = ".";
	  
	  if ( annot.size() == annot_trans.size() && annot.size() == annot_prot.size() )
	    {				
	      for( int i=0;i<annot.size(); i++)
		{		  
		  if ( annot_trans[i] == vars.name() ) 
		    {		      
		      a = annot[i]; ap = annot_prot[i];
		    }
		}
	    }
	  
	  ss2 << a ;
	  if ( exonic ) 
	    {	  		  
	      ss2 << "|" << ap;
	    }
	  
	      
	  
	  ss2 << "]";
	  
	  if ( varannot != "" ) varannot += " ";
	  varannot += ss2.str();
	  
	}
  

      //
      // append ref-variant information
      //

      if ( printing && erefvars.find( searchbp ) !=  erefvars.end() )
	{
	  const RefVariant * refvar = erefvars.find( searchbp )->second;
	  std::stringstream ss2;
	  ss2 << *refvar;
	  if ( refannot != "" ) refannot += " ";
	  refannot += ss2.str();
	}
      
      // cpos == 2 <==> last base in codon (since still on 0-based base count within exon: 0,1,2), so only annotate once
      // [Can change to 'cpos == 0' if want to annotate the first part of the 'split codon' instead of the last part]:
      if ( printing && aux->protdb && cds && cpos == 2 && pdm.find( apos ) != pdm.end() )
	{
	  refannot += pdm[ apos ];
	}

      
      //
      // Only show CDS
      //
      

      if ( printing )
	{	  
	  if ( gpos == 0 ) 
	    if ( cds )
	      ss << "exon " << exon+1 << ( exon < 9 ? " " : "" ) << "     " << Helper::chrCode( chr ) << ":" << searchbp;  	  
	    else 
	      {	      
		if ( last_exon == 0 ) 
		  ss << "intron */1 " << ( exon < 9 ? " " : "" );
		else if ( exon == num_cds_exons - 1 ) 
		  ss << "intron " << last_exon << "/* " << ( exon < 9 ? " " : "" );
		else 
		  ss << "intron " << last_exon << "/" << last_exon+1 << " " << ( exon < 9 ? " " : "" );
		
		ss << Helper::chrCode( chr ) << ":" << searchbp;
		
	      }
	}
      
      
      if ( printing )  
	{
	  
	  ++gpos;

	  if ( cds )
	    {
	      codon        += g.seqdb.lookup( chr , bp ) ; 
	      prt_codon    += g.seqdb.lookup( chr , bp ) ; 
	      prt_intronic += ".";
	      ++cpos;
	    }
	  else
	    {
	      prt_codon    += ".";	  
	      prt_intronic += g.seqdb.lookup( chr , bp ) ; 	      
	    }
	  

	  //
	  // Print row
	  //
	  
	  if ( gpos == 3 )
	    {
	      
	      gpos = 0;	      
	      
	      std::string aa_code = cds && cpos == 3 ? aa.substr( apos++ , 1 ) : ( split_codon ? ">" : "." ) ;
	      std::string aa_name = cds && cpos == 3 ? Annotate::aa[ aa_code ] : ( split_codon ? ">>>" : " . " );
	      
	      ss << "\t" << prt_intronic << " " << prt_codon
		 << " " << aa_code << " " << aa_name << " " << ( cds ? Helper::int2str(apos) : ( split_codon ? "> " : ". " ) ) 
		 << "\t" << (refannot == "" ? "." : refannot ) 
		 << "\t" << (varannot == "" ? "." : varannot ) 		       
		 << "\n";
	      	      
	      // print? 
	      if ( (!aux->only_variant_sites) || has_evar ) 
		pout << ss.str();
	      
	      // reset all codon-specific stuff
	      if ( ! split_codon ) codon = "";
	      prt_codon = "";
	      prt_intronic = "";
	      refannot = "";
	      varannot = "";	      
	      has_evar = false;
	      split_codon = false;
	      
	      // and clear stream
	      ss.str( std::string() );
	      
	    }
	}
      
      if ( gpos == 3 ) gpos = 0;
      if ( cpos == 3 ) cpos = 0;

    }




  //
  // Create R plot?
  //

  if ( aux->R_plot ) 
    {

      //
      //  Upfront stuff 
      //

      *rout << "## Populate data-structures for " << vars.name() << "\n"
	    << "\n"
	    << "transname = \"" << vars.name() << "\"\n"
	    << "genename = \"" << g.locdb.alias( vars.name() , false ) << "\"\n" 
	    << "chrcode = \"" << Helper::chrCode( vars(0).chromosome() ) << "\"\n"
	    << "refname = \"" << ( aux->ref ? args.as_string( "ref" ) : "" ) << "\"\n"
	    << "var <- list();ref <- list();dom <- list();exon  <- list() \n"

	    << "# main transcript \n"

	    << "trans <<- c( " << region.start.position() << " , " << region.stop.position() << ") \n" 

	    << "# strand for main (and only) transcript (+1, -1, 0) \n"

	    << " strand <<- " << ( positive_strand ? "+1" : "-1" ) << "\n" 
	
	    << "# determine border \n"	
	    << " translen <<- trans[2] - trans[1] + 1 \n"
	    << " total <<- round(c( trans[1] - 1000 , trans[2] + 1000 ) ) \n"
	    << " totallen <<- total[2] - total[1] + 1  \n"

	    << " ## exon structure \n";
      
      

      //
      // look at all CDS exons
      //

      int exc = 1;
      for (int s = 0 ; s < region.subregion.size(); s++) 
	{	  
	  if ( ! region.subregion[s].CDS() ) continue;      	  
	  *rout << " exon[[" << exc << "]] <- c( " << region.subregion[s].start.position() << " , " << region.subregion[s].stop.position() << " ) \n";
	  ++exc;
	}
      
      *rout << " cdslength <<- " << cds_bp << " \n";
      

      //
      // non-CDS exons
      //

      exc = 1;
      *rout << "exon_notcds <- list() \n";
      for (int s = 0 ; s < region.subregion.size(); s++) 
	{	  
	  if ( ! region.subregion[s].exon() ) continue;      	  
	  *rout << " exon_notcds[[" << exc << "]] <- c( " << region.subregion[s].start.position() << " , " << region.subregion[s].stop.position() << " ) \n";
	  ++exc;
	}
      

      //
      // All other overlapping transcripts
      //

      std::string g_chr_code = Helper::chrCode( region.start.chromosome() );
      int g_chr = region.start.chromosome();
      int g_bp1 = region.start.position();
      int g_bp2 = region.stop.position();

      std::set<Region> others = g.locdb.get_regions( g.locdb.lookup_group_id( PLINKSeq::DEFAULT_LOC_GROUP() ) , g_chr , g_bp1 , g_bp2 );      



      *rout << "others <- list() \n";

      if ( others.size() > 0 ) 
	{
	  int cnt1 = 1;
	  std::set<Region>::iterator ii = others.begin();
	  while ( ii != others.end() )
	    {
	      // for this transcript, exclude self
	      if ( ii->subregion.size() > 0 && ii->name != vars.name() )
		{		  
		  *rout << "others[[" << cnt1 << "]] <- list( name = \"" << ii->name << "\" , \n ";
		  
		  *rout << "exon_cds = list( " ;

		  int cnt2 = 1;
		  for (int s = 0 ; s < ii->subregion.size(); s++) 
		    {	  
		      if ( ii->subregion[s].CDS() ) { *rout << (cnt2>1?",":"") << " c( " << ii->subregion[s].start.position() << " , " << ii->subregion[s].stop.position() << " ) \n"; ++cnt2; }
		    }		  
		  
		  *rout << " ) , \n";
		  
		  *rout << "exon_notcds = list( " ;
		  
		  cnt2 = 1;
		  for (int s = 0 ; s < ii->subregion.size(); s++) 
		    if ( ii->subregion[s].exon() ) { *rout << (cnt2>1?",":"") << " c( " << ii->subregion[s].start.position() << " , " << ii->subregion[s].stop.position() << " ) \n"; ++cnt2; }
		  *rout << " ) \n";
		  
		  *rout << " ) \n";

		  ++cnt1;
		}
	      
	      ++ii;
	    }
	}
      

      //
      // Variants;
      //

      *rout << " ## variants, with MAC (by cases), with annotation \n";
      
      for (int v=0;v<vars.size();v++)
	{
	  
	  std::string annot = ".";
	  std::string annotdet = ".";
	  int cnta , cntu;

	  int case_count = 0 , case_tot = 0; 
	  int control_count = 0 , control_tot = 0;	  
	  bool ma, control_ma;
	  
	  Variant * pvar = &vars(v);

	  if ( aux->pheno )       
	    {
	      ma = pvar->n_minor_allele( &case_count , &case_tot , NULL , CASE );
	      control_ma = pvar->n_minor_allele( &control_count , &control_tot , NULL , CONTROL );
	    }
	  else    
	    ma = pvar->n_minor_allele( &case_count , &case_tot );
	  
	  // always report non-ref
	  if ( ! ma ) case_count = case_tot - case_count;
	  if ( ! control_ma ) control_count = control_tot - control_count;
	  
	  if ( ! aux->pheno ) 
	    {
	      cntu = -ma;
	    }
	  else
	    {
	      cnta = case_count;
	      cntu = control_count;
	    }
	  	  
	  bool exonic = Annotate::annotate( (Variant&)*pvar , region );

	  // pull out correct transcript annotation
	  std::vector<std::string> annot_func = pvar->meta.get_string( PLINKSeq::ANNOT_TYPE() );
	  std::vector<std::string> annot_trans = pvar->meta.get_string( PLINKSeq::ANNOT_GENE() );
	  std::vector<std::string> annot_prot = pvar->meta.get_string( PLINKSeq::ANNOT_PROTEIN() );	  
	  if ( annot_func.size() == annot_trans.size() && annot_func.size() == annot_prot.size() )
	    {				
	      for( int i=0;i<annot_func.size(); i++)
		{		  
		  if ( annot_trans[i] == vars.name() ) 
		    {		      
		      annot = annot_func[i]; 
		      annotdet = annot_prot[i];
		    }
		}
	    }


	  *rout << "var[[" << v+1 << "]] <- list( pos = c( " << vars(v).position() << " , " <<  vars(v).stop() 
		<< " ) , name = \"" << vars(v).name() 
		<< "\" , annot = \"" << annot << "\" , annotdet = \"" << annotdet 
		<< "\" , a = " << cnta << " , u = " << cntu << " ) \n";

	  
	}

      
      //
      // Reference variants
      //
      
      if ( aux->ref ) 
	{
	  *rout << "## Reference variants\n";

	  *rout << "ref[[1]] <- list( group = \"" << args.as_string( "ref" ) << "\" , det = list() ) \n";
	  
	  std::set<RefVariant>::iterator i = rvars.begin();
	  int rfc = 1;
	  while ( i != rvars.end() ) 
	    {      
	      int pos = i->start();
	      *rout << "ref[[1]]$det[[" << rfc++ << "]] <- list( pos = c( " << i->start() << " , " << i->stop() << " ) , name = \"" << i->name() << "\" ) \n";
	      ++i;
	    }
	}


      //
      // Protein domains
      //
      
      if ( aux->protdb )
	{

	  // just repeat this lookup from above for now...

	  *rout << "## Protein domains \n";
	  
	  std::set<Feature> features = aux->protdb->fetch( vars.name() );
	  
	  bool all_prot = aux->protdom.find( "*" ) != aux->protdom.end() 
	    || aux->protdom.find( "ALL" ) != aux->protdom.end() 
	    || aux->protdom.find( "all" ) != aux->protdom.end(); 
	  
	  int pdct = 1;
	  
	  std::map<std::string,int> sourcemap;
	  std::map<std::string,int> sourcecnt;
	  
	  std::set<Feature>::iterator ii = features.begin();
	  while ( ii != features.end() )
	    {
	      if ( all_prot || aux->protdom.find( ii->source_id ) != aux->protdom.end() )
		{

		  if ( sourcemap.find( ii->source_id ) == sourcemap.end() )
		    {
		      int t = sourcemap.size() + 1 ;
		      sourcemap[ ii->source_id ] = t;
		      *rout << "dom[[" << t << "]] <- list( group = \"" << ii->source_id << "\" , det = list() ) \n";    
		    }
		  

		  // 1-based count
		  		  
		  *rout << "dom[[" << sourcemap[ ii->source_id ] << "]]$det[[" << ++sourcecnt[ ii->source_id ] 
			<< "]] <- list( name = \"" << ii->feature_id << ":" << ii->feature_name << "\" , "
			<< "pos = c( " << ii->gstart << " , " << ii->gstop << " ) , aa = c( " << ii->pstart << " , " << ii->pstop << " ) ) \n ";
		}
	      ++ii;	  
	    }
	}
      

  // Perform actual plot, redircted to PDF
  
  *rout << "\n\n#create plot\n"
	<< "pdf( width = 14 , height = 10 , file=\"plot" << vars.name() << ".pdf\") \n "
	<< "doplot() \n"
	<< "dev.off() \n\n"
	<< "#---------------------------------------------------------------------------\n\n\n";
 
}



  
  
  pout << "------------------------------------------------------------\n\n";
}



void geneseq_define_Rfunc( Out * rout )
{

  *rout << "## Map of gaps\n"
	<< "gapmap <- function( x , e , st = 10 ) "
	<< "{ e <- c( list( c(x[1],x[1]) ) , e , list( c( x[2],x[2] ) ) )\n"
//	<< " t  <- (x[2] - x[1] + 1) * 0.005 \n"
	<< " t  <- 20 \n"   // simple 20 base rule for intron-skipping
	<< " g <- list() \n"
	<< " for (i in 2:length(e)) \n"
	<< " { if ( e[[i]][1] - e[[i-1]][2] > t ) g <- c( g , list( c( e[[i-1]][2] + st , e[[i]][1] - st  ) ) ) } \n"
	<< " return(g) \n"
	<< "} \n\n";
  
  *rout << "## Mapping f()\n"
	<< " f <- function(x) { \n"
	<< " if ( length(m)==0) return(x) \n"
	<< " if ( strand == +1 ) { \n "
	<< " if ( x < total[1] | x > total[2] ) return(0) \n"
	<< " p <- 0 \n"
	<< " for (i in 1:length(m)) \n"
	<< " {\n"
	<< " if ( x < m[[i]][1] ) return( x - p - total[1] + 1 ) \n"
	<< " if ( x >= m[[i]][1] & x <= m[[i]][2] ) return(0) \n"
	<< " if ( x > m[[i]][2] ) p <- p + ( m[[i]][2] - m[[i]][1] + 1 ) \n"
	<< " } \n"
	<< " return( x - p - total[1] + 1 ) \n"
	<< " } \n"

    // as above, but for the negative strand, flip things
	<< " else { \n"
	<< " if ( x < total[1] | x > total[2] ) return(0) \n"
	<< " if ( x == total[2] ) return(1) \n"   // kludge
	<< " p <- 0 \n"
	<< " for (i in length(m):1) \n"
	<< " {\n"
	<< " if ( x > m[[i]][2] ) return( total[2] - x - p + 1 ) \n"
	<< " if ( x >= m[[i]][1] & x <= m[[i]][2] ) return(0) \n"
	<< " if ( x < m[[i]][1] ) p <- p + ( m[[i]][2] - m[[i]][1] + 1 ) \n"
	<< " } \n"
	<< " return( total[2] - x - p + 1 ) \n"
	<< " } \n"
    // end of f()
	<< " } \n"
    // make a vectorized version
	<< " vf <- Vectorize(f) \n";

  *rout << "# is a given single point in the map range?\n"
	<< "inmap <- function( x ) \n"
	<< " { for (i in 1:length(m)) { if ( m[[i]][1] <= x & m[[i]][2] >= x ) return(F) } \n"
	<< " return(T) \n"
	<< "} \n\n";

  *rout << " genomic <- function(x) \n"
	<< " { return( tmin + ( ( x - total[1] ) / ( total[2] - total[1] ) ) * tmax ) } \n"
	<< " vgenomic <- Vectorize(genomic) \n\n";


  *rout << "doplot <- function() { \n"

	<< "par(lend=2,ljoin=3, oma=c(0,0,0,0) , mar=c(1,1,1,1) ) \n"

	<< "## Make reduced transcript mapping -- used by f()\n"
	<< "m <<- gapmap( total , exon ) \n"
	<< "## Reverse mapping (i.e. to plot genomic scale on actual grid) \n"
	<< " tmin <<- f(total[ifelse( strand == +1 , 1 , 2 ) ]) \n"
	<< " tmax <<- f(total[ifelse( strand == +1 , 2 , 1 ) ]) \n"
	
// ## vertical tracks:

// ## 9 chromosome/BP lines                1
// ## 8 intronic ref-variants              length(ref)
// ## 7 intronic-variants                  1
// ## 6 genomic-scale                      1
// ## 5  lines                             1
// ## 4 transcript-scale (introns skipped) 1
// ## 3 protein-domains                    length(dom)
// ## 2 exonic variants                    2 (allow more space for encoding)
// ## 1 exonic ref-variants                length(ref)

	<< " tracks = 8 + round(length(others)/8) + length(ref)*2 + length(dom) \n"
    
	<< " tr_exref = 1 \n"
	<< " tr_exvar = 1 + length(ref) \n"   
	<< " tr_dom   = tr_exvar + 2 \n"      // i.e. exvar is two deep
	<< " tr_trans = tr_exvar + 1 + length(dom) + 1 \n"
	<< " tr_lines = tr_trans + 0 \n"
	<< " tr_genom = tr_lines + 1 \n"
	<< " tr_alt   = tr_genom + 1 \n" // alternative transcripts
	<< " tr_invar = tr_alt   + round(length(others)/8) + 1 \n"
	<< " tr_inref = tr_invar + 1 \n"
	<< " tr_scale = tr_invar + length(ref) + 1 \n\n"

	<< "## Main box \n"
    
	<< "plot( c( f( total[1] ) , f( total[2] ) ) , c( 0.5 , tracks+0.5 ) , ylab=\"\" , xlab=\"\" , yaxt= \"n\" , xaxt=\"n\", type = \"n\" )  \n";


  // draw in 'gaps' ; in domains; in trans in exvar
  
  *rout << "## indicate where gaps are in coding sequence \n"
	<< " for( i in 1:length(m)) \n"
	<< " { \n"
 	<< "   lines( c( f( m[[i]][1]-1 ) , f( m[[i]][1]-1 ) ) , c( tr_exvar , tr_trans + 1 ) , col=\"gray\" , lty=2 , lwd= 1 ) \n "
	<< " } \n\n";
  
  
  // draw CDS and other exons for main transcript

   *rout << "## Exons \n"

	 << " text( tmin , tr_trans + 1 , labels= ifelse( strand == +1 , \"+ve coding\" , \"-ve coding\" ) , pos=3,offset=0.2,col=\"gray\" , cex=0.6 ) \n";
   
     
     // CDS
   *rout << " if ( length(exon)>0 ) \n"
	 << " { lines( vf( c( exon[[1]][1] , exon[[length(exon)]][2] ) ) , c( tr_trans+0.5 , tr_trans+0.5 ) ,col=\"blue\",lwd=1,lend=2 ) } \n"
	 << "  for (i in 1:length(exon)) { \n "
	 << " yp <- ifelse( i %% 2 == 1 , tr_trans+ 0.5 , tr_trans + 0.2 ) \n"
	 << " p1 <- f(exon[[i]][1]) \n"
	 << " p2 <- f(exon[[i]][2]) \n"
     
	 << "  rect( p1 , yp , p2 , yp+0.3, col=\"lightblue\",border=\"blue\") \n"
     
     // exon number
	 << " text( min(p1,p2)+10 , ifelse( i %% 2 == 1 , tr_trans+ 0.25 , tr_trans + 0.55 ) , labels = ifelse( strand == +1 , i , length(exon)-i+1 ) , pos=3, offset=0, cex=0.6 , col = \"black\" ) \n "  
     
     // now transcript always on positive strand 
	 << " arrows( min(p1,p2) , yp+0.15 , min(p1,p2)+12 , yp+0.15 , length=0.05 , col=\"blue\" , ljoin=1 , lend=1) \n"

	 << " } \n\n";
     
     

     //
     // legend at top (genomic sequence scale)
     //

     *rout << "## Sequence info at top \n"
       
	   << "lines( c(tmin,tmin),c( tr_scale + 0.2, tr_scale + 0.4),col=\"black\") \n"
	   << "lines( c(tmax,tmax),c( tr_scale + 0.2, tr_scale + 0.4),col=\"black\") \n"
	   << "msg <- paste( chrcode , \":\" , round( total[1]/1e6 , 2 ),\"Mb\",sep=\"\")  \n" 
	   << "msg2 <- paste( genename , \";\", transname , \"; total transcript \", round(translen/1e3),\"kb; coding \" , cdslength, \"bp; \", ifelse( strand == +1 , \"+ve strand\" , \"-ve strand\" ) , sep=\"\") \n"
	   << "text( tmin , tr_scale + 0.4 , labels = msg , offset=0.2 , pos=3,cex=0.6) \n"
	   << "text( (tmax-tmin)/2 , tr_scale + 0.5 , labels = msg2 , offset=0 , pos=3,cex=0.7) \n"
	   << "text( tmax , tr_scale +0.4 , labels = paste( round( total[2]/1e6 , 2 ),\"Mb\" , sep=\"\" ) , offset=0.2 , pos=3,cex=0.6) \n"
       
	   << "# ticks every 1/10kb \n"
	   << "kb1 <- seq( total[1] , total[2] , length.out = round(translen/1e3) ) \n"
	   << "kb10 <- seq( total[1] , total[2] , length.out = round(translen/1e3)/10 ) \n"
	   << "for (i in kb1) { lines( c(genomic(i),genomic(i)) , c(tr_scale + 0.25,tr_scale + 0.35) , col=\"gray\" , lend=2) } \n"
	   << "for (i in kb10) { lines( c(genomic(i),genomic(i)) , c( tr_scale + 0.2, tr_scale + 0.4) , col=\"gray\" , lend=2) } \n"
       

       //
       // transcripts on genomic scale
       //

	   << "## Transcript on genomic-scale \n"

       // also show non-CDS exons here
	   << " if ( length(exon_notcds)>0 ) \n"
	   << " { lines( vf( c( exon_notcds[[1]][1] , exon_notcds[[length(exon_notcds)]][2] ) ) , c( tr_genom+0.5 , tr_genom+0.5 ) ,col=\"gray\",lwd=1,lend=2 )  \n"
	   << "  for (i in 1:length(exon_notcds)) { \n "	   
	   << " p1 <- genomic(exon_notcds[[i]][1]) \n"
	   << " p2 <- genomic(exon_notcds[[i]][2]) \n"
	   << "  rect( p1 , tr_genom+0.6 , p2 , tr_genom+1 , col=\"gray\",border=\"gray\") \n"
	   << " } }\n"


	   << " text( tmin , tr_genom + 1 , labels=\"+ve genomic\" , pos=3,offset=0.2,col=\"gray\" , cex=0.6 ) \n"
	   << " lines( c( genomic( trans[1] ) , genomic( trans[2] ) ) , c( tr_genom + 0.8 ,  tr_genom + 0.8 ) , col=\"gray\" , lwd=1 , lend=2)  \n"

     // and CDS

	   << " for (i in 1:length(exon)) \n"
	   << " { \n"
	   << " yp <- tr_genom \n"
	   << " xx <-  ( f(exon[[i]][1]) + f(exon[[i]][2]) ) / 2.0 \n"
	   << " x2 <-  ( genomic( exon[[i]][1] )  + genomic( exon[[i]][2] ) ) / 2.0 \n"
       
	   << " lines( c( xx , x2 ) , c(yp,yp+.6) , col=\"gray\" , lty=1,lwd=.5,lend=2)  \n"
       
	   << " if ( i %% 2 == 0 )  \n"
	   << " {   lines( c(x2,x2) , c(yp+.6,yp+.8) , col=\"gray\" , lty=1,lwd=.5,lend=2)  \n"
	   << "      lines( c(xx,xx) , c(yp,yp-.5) , col=\"gray\" , lty=1,lwd=.5,lend=2)   } \n"
	   << " else \n"
	   << "  { lines( c(x2,x2) , c(yp+.6,yp+.8) , col=\"gray\" , lty=1,lwd=.5,lend=2)  \n"
	   << "    lines( c(xx,xx) , c(yp,yp-.2) , col=\"gray\" , lty=1,lwd=.5,lend=2)  \n"
	   << " } \n"
       
	   << "  offset <- ifelse( i %% 2 == 1 , -0.2 , 0 )  \n"
	   << " rect( genomic( exon[[i]][1] ) ,  tr_genom + 0.8 + offset ,  \n"
	   << " genomic( exon[[i]][2] ) , tr_genom + 0.8 + offset + 0.2 , col=\"lightblue\" , border=\"blue\" , lend=2) \n"
	   << " } \n\n";


       //
       // Other transcripts on genomic scale
       //
     
     *rout << "## other transcripts on genomic scale\n"
	   << "if(length(others)>0) { \n"
	   << " yp <- tr_alt+0.4 \n "
	   << " for (t in 1:length(others)) { \n"
	   << "  text( tmin , yp+0.06 , labels= others[[t]]$name , cex=0.3 ) \n"  // labels
	   << "  lines( c( genomic( others[[t]]$exon_notcds[[1]][1]) , genomic( others[[t]]$exon_notcds[[length(others[[t]]$exon_notcds)]][2]) ) , c( yp,yp ) , lwd=0.25 ) \n "
	   << "  for (a1 in 1:length(others[[t]]$exon_notcds)) { \n"
	   << "    lines( c( genomic( others[[t]]$exon_notcds[[a1]][1]) , genomic(others[[t]]$exon_notcds[[a1]][2] ) ) , c( yp,yp ) , lwd=2 , lend=2) \n "
	   << "  } \n"
 	   << "  if ( length(others[[t]]$exon_cds) > 0 ) for (a1 in 1:length(others[[t]]$exon_cds)) { \n"
 	   << "    lines( c( genomic( others[[t]]$exon_cds[[a1]][1]) , genomic(others[[t]]$exon_cds[[a1]][2] ) ) , c( yp,yp ) , lwd=4 , lend=2) \n "
 	   << "  } \n"
	   << " yp <- yp + 0.15 \n "
	   << " } \n"
	   << " }\n\n";
	   


       //
       // Domains
       //

	 *rout << "## Domains \n"
	 
	       << " domcol <- c(\"darkgreen\",\"brown\",\"orange\",\"purple\",\"darkred\") \n"
	       << " domcol2 <- c(\"lightgreen\",\"beige\",\"yellow\",\"pink\",\"red\") \n"
	 
	       << " if ( length(dom) > 0 ) { \n"
	       << " for (i in 1:length(dom)) \n"
	       << " { \n"
	       << "  yp <- tr_dom + i - 1 \n "
	       << "  for (j in 1:length(dom[[i]]$det)) \n"
	       << " { \n"
	       << " ypact <- ifelse( j %% 2 == 1 , yp + 0.20 , yp + 0.3 ) \n"
	       << "  rect( f( dom[[i]]$det[[j]]$pos[1]) ,ypact ,  \n"
	       << "        f( dom[[i]]$det[[j]]$pos[2]) ,ypact + 0.2 ,   \n"
	       << "   col=domcol2[i %% 5 +1 ] , border= domcol[ i %% 5 +1 ] ) \n"
	 
	       << " ypact <- ifelse( j %% 2 == 1 , yp + 0.05 , yp + 0.6 ) \n"
	       << " text( min( vf( c( dom[[i]]$det[[j]]$pos[1] , dom[[i]]$det[[j]]$pos[2] ) ) ) , ypact , cex=.6, labels= dom[[i]]$det[[j]]$name , pos=4,offset=0)  \n"
	       << " } \n"
	       << "  text( 1 , yp + 0.6 , labels=dom[[i]]$group , cex=.6,pos=4) \n"
	       << " } } \n"
	 

       //
       // all reference variants
       //

	<< "## Reference variants \n"

	 << " if ( length(ref) > 0 ) { \n"
	 << " text( tmin , tr_exref-0.2  , labels= refname , pos=3,offset=0.2,col=\"black\" , cex=0.6 ) \n"     
	 << " text( tmin , tr_inref+0.5  , labels= refname , pos=3,offset=0.2,col=\"black\" , cex=0.6 ) \n"     
	 << " for (i in 1:length(ref)) \n"
	 << " { \n" 
	
	<< " ypin <- tr_inref + i - 1 \n"
	<< " ypex <- tr_exref + i - 1 \n"

	<< "  if ( length(ref[[i]]$det) > 0 ) { \n" 
	<< " for (j in 1:length(ref[[i]]$det)) \n"
	<< "  { \n" 
	<< "  ex1 <- inmap( ref[[i]]$det[[j]]$pos[1] ) \n"
	<< "  ex2 <- inmap( ref[[i]]$det[[j]]$pos[2] ) \n"

	<< "  if ( ex1 & ex2 ) \n"
	<< "  { \n"
	<< "    lines( vf( ref[[i]]$det[[j]]$pos ) , c( ypex - 0.2 , ypex + 0.0) , col=\"darkgreen\" ) \n" 
    //	<< "  text( f( ref[[i]]$det[[j]]$pos[1] ) + 5 , ypex  , col=\"darkgray\" , pos=1 , labels= ref[[i]]$det[[j]]$name , cex=0.6 , srt=270) \n"
	<< " } \n"
	<< " lines( vgenomic( ref[[i]]$det[[j]]$pos ) , c( ypin + 0.5, ypin + 0.7) , col=\"darkgreen\" ) \n"
	<< "  } \n"
	<< " }}} \n"


     
     // all actual, observed variant data

	<< "## Variants \n"

	<< " lines( vgenomic(trans) , c( tr_invar  + 0.5, tr_invar + 0.5) , col=\"gray\" , lwd=1) \n"

     // note ; kludge below, not the f() function misses the edge case for -ve strand genes, so just shift back 1 base
	<< " lines( c( f(trans[1]) , f(trans[2]-1) ) , c( tr_exvar + 0.5, tr_exvar + 0.5) , col=\"gray\" , lwd=1) \n"
     
     //<< "cat( \" vf trans \" , f(trans[1]-1) , f(trans[2]-1) , \"\n\" ) \n" 
     
// # silent   = gray
// # missense = black
// # nonsense = green

	<< " for (i in 1:length(var)) \n"
	<< " { \n" 
	 
     //	 << " cat( \"in var\" , i , \"\n\" ) \n "
     
	 << " p1 <- var[[i]]$pos[1] \n"
	 << " p2 <- var[[i]]$pos[2] \n"
	 << " n  <- var[[i]]$name \n"
	 << " a  <- var[[i]]$a \n"
	 << " u  <- var[[i]]$u \n"
	 << " cc <- u >= 0 \n"
	 << " ex1 <- inmap( p1 ) \n"
	 << " ex2 <- inmap( p2 )  \n"
	 << " cl <- \"gray\" \n"
	 << "  if ( var[[i]]$annot == \"missense\" ) cl <- \"black\" \n"
	 << "  if ( var[[i]]$annot == \"nonsense\" ) cl <- \"red\" \n"
     
	 << " if ( ex1 & ex2 ) \n"
	 << " { \n"
	 << "   rect( f(p1)-1 , tr_exvar + 0.45 , f(p2)+1 , tr_exvar + 0.55 , col=\"white\" , border=cl , lend=2) \n"

	 << " if ( cc ) {  \n"
	 << " text( f( p1 )  , tr_exvar + 0.7 , col=\"red\" , labels=a , cex=0.6)  \n"
	 << " text( f( p1 )  , tr_exvar + 0.3 , col=\"blue\" , labels=u , cex=0.6 ) \n"
     
     // circle plot for CC counts
 	 << " if ( a+u < 20 ) {  \n"
 	 << "  yt <- tr_exvar + 0.9 \n"
 	 << "  if ( u>0 ) { for (j in 1:u ) { \n "
	 << " points( f( p1 ) , yt , col = \"blue\" , pch=21 , bg = \"white\" , cex=0.6 ) \n"
         << " yt <- yt + 0.1 \n"
	 << "} }\n "
 	 << " if ( a > 0 ) { for (j in 1:a ) { \n"
	 << " points( f( p1 ) , yt , col = \"red\" , pch=21 , bg = \"white\" , cex=0.6 ) \n"
	 << " yt <- yt + 0.1 \n"
	 << " } } \n "
 	 << " }  \n"
     
     //	 << " else {   }  \n "

     // if not C/C data
	 << " } \n"
	 << " else { text( f( p1 )  , tr_exvar + 0.2 , col=\"gray\" , labels=-u , cex=0.6 ) } \n"

     
	 << " if ( var[[i]]$annotdet != \".\" ) text( f(p1) + 5 , tr_exvar - 0.3 , col=\"black\" , labels=var[[i]]$annotdet , cex=0.5 , pos = 1 , offset = 0 , srt=270 )  \n"
	 << " } \n" 

     // also show in genomic space
	<< " rect( genomic(p1) , tr_invar + 0.4 , genomic(p2) , tr_invar + 0.6 , col=cl , border=cl , lend=2)  \n"

     // do not try to show counts for genomic variants
//  	<< " if ( cc ) { \n"
//  	<< " text( genomic( p1 )  , tr_invar + 0.7 , col=\"red\" , labels=a , cex=0.6)  \n"
//  	<< " text( genomic( p1 ) , tr_invar + 0.3 , col=\"blue\" , labels=u , cex=0.6 ) \n"
//  	<< " } \n"
//  	<< " else { text( genomic( p1 ) , tr_invar + 0.2 , col=\"gray\" , labels=-u , cex=0.6 ) } \n"
     
     // genomic circle plot for CC counts
 	 << " if ( cc & a+u < 20 ) {  \n"
 	 << "  yt <- tr_invar + 0.7 \n"
 	 << "  if ( u>0 ) { for (j in 1:u ) { \n "
	 << " points( genomic( p1 ) , yt , col = \"blue\" , pch=21 , bg = \"white\" , cex=0.6 ) \n"
         << " yt <- yt + 0.1 \n"
	 << "} }\n "
 	 << " if ( a > 0 ) { for (j in 1:a ) { \n"
	 << " points( genomic( p1 ) , yt , col = \"red\" , pch=21 , bg = \"white\" , cex=0.6 ) \n"
	 << " yt <- yt + 0.1 \n"
	 << " } } \n "
 	 << " }  \n"


     // next variant
	<< " } \n";


  // end of do plot

  *rout << "}\n\n";


}






