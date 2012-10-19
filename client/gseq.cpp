#include "views.h"
#include "plinkseq/protdb.h"

#include <map>
#include <set>
#include <vector>

extern GStore g;


void g_geneseq( VariantGroup & vars , void * p )
{
  std::cout << "in here " << vars.name() << "\n";

  Out & pout = Out::stream( "gsview" );
  
  Opt_geneseq * aux = (Opt_geneseq*)p;
  
  Out * Rplot = aux->R_plot ? &Out::stream( "gsview.R" ) : NULL ;   
  
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

      std::cout << "found " << pos << "\n";

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
	      for (int aa= ii->pstart; aa<= ii->pstop; aa++)	    
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
      
      if ( cds && printing && aux->protdb && cpos == 0 && pdm.find( apos ) != pdm.end() )
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
  
  
  pout << "------------------------------------------------------------\n\n";
}

