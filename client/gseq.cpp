#include "views.h"
#include "plinkseq/protdb.h"

#include <map>
#include <set>
#include <vector>

extern GStore g;


void g_geneseq( VariantGroup & vars , void * p )
{

  Out & pout = Out::stream( "gsview" );

  Opt_geneseq * aux = (Opt_geneseq*)p;
  
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
  
  

  // exons
  int num_cds_exons = 0;
  int cds_bp = 0;
  for (int s =0 ; s < region.subregion.size(); s++) 
    {
      // only consider CDS exons      
      if ( ! region.subregion[s].CDS() ) continue;      

      cds_num[s] = num_cds_exons;
//       int pos1 = positive_strand ? region.subregion[s].start.position() : -region.subregion[s].start.position(); 
//       int pos2 = positive_strand ? region.subregion[s].stop.position() : -region.subregion[s].stop.position(); 

      int pos1 = region.subregion[s].start.position();
      int pos2 = region.subregion[s].stop.position();

      elocstart[ pos1 ] = num_cds_exons;
      elocstop[ pos2 ] = num_cds_exons;

      ++num_cds_exons;  
      
      cds_bp += region.subregion[s].stop.position() - region.subregion[s].start.position()  + 1;
    }
  
  // if -ve strand, flip numbers;
  

  // reference variant
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
    }
  
  pout << vars.name() << " | "
       << num_cds_exons << " CDS exons | " 
       << ( positive_strand ? "+ve strand | " : "-ve strand | " )
       << vars.size() << " variants | ";
  if ( rvars.size() > 0 ) 
    pout << rvars.size() << " refvars | ";      
  pout << ( region.stop.position() - region.start.position() + 1 )/1000.0 << " kb | "      
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
      
      std::set<Feature>::iterator ii = features.begin();
      while ( ii != features.end() )
	{
	  if ( ii->source_id == "TMHMM" )
	    {
	      for (int aa= ii->pstart; aa<= ii->pstop; aa++)	    
		{		  
		  //		  if ( pdm[aa] != "" ) pdm[aa] += " ";
		  pdm[ aa ] += ii->source_id + "::" + ii->feature_id + ":" + ii->feature_name + " ";
		}
	    }
          ++ii;	  
	}
    }

      
  
  //
  // Display full AA sequence
  //

  int pmin = region.start.position();
  int pmax = region.stop.position();
  
  if ( ! positive_strand )
    {
      int t = pmin;
      pmin = pmax;
      pmax = t;
    }

  int step = positive_strand ? +1 : -1;
  
  // codon position, cycle 0,1,2
  // transcript CDS should always start at 0
  int cpos = 0;  
  
  bool cds = false;
  
  if ( ( positive_strand && elocstart.find( pmin ) != elocstart.end() ) ||
       ( (!positive_strand) && elocstop.find( pmin ) != elocstop.end() ) )
    cds = true;
  
  
  int exon = 0;
  
  std::string codon = "";
  std::string prt_codon = "";
  int chr = region.start.chromosome();
  
  int apos = 0;
  std::string refannot = "";
  std::string varannot = "";
  bool print_pos = true;
  int stop_here = positive_strand ? pmax+1 : pmax-1;
  bool hide = false;
  
  for ( int bp = pmin ; bp != pmax+step ; bp += step )
    {
      
      //      const int searchbp = positive_strand ? bp : -bp;
      
      const int searchbp = bp;
      
      //      std::cout << "\n---------------------------------\nconsidering " << searchbp << "\n";

      //
      // in CDS?
      //
      
      if ( ( positive_strand && elocstart.find( searchbp ) != elocstart.end() ) ||
	   ( (!positive_strand) && elocstop.find( searchbp ) != elocstop.end() ) )
	{

	  cds = true;

	  exon = positive_strand ? elocstart[ searchbp ] : elocstop[ searchbp ];
	  
	  if ( ! positive_strand ) exon = num_cds_exons - (exon+1) ;
	  
	  //	  pout << "---- start exon " << exon+1 << ( positive_strand ? "(+)" : "(-)" ) << " ---- \n";
	  pout << "\n";

	  print_pos = true;
	}

      if ( ( positive_strand && elocstop.find( searchbp - step ) != elocstop.end() ) 
	   ||
	   ( (!positive_strand) && elocstart.find( searchbp - step ) != elocstart.end() ) 
	   ) // look 1-past end
	{
	  
	  cds = false;

	  if ( elocstop[ searchbp - step ] != exon ) plog.warn("hmm");
	  
	  // are we mid-frame? 
	  if ( cpos == 1 ) 
	    {
	      if ( ! hide ) pout << "\t" << codon << "..\n" ; 
	      prt_codon = ".";
	    }
	  else if ( cpos == 2 ) 
	    {
	      if ( ! hide ) pout << "\t" << codon << ".\n" ; 
	      prt_codon = "..";
	    }
	  
	  // if ( ! hide ) pout << "----   end exon " << exon+1 << "    ---- \n";
	}
      
      if ( bp == stop_here ) { continue; } 
      

   

//       if ( print_pos ) 
// 	{
// 	  std::set<int>::iterator ne = events.upper_bound( bp );
// 	  bool ohide = hide;
// 	  hide = true;
// 	  if ( ne != events.end() && *ne - bp <= ( cds ? 3 : 2 ) ) 
// 	    {
// 	      hide = false;
// 	      //std::cout << "hide true = " << bp << " " << *ne << " (1)\n";
// 	    }
// 	  if ( ne != events.begin() )
// 	    {
// 	      --ne;
// 	      if ( bp - *ne <= ( cds? 3 : 2 )  ) hide = false;
// 	      //std::cout << "hide true = " << bp << " " << *ne << " (2)\n";
// 	    }
	  
// // 	  if ( hide && ! ohide ) // going into hiding...
// // 	    pout << "  ...\n";	  

// 	}

//       // Over-ride 'hide' feature 

//       hide = false;




      //
      // append variant information
      //

      if ( evars.find( searchbp ) !=  evars.end() )
	{
	  const Variant * pvar = evars.find( searchbp )->second;
	  std::stringstream ss;
	  
	  ss << "[" ;
	  
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
	    ss << "" << pvar->alternate() << "..";
	  else if ( cpos == 1 ) 
	    ss << "." << pvar->alternate() << ".";
	  else
	    ss << ".." << pvar->alternate() ;
	  
	  if ( aux->pheno )
	    ss << "|A/U=" << case_count << ":" << control_count ;
	  else
	    ss << "|MAC=" << case_count ;
	  
	  ss << "|";
	  
	  if ( ! pvar->meta.has_field( PLINKSeq::META_ANNOT() ) )
	    {
	      bool exonic = Annotate::annotate( (Variant&)*pvar , &region );
	      ss << pvar->meta.get1_string( PLINKSeq::ANNOT_TYPE() );
	      if ( exonic ) 
		{	  
		  //ss << pvar->meta.get1_string( PLINKSeq::ANNOT_GENE() );
		  ss << "|" << pvar->meta.get1_string( PLINKSeq::ANNOT_PROTEIN() );
		}
	    }
	  else
	    ss << pvar->meta.get1_string( PLINKSeq::META_ANNOT() );
	  
	  ss << "]";

	  // CURRENTLY, ONLY APPEND CDS VARIANTS:

	  if ( cds )
	    {
	      
	      if ( varannot != "" ) varannot += " ";
	      varannot += ss.str();
	    }
	}
  
      
      // append ref-variant information
      
      if ( cds && erefvars.find( searchbp ) !=  erefvars.end() )
	{
	  const RefVariant * refvar = erefvars.find( searchbp )->second;
	  std::stringstream ss;
	  ss << *refvar;
	  if ( refannot != "" ) refannot += " ";
	  refannot += ss.str();
	}

      if ( cds && aux->protdb && cpos == 0 && pdm.find( apos ) != pdm.end() )
	{
	  refannot += pdm[ apos ];
	}

      
      //
      // Only show CDS
      //

      if ( print_pos && cds ) 
	{	  
	  if ( ! hide ) pout << ( cds ? "exon " + Helper::int2str(exon+1) + " " : "       " ) + "chr" << chr << ":" << searchbp;  
	  print_pos = false;
	}
      
      

      // display in coding region
      //      std::cout << "\n -- " << bp << " " << cds << " " << print_pos << "\n";

      if ( cds )  
	{
	  
	  //std::cout << "\nadding " << cpos << "\t" << bp << " " << g.seqdb.lookup( chr , bp ) << "\n";
	  
	  codon += g.seqdb.lookup( chr , bp ) ; 
	  prt_codon += g.seqdb.lookup( chr , bp ) ; 
	  
	  ++cpos;
	  
	  if ( cpos == 3 ) 
	    {

	      cpos = 0;
	      std::string aa_code = aa.substr( apos++ , 1 ) ;	     
	      std::string aa_name = Annotate::aa[ aa_code ];

	      if ( ! hide ) 
		{
		  pout << "\t" << prt_codon << " " << aa_code << " " << aa_name << " " <<apos ;
		  
		  pout << "\t" << (refannot == "" ? "." : refannot ) 
		       << "\t" << (varannot == "" ? "." : varannot ) 		       
		       << "\n";
		}

	      codon = "";
	      prt_codon = "";
	      refannot = "";
	      varannot = "";
	      print_pos = true;
	    }
	}

      if ( false && ! cds )
	{
	  if ( ! hide )
	    {
	      pout << "\t" << g.seqdb.lookup( chr , bp ) 
		   << "   . .";
	      
	      pout << "\t" << (varannot == "" ? "." : varannot ) 
		   << "\t" << (refannot == "" ? "." : refannot ) 
		   << "\n";
	    }
	  
	  refannot = "";
	  varannot = "";
	  print_pos = true;	  
	}
      
    }

//   std::map<int,const Variant*> evars;
//   std::map<int,int> elocstart;
//   std::map<int,int> elocstop;
//   std::map<int,const RefVariant*> erefvars;

  // foooter
  pout << "\n------------------------------------------------------------\n\n";
}

