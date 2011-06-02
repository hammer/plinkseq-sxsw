#include "views.h"

#include <map>
#include <set>
#include <vector>

extern GStore g;

void g_geneseq( VariantGroup & vars , void * p )
{
  Opt_geneseq * aux = (Opt_geneseq*)p;

  // assume refseq for now... 
  Region region = g.locdb.get_region( "refseq" , vars.name() ) ;  
  if ( region.subregion.size() == 0 ) return;
  bool positive_strand = region.subregion[0].meta.get1_string( PLINKSeq::TRANSCRIPT_STRAND() ) != "-";
  
  
  // get list of all 'events' (variants, intron/exon boundaries, reference variants)
  
  std::set<int> events;  
  std::map<int,const Variant*> evars;
  std::map<int,int> elocstart;
  std::map<int,int> elocstop;
  std::map<int,const RefVariant*> erefvars;
  
  for (int v=0; v<vars.size(); v++)
    {
      int pos = positive_strand ? vars(v).position() : -vars(v).position() ;
      events.insert( pos ) ;
      evars[ pos ] = &vars(v);
    }


  // exons
  int cds_bp = 0;
  for (int s =0 ; s < region.subregion.size(); s++) 
    {
      int pos1 = positive_strand ? region.subregion[s].start.position() : -region.subregion[s].start.position(); 
      int pos2 = positive_strand ? region.subregion[s].stop.position() : -region.subregion[s].stop.position(); 
//      events.insert( pos1 );
      elocstart[ pos1 ] = s;
//      events.insert( pos2 );
      elocstop[ pos2 ] = s;
      cds_bp += region.subregion[s].stop.position() - region.subregion[s].start.position()  + 1;
    }
  
  // reference variant
  std::set<RefVariant> rvars = g.refdb.lookup( region , aux->ref );
  std::set<RefVariant>::iterator i = rvars.begin();
  while ( i != rvars.end() ) 
    {
      // for now, just place 'start' of reference-variant events
      int pos = positive_strand ? i->start() : - i->start() ; 
      events.insert( pos );
      erefvars[ pos ] = &(*i);
      ++i;
    }
  
  plog << vars.name() << " | "
       << region.subregion.size() << " exons | "
       << ( positive_strand ? "+ve strand | " : "-ve strand | " )
       << vars.size() << " variants | ";
  if ( rvars.size() > 0 ) 
    plog << rvars.size() << " refvars | ";      
  plog << ( region.stop.position() - region.start.position() + 1 )/1000.0 << " kb | "      
       << cds_bp << " coding bases\n";
  
  
  // attach SEQDB

  LocDBase * locdb = g.resolve_locgroup( "refseq" ) ;
  if ( ! locdb ) return;
  
  int gid = locdb->lookup_group_id( "refseq" );  
  if ( gid == 0 ) return;
  
  Annotate::setDB( locdb , &g.seqdb );
  
  std::string aa = Annotate::translate_reference( region , false );
  
  // display in 3-base 
  int pmin = *events.begin();
  int pmax = *(--events.end());

  // Assume reading frame starts for gene either at pmin (+ve strand) 
  // or at pmax (-ve strand) and that the full  
  
  // start of gene is always in frame;
  // but if on =ve strand, this means pmax

  if ( ! positive_strand ) 
    {
      int t = pmin;
      pmin = pmax;
      pmax = t;
    }

  int step = positive_strand ? +1 : -1;

  int cpos = 0;  // codon position, cycle 0,1,2
  bool cds = true;
  int exon = 0;
  std::string codon = "";
  std::string prt_codon = "";
  int chr = region.start.chromosome();
  int apos = 0;
  std::string refannot = "";
  std::string varannot = "";
  bool print_pos = true;
  int stop_here = positive_strand ? pmax+1 : pmin-1;
  bool hide = false;

  
  for ( int bp = pmin ; bp <= (pmax+1) ; bp += step )
    {
      const int searchbp = positive_strand ? bp : -bp;

      // are we a long way to the next interesting event ? if so, skip
      
      if ( print_pos ) 
	{
	  std::set<int>::iterator ne = events.upper_bound( bp );
	  bool ohide = hide;
	  hide = true;
	  if ( ne != events.end() && *ne - bp <= ( cds ? 3 : 2 ) ) 
	    {
	      hide = false;
	      //std::cout << "hide true = " << bp << " " << *ne << " (1)\n";
	    }
	  if ( ne != events.begin() )
	    {
	      --ne;
	      if ( bp - *ne <= ( cds? 3 : 2 )  ) hide = false;
	      //std::cout << "hide true = " << bp << " " << *ne << " (2)\n";
	    }
	  
	  if ( hide && ! ohide ) // going into hiding...
	    plog << "  ...\n";	  
	}


      // in CDS?
      if ( elocstart.find( searchbp ) != elocstart.end() ) 
	{
	  cds = true;
	  exon = elocstart[ searchbp ];
	  if ( ! hide ) 
	    plog << "---- start exon " << exon+1 << ( positive_strand ? "(+)" : "(-)" ) << " ---- \n";
	  print_pos = true;
	}
      
      if ( elocstop.find( searchbp - step ) != elocstop.end() ) // look 1-past end
	{
	  cds = false;
	  if ( elocstop[ searchbp -step ] != exon ) plog.warn("hmm");
	  
	  // are we mid-frame? 
	  if ( cpos == 1 ) 
	    {
	      if ( ! hide ) plog << "\t" << codon << "..\n" ; 
	      prt_codon = ".";
	    }
	  else if ( cpos == 2 ) 
	    {
	      if ( ! hide ) plog << "\t" << codon << ".\n" ; 
	      prt_codon = "..";
	    }
	  
	  if ( ! hide ) plog << "----   end exon " << exon+1 << "    ---- \n";
	}
      
      if ( bp == stop_here ) continue;
      
      // append variant information
      
      if ( evars.find( searchbp ) !=  evars.end() )
	{
	  const Variant * pvar = evars.find( searchbp )->second;
	  std::stringstream ss;

	  ss << pvar->position();
	  
	  int case_count = 0 , case_tot = 0; 
	  int control_count = 0 , control_tot = 0;	  
	  bool ma, control_ma;
	  
	  if ( aux->pheno )       
	    {
	      ma = pvar->n_minor_allele( case_count , case_tot , CASE );
	      control_ma = pvar->n_minor_allele( control_count , control_tot , CONTROL );
	    }
	  else    
	    ma = pvar->n_minor_allele( case_count , case_tot );
	  
	  // always report non-ref
	  if ( ! ma ) case_count = case_tot - case_count;
	  if ( ! control_ma ) control_count = control_tot - control_count;
	  
	  // print minor/major allele(s)
	  
	  if ( cpos == 0 ) 
	    ss << "|" << pvar->alternate() << "..";
	  else if ( cpos == 1 ) 
	    ss << "|." << pvar->alternate() << ".";
	  else
	    ss << "|.." << pvar->reference() ;
	  
	  if ( aux->pheno )
	    ss << "|" << case_count << ":" << control_count ;
	  else
	    ss << "|" << case_count ;
	  
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
	
	  if ( varannot != "" ) varannot += " ";
	  varannot += ss.str();
	}
  
      
      // append ref-variant information
      
      if ( erefvars.find( searchbp ) !=  erefvars.end() )
	{
	  const RefVariant * refvar = erefvars.find( searchbp )->second;
	  std::stringstream ss;
	  ss << *refvar;
	  if ( refannot != "" ) refannot += " ";
	  refannot += ss.str();
	}

      if ( print_pos ) 
	{	  
	  if ( ! hide ) plog << ( cds ? "exon "+Helper::int2str(exon+1) + " " : "       " ) + "chr" << chr << ":" << bp;	      
	  print_pos = false;
	}
      
      
      // display in coding region
      
      if ( cds )  
	{
	  	  
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
		  plog << "\t" << prt_codon << " " << aa_code << " " << aa_name;
		  
		  plog << "\t" << (varannot == "" ? "." : varannot ) 
		       << "\t" << (refannot == "" ? "." : refannot ) 
		       << "\n";
		}

	      codon = "";
	      prt_codon = "";
	      refannot = "";
	      varannot = "";
	      print_pos = true;
	    }
	}

      if ( ! cds )
	{
	  if ( ! hide )
	    {
	      plog << "\t" << g.seqdb.lookup( chr , bp ) 
		   << "   . .";
	      
	      plog << "\t" << (varannot == "" ? "." : varannot ) 
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
  plog << "\n------------------------------------------------------------\n\n";
}

