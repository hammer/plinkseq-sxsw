#include "plinkseq/annot.h"
#include "plinkseq/variant.h"
#include "plinkseq/gstore.h"
#include "plinkseq/filemap.h"
#include "plinkseq/output.h"

#include <algorithm>
#include <cmath>

extern GStore * GP;

using namespace std;
using namespace Helper;

std::map<uint64_t,Region> Annotate::rmap;
uint64_t Annotate::transcript_group_id = 0;

std::map<seq_annot_t,std::string> populate_seqinfo()
{
  std::map<seq_annot_t,std::string> m;
  m[UNDEF]    		= ".";
  m[MONO]     		= "monomorphic";
  m[IGR]     	        = "intergenic-region";
  m[INTRON]  	        = "intronic";
  m[UTR5]     		= "5-UTR";
  m[UTR3]     		= "3-UTR";
  m[SYN]      		= "silent";
  m[MIS]     		= "missense";
  m[INDEL] 	  	= "indel";
  m[CODONDELETION]      = "codon-deletion";
  m[STOPDELETION]       = "stop-deletion";
  m[OOFCODONINSERTION]  = "out-of-frame-codon-insertion";
  m[STOPINSERTION]      = "stop-insertion";
  m[CODONINSERTION]     = "codon-insertion";
  m[OOFCODONDELETION]   = "out-of-frame-codon-deletion";
  m[NON]      		= "nonsense";
  m[SL]	      		= "start-lost";
  m[PART]     		= "partial-codon";
  m[DONORIN2] 		= "splice-donor-in2";
  m[DONOREX2AG]         = "splice-donor-ex2ag";
  m[ACCEPTOREX1G]       = "splice-acceptor-ex1g";
  m[ACCEPTORIN2]        = "splice-acceptor-in2";
  m[DONORIN45AG]        = "splice-donor-in45ag";
  m[SPLICE]   		= "splice";
  m[EXONIC_UNKNOWN]     = "exonic-unknown";
  m[FRAMESHIFT]         = "frameshift";
  m[RT]                 = "read-through";
  m[SPLICEDEL]          = "splice-deletion";
  m[npcRNA]             = "npcRNA";
  return m;
}

std::map<std::string,std::string> populate_t()
{
  // stop codon is '*' 
  std::map<std::string,std::string> m;
  m["TAA"] = "*"; m["TAG"] = "*"; m["TGA"] = "*"; m["GCT"] = "A"; m["GCC"] = "A";
  m["GCA"] = "A"; m["GCG"] = "A"; m["TGT"] = "C"; m["TGC"] = "C"; m["GAT"] = "D";
  m["GAC"] = "D"; m["GAA"] = "E"; m["GAG"] = "E"; m["TTT"] = "F"; m["TTC"] = "F";
  m["GGT"] = "G"; m["GGC"] = "G"; m["GGA"] = "G"; m["GGG"] = "G"; m["CAT"] = "H";
  m["CAC"] = "H"; m["ATT"] = "I"; m["ATC"] = "I"; m["ATA"] = "I"; m["AAA"] = "K";
  m["AAG"] = "K"; m["TTG"] = "L"; m["TTA"] = "L"; m["CTT"] = "L"; m["CTC"] = "L";
  m["CTA"] = "L"; m["CTG"] = "L"; m["ATG"] = "M"; m["AAT"] = "N"; m["AAC"] = "N";
  m["CCT"] = "P"; m["CCC"] = "P"; m["CCA"] = "P"; m["CCG"] = "P"; m["CAA"] = "Q";
  m["CAG"] = "Q"; m["CGT"] = "R"; m["CGC"] = "R"; m["CGA"] = "R"; m["CGG"] = "R";
  m["AGA"] = "R"; m["AGG"] = "R"; m["TCT"] = "S"; m["TCC"] = "S"; m["TCA"] = "S";
  m["TCG"] = "S"; m["AGT"] = "S"; m["AGC"] = "S"; m["ACT"] = "T"; m["ACC"] = "T";
  m["ACA"] = "T"; m["ACG"] = "T"; m["GTT"] = "V"; m["GTC"] = "V"; m["GTA"] = "V";
  m["GTG"] = "V"; m["TGG"] = "W"; m["TAT"] = "Y"; m["TAC"] = "Y";
  return m;
}

std::map<std::string,std::string> populate_aa()
{
  // stop codon is '*'
  std::map<std::string,std::string> m;
  m["A"] ="Ala"; m["L"] ="Leu"; m["R"] ="Arg"; m["K"] ="Lys"; m["N"] ="Asn";
  m["M"] ="Met"; m["D"] ="Asp"; m["F"] ="Phe"; m["C"] ="Cys"; m["P"] ="Pro";
  m["Q"] ="Gln"; m["S"] ="Ser"; m["E"] ="Glu"; m["T"] ="Thr"; m["G"] ="Gly";
  m["W"] ="Trp"; m["H"] ="His"; m["Y"] ="Tyr"; m["I"] ="Ile"; m["V"] ="Val";
  m["*"] ="***";
  return m;
}

std::map<seq_annot_t,std::string> SeqInfo::types = populate_seqinfo();
std::map<std::string,std::string> Annotate::t = populate_t();
map<std::string,std::string> Annotate::aa = populate_aa();


SeqDBase * Annotate::seqdb = NULL;
LocDBase * Annotate::db = NULL;

void Annotate::setDB( const fType t )
{
  // by default, try point to LOCDB
  if ( t == SEGDB ) db = &(GP->segdb);
  else db = &(GP->locdb);
}

void Annotate::setDB( LocDBase * p , SeqDBase * s )
{
  db = p;
  seqdb = s;
}

bool Annotate::set_transcript_group( const std::string & grp )
{
  init();
  setDB( LOCDB );
  if ( ( ! db ) || ( ! db->attached() ) ) return false;
  uint64_t id = db->lookup_group_id( grp );
  if ( id == 0 ) return false;
  transcript_group_id = id;
  return true;
}

Region * Annotate::pointer_to_region( const std::string & name )
{
  if ( ! db ) return NULL;
  uint64_t loc_id = db->get_region_id( transcript_group_id , name );
  return from_cache( loc_id );
}

Region * Annotate::from_cache( uint64_t id )
{
  std::map<uint64_t,Region>::iterator ii = rmap.find( id );
  if ( ii != rmap.end() ) return &ii->second;
  std::vector<uint64_t> t(1,id);
  add_transcripts( t );
  // now it should exist
  return &rmap[ id ];
}

void Annotate::add_transcripts( const std::vector<uint64_t> & id )
{

  // stop cache from expanding indefinitely
  if ( rmap.size() > 100000 ) rmap.clear();
  
  for (int i=0;i<id.size();i++)
    {
      // is this transcript in the cache already?
      if ( rmap.find( id[i] ) != rmap.end() ) continue;
      
      // if not, pull from LOCDB
      if ( ( ! db ) || ( ! db->attached() ) || id[i] == 0 ) 
	Helper::halt( "no LOCDB attached, or group not found in LOCDB, so cannot annotate" );
      
      // add actual region to cache
      rmap[ id[i] ] = db->get_region( id[i] );
    }
}


bool Annotate::annotate(Variant & var , const Region & region )
{  
  std::set<Region> t;  
  t.insert( region );
  return annotate( var , t );
}


bool Annotate::annotate(Variant & var , const std::set<Region> & pregion )
{
  // if giving a list of predefined regions, assume these have valud 'id' that
  // corresponds to the LOCDB/group specified, and add to cache
  
  std::vector<uint64_t> ids;
  std::set<Region>::const_iterator ii = pregion.begin();
  while ( ii != pregion.end() )
    {
      ids.push_back( ii->id );
      rmap[ ii->id ] = *ii;
      ++ii;
    }
  
  // and now call standard Annotate function
  
  return annotate( var , ids );
}


bool Annotate::annotate(Variant & var)
{
  // lookup IDs from LOCDB
  std::vector<uint64_t> ids = db->get_region_ids( transcript_group_id , var.chromosome() , var.position() , var.stop() );
  return annotate( var , ids );
}


bool Annotate::annotate(Variant & var , const std::vector<uint64_t> & ids )
{
  
  std::set<SeqInfo> s;

  // do the actual annotation, for all transcripts indicated in 'ids')
  s = annotate( var.chromosome() ,
		var.position() ,
		var.alternate() ,
		var.reference() , 
		ids ) ;
  
  bool annot = false;

  // summary of 'worst' annotation (int, syn, nonsyn)

  int is_silent 		= 0;
  int is_missense 		= 0;
  int is_codondeletion 	        = 0;
  int is_codoninsertion         = 0;
  int is_frameshift 	        = 0;
  int is_indel 			= 0;
  int is_splice 		= 0;
  int is_esplice 		= 0;
  int is_csplice 		= 0;
  int is_nonsense 		= 0;
  int is_startlost		= 0;
  int is_readthrough 	        = 0;
  int is_intergenic 	        = 0;
  int is_intronic 		= 0;
  int is_exonic_unknown         = 0;
  int is_splicedel              = 0;
  int is_utr3                    = 0;
  int is_utr5                    = 0;
  int is_npcRNA                 = 0;

  std::set<SeqInfo>::iterator i = s.begin();
  while ( i != s.end() )
    {
      // track whether this is coding, for a 'single' return code
      if ( i->synon() ) ++is_silent;
      if ( i->frameshift() ) {++is_frameshift; ++is_indel;}
      if ( i->codondeletion() ) {++is_codondeletion; ++is_indel;}
      if ( i->codoninsertion() ) {++is_codoninsertion; ++is_indel;}
      if ( i->missense() ) ++is_missense;
      if ( i->nonsense() ) ++is_nonsense;
      if ( i->startlost() ) ++is_startlost;
      if ( i->splice() ) ++is_splice;
      if ( i->esplice() ) ++is_esplice;
      if ( i->csplice() ) ++is_csplice;
      if ( i->readthrough() ) ++is_readthrough;
      if ( i->intergenic() ) ++is_intergenic;
      if ( i->intronic() ) ++is_intronic;
      if ( i->indel() ) ++is_indel;
      if ( i->exonic_unknown() ) ++is_exonic_unknown;
      if ( i->npcRNA() ) ++is_npcRNA;
      if ( i->utr3() ) ++is_utr3;
      if ( i->utr5() ) ++is_utr5;

      // add annotations
      var.meta.add( PLINKSeq::ANNOT_TYPE() , i->status() );
      var.meta.add( PLINKSeq::ANNOT_GENE() , i->gene_name() );
      var.meta.add( PLINKSeq::ANNOT_CODING() , (int)i->coding() );
      
      // var.meta.add( PLINKSeq::ANNOT_EXONIC() , i->exonic() );

      var.meta.add( PLINKSeq::ANNOT_CHANGE() , i->genomic() );
      var.meta.add( PLINKSeq::ANNOT_CODON() , i->codon() );
      var.meta.add( PLINKSeq::ANNOT_PROTEIN() , i->protein() );
      
      // add any details...
      i->details( var );

      ++i;
    }
  
  //
  // what is the 'worst' annotation?
  //

  std::string aworst = "";
  if ( is_frameshift ) aworst = "frameshift";
  else if ( is_nonsense ) aworst = "nonsense";
  else if ( is_startlost ) aworst = "start-lost";
  else if ( is_esplice ) aworst = "esplice";
  else if ( is_splicedel ) aworst = "esplice";
  else if ( is_csplice ) aworst = "csplice";
  else if ( is_readthrough ) aworst = "readthrough";
  else if ( is_codoninsertion ) aworst = "codon-insertion";
  else if ( is_codondeletion ) aworst = "codon-deletion";
  else if ( is_missense ) aworst = "missense";
  else if ( is_indel ) aworst = "indel";
  else if ( is_splice ) aworst = "splice";
  else if ( is_exonic_unknown ) aworst = "exonic-unknown";
  else if ( is_silent ) aworst = "silent";
  else if ( is_utr5 ) aworst = "UTR5";
  else if ( is_utr3 ) aworst = "UTR3";
  else if ( is_npcRNA ) aworst = "npcRNA";
  else if ( is_intronic ) aworst = "intronic";
  else aworst = "intergenic";

  var.meta.set( PLINKSeq::ANNOT() , aworst );

  //
  // what is the 'summary/consensus' annotation?
  //

  int acount = 0;
  if ( is_silent ) ++acount;
  if ( is_missense ) ++acount;
  if ( is_nonsense ) ++acount;
  if ( is_startlost ) ++acount;
  if ( is_esplice ) ++acount;
  if ( is_csplice ) ++acount;
  if ( is_splice ) ++acount;
  if ( is_readthrough ) ++acount;
  if ( is_intergenic ) ++acount;
  if ( is_frameshift ) ++acount;
  if ( is_codoninsertion ) ++acount;
  if ( is_codondeletion ) ++acount;
  if ( is_exonic_unknown ) ++acount;
  if ( is_splicedel ) ++acount;
  if ( is_utr3 ) ++acount;
  if ( is_utr5 ) ++acount;
  if ( is_npcRNA ) ++acount;

  std::string annot_summary = aworst;
  if ( acount > 1 ) annot_summary = "mixed";

  annot_summary += ",NON=" + Helper::int2str( is_nonsense );
  annot_summary += ",SL="  + Helper::int2str( is_startlost );
  annot_summary += ",MIS=" + Helper::int2str( is_missense );
  annot_summary += ",SYN=" + Helper::int2str( is_silent );
  annot_summary += ",ESPL=" + Helper::int2str( is_esplice );
  annot_summary += ",CSPL=" + Helper::int2str( is_csplice );
  annot_summary += ",SPL=" + Helper::int2str( is_splice );
  annot_summary += ",RTH=" + Helper::int2str( is_readthrough );
  annot_summary += ",INT=" + Helper::int2str( is_intronic );
  annot_summary += ",IGR=" + Helper::int2str( is_intergenic );
  annot_summary += ",FRAMESHIFT="  + Helper::int2str( is_frameshift );
  annot_summary += ",INDEL="  + Helper::int2str( is_indel );
  annot_summary += ",CODONINSERTION="  + Helper::int2str( is_codoninsertion );
  annot_summary += ",CODONDELETION="  + Helper::int2str( is_codondeletion );
  annot_summary += ",UNKNOWN="  + Helper::int2str( is_exonic_unknown );
  annot_summary += ",SPLICEDEL=" + Helper::int2str( is_splicedel );
  annot_summary += ",npcRNA=" + Helper::int2str( is_npcRNA );
  annot_summary += ",UTR3=" + Helper::int2str( is_utr3 );
  annot_summary += ",UTR5=" + Helper::int2str( is_utr5 );
  var.meta.set( PLINKSeq::ANNOT_SUMMARY() , annot_summary );

  //
  // did we receive any annotation?
  //

  return is_silent 
    || is_startlost 
    || is_nonsense 
    || is_missense 
    || is_splice 
    || is_esplice 
    || is_csplice 
    || is_readthrough 
    || is_frameshift 
    || is_codondeletion 
    || is_indel 
    || is_codoninsertion
    || is_exonic_unknown
    || is_splicedel
    || is_utr3
    || is_utr5
    || is_npcRNA;

}

std::set<SeqInfo> Annotate::annotate( int chr,
				      int bp1 ,
				      const std::string & alternate ,
				      const std::string & reference ,
				      const std::vector<uint64_t> & pregions )
{

  //
  // Store annotations generated in here
  //

  std::set<SeqInfo> annot;

  //
  // Check SEQDB is attached
  //

  if ( ( ! seqdb ) || ( ! seqdb->attached() ) ) 
    return annot;
  
  
  // Pull regions to use for annotation from database, if they weren't supplied directly

  //
  // Variant does not fall into any transcript 
  //

  if ( pregions.size() == 0 )
    {
      annot.insert( SeqInfo( IGR ) );
      return annot;
    }

 
  //
  // Consider each alternate allele, one at a time
  //
    
  std::set<std::string> alt = Helper::parseCommaList( alternate );
  std::set<std::string>::iterator a = alt.begin();
  while ( a != alt.end() )
    {

      std::string act_ref = reference;
      std::string act_alt = *a;
      int act_bp1 = bp1;
      int act_bp2 = act_alt.size() > act_ref.size() ? 
        act_bp1 + act_alt.size() - 1 :
        act_bp1 + act_ref.size() - 1 ;
      int bp2 = act_bp2;


      //
      // Trim actual REF/ALT range to only reflect variant bases
      //

      int min_size = act_ref.size() < act_alt.size() ? act_ref.size() : act_alt.size();
      int o1 = 0 , o2 = 0;

      for (int i=0;i<min_size;i++)
        {
          if ( act_ref[i] == act_alt[i] ) ++o1;
          else break;
        }
      if ( min_size > 1 && act_ref.size() == act_alt.size() ) 
        {
          for (int i=min_size-1;i>=0;i--)
            {
              if ( act_ref[i] == act_alt[i] ) ++o2;
              else break;
            }
        }


      // need to trim?
      if ( o1 || o2 ) 
        {
          act_bp1 += o1;
          act_bp2 -= o2;

          if ( act_bp1 > act_bp2 ) 
            {
              ++a;
              plog.warn("found invalid allele");
              continue;
            }
          
          if ( o1 > act_ref.size() ) act_ref = ""; else act_ref = act_ref.substr( o1 , act_ref.size() - o2 - o1 ) ; 
          if ( o1 > act_alt.size() ) act_alt = ""; else act_alt = act_alt.substr( o1 , act_alt.size() - o2 - o1 ) ; 
        }

      //
      // Consider each transcript supplied, that should overlap this position
      //
      
      std::vector<uint64_t>::const_iterator ii = pregions.begin();
      while ( ii != pregions.end() )
	{
	  
	  //
	  // Get actual region, from rmap
	  //
	  
	  Region * r = Annotate::from_cache( *ii );
	
	  //
	  // We assume the encoding of CDS/exons/start/stop following the LOCDB convention; this 
	  // will automatically be the case if a GTF file was loaded into the LOCDB, so this is
	  // effectively a requirement to use annotation functions (i.e. that the LOCDB was populated 
	  // with GTFs)
	  //
	  
	  if ( r->subregion.size() == 0 )
	    {
	      ++ii;
	      continue;
	    }

	  //
	  // Distinguish here between CDS and non-CDS exons
	  //

	  Region r_exon; // contains all exons (i.e., will include 3UTR , 5UTR, but not 'stop')
	  Region r_cds;  // only contains CDS (and extra stop-codon)
	  
	  for ( int ss=0;ss < r->subregion.size();ss++)
	    {

	      if ( r->subregion[ss].CDS() || r->subregion[ss].stop_codon() )
		r_cds.subregion.push_back( r->subregion[ss] );
	      
	      else if ( r->subregion[ss].exon() )
		r_exon.subregion.push_back( r->subregion[ss] );
	      
	    }

	  //
	  // Which exon(s) does this mutation impact?  Pull in
	  // neighbouring exon if needed. Assume all subregions are on
	  // the same chromosome.
	  //
	  std::set<int> CDS_exons;
	  int pos = 1;
	  
	  //
	  // Verify whether or not splice variant is in frame
	  //

	  int notinframe = 0;
	  
	  //
	  // Strand and exon status
	  //
	  
	  bool negative_strand = false;
	  bool positive_strand = false; // to check at least one strand is given, test both

	  for (int ss=0;ss< r_cds.subregion.size(); ss++)
	    {
	      int s = r_cds.subregion[ss].meta.get1_int( PLINKSeq::TRANSCRIPT_STRAND() );
	      if ( s == 0 ) continue;
	      negative_strand = s < 0 ;
	      positive_strand = s > 0 ;
	      break;
	    }

	  //
	  // skip if both strands are represented 	  
	  //

	  if ( negative_strand && positive_strand )
	    {
	      ++ii;
	      continue;
	    }

	  int first_exon = negative_strand ? 0 : r_exon.subregion.size()-1;
	  int last_exon = negative_strand ? r_exon.subregion.size()-1 : 0;

	  //
	  // Does variant fall within an exon, or near an intron/exon
	  // splice-site boundary size of exons before exon skipping
	  // ... Here we assume exon skipping for splice variants. We
	  // want to distinguish between splice variants that lead to
	  // NMD and those that escape NMD but maintain proper
	  // translation frame.
	  //

	  int transtruncsize = 0;
	  
	  //
	  // want to keep track of the size of the exon we are interested in
	  //

	  int sizeexonint = 0;

	  //
	  // Keep track of the size of the transcript where
	  // final splicing occurs. Cases where this will not work:
	  // Exon Number = 1. Have not investigated NMD in genes with
	  // exon size = 1 -- under investigation.
	  //

	  int sizepenult = 0;
	  int inCDS = -1;
	  int inExon = -1;
	  for ( unsigned int s  = 0 ; s < r_cds.subregion.size(); s++)
	    {
	      //
	      // let us keep all the exons to compute protein
	      // truncated variant size and also to evaluate effect of
	      // splice variants on transcript. 
	      //
	      
	      if ( ( act_bp1 >= r_cds.subregion[s].start.position() && act_bp1 <= r_cds.subregion[s].stop.position() ) || ( act_bp2 >= r_cds.subregion[s].start.position() && act_bp2 <= r_cds.subregion[s].stop.position() ))
		{
		  inCDS = s;
		  for( unsigned int stmp = 0 ; stmp < r_cds.subregion.size(); stmp++ )
		    CDS_exons.insert(stmp);
		    }
	    }
	   
	  //
	  // Define if variant is within/overlapping an exon and which one
	  //

	  for ( unsigned int s  = 0 ; s < r_exon.subregion.size(); s++)
            {
	      if ( ( act_bp1 >= r_exon.subregion[s].start.position() && act_bp1 <= r_exon.subregion[s].stop.position() ) || ( act_bp2 >= r_exon.subregion[s].start.position() && act_bp2 <= r_exon.subregion[s].stop.position() ))
		inExon = s;
	      
	    }

	  //                                                                                                                                                                     
	  // identify UTR mutations                                                                                                                                              
	  //                                                                                                                                                                     
	  /******
	  if ( inExon > -1 && inCDS == -1 && r_cds.subregion.size() > 0 )
	    {
	      int cds_start = r_cds.subregion[0].start.position();
	      int cds_end = r_cds.subregion[r_cds.subregion.size()-1].stop.position();
	      
	      if ( ( act_bp1 < cds_start && positive_strand ) || ( act_bp2 > cds_end && negative_strand ) )
		annot.insert( SeqInfo( r->name , UTR5 ) );
	      if ( ( act_bp1 < cds_start && negative_strand ) || ( act_bp2 > cds_end && positive_strand ) )
		annot.insert( SeqInfo( r->name , UTR3 ) ); 
	    }
	  ******/
	  //
	  // Is this a SPLICE-SITE?	  
	  //
	  
	  for ( unsigned int s = 0 ; s < r_exon.subregion.size(); s++ )
	    {

	      sizeexonint =     r_exon.subregion[s].stop.position() - r_exon.subregion[s].start.position() + 1;
	      transtruncsize += r_exon.subregion[s].stop.position() - r_exon.subregion[s].start.position() + 1;
	      pos +=            r_exon.subregion[s].stop.position() - r_exon.subregion[s].start.position() + 1;
	      
	      if( s == r_exon.subregion.size() - 2 ) sizepenult = transtruncsize;
	      	      
	      bool splice = false;
	      int splicedist = 0;
	      
	      //
	      // Changed splicedist to be defined as - if intronic splice
	      // Changed to + if splice variant is DONOR  EXON 1 | INTRON | EXON2
	      //                                            -  321     321
	      //                                            +     123     123

	      //
	      // define splice deletion when a deletion spans the splice junction outside of the beginning of the first or the end of the last exon
	      //
	      	      
	      if ( act_ref.size() > act_alt.size() && positive_strand && ( ( s != 0 && act_bp1 < r_exon.subregion[s].start.position() && act_bp2 > r_exon.subregion[s].start.position() ) || ( s < r_exon.subregion.size()-1 && act_bp1 < r_exon.subregion[s].stop.position() && act_bp2 > r_exon.subregion[s].stop.position() ) ) )
		{
		  SeqInfo si = SeqInfo( r->name , SPLICEDEL );
		  si.alt = getc(*a);
		  annot.insert(si);
		}
	      else if ( s != first_exon  && ( r_exon.subregion[s].start.position() - act_bp1 <= 5 && r_exon.subregion[s].start.position() - act_bp1 >= -4 ) )

		{
		  // This is the modulo 3 of the transcript of all
		  // exons leading to the exon where splicing will
		  // occur.

	    	  int splicedtransc = (transtruncsize - sizeexonint) % 3;
	    	  
		  // This is the size of the spliced exon assuming
		  // exon skipping.

		  int splicedexon = sizeexonint % 3;

		  //  If they are equal to each other then we assume
		  //  out of frame protein truncated variant equals 0,
		  //  i.e. ofptv = 0. Else ofptv = 1.
		  
		  if ( splicedtransc != splicedexon )
		    {
		      // debugging
		      notinframe = 1;
		    }

		  if ( negative_strand )
		    {
		      int in_exonsp = r_exon.subregion.size() - (s);
		      
		      // donor intronic 2bp splice variants 100%
		      // conserved GT
		      
		      std::string iseq = seqdb->lookup( chr , r_exon.subregion[s].start.position() - ( 5 + act_ref.size() ) , r_exon.subregion[s].start.position() - 1 );
		      std::string eseq = seqdb->lookup( chr , r_exon.subregion[s].start.position() , r_exon.subregion[s].start.position() + ( 4 + act_ref.size() )  );
		      Helper::str2upper(iseq);
		      Helper::str2upper(eseq);
		      std::string var_iseq = iseq;
		      std::string var_eseq = eseq;

		      int d = r_exon.subregion[s].start.position() - act_bp1;
		      
		      if ( d > 0 )
			var_iseq.replace( var_iseq.length() - d , act_ref.size() , act_alt );
		      else
			var_eseq.replace(  abs(d) , act_ref.size() , act_alt );

		      iseq = getc( iseq ); var_iseq = getc( var_iseq );
                      eseq = getc( eseq ); var_eseq = getc( var_eseq );

		      // Donor essential conserved GT
		      if( ( var_iseq.at(var_iseq.length()-2) != 'T' || var_iseq.at(var_iseq.length()-1) != 'G' ) && iseq.at(iseq.length() - 2  ) == 'T' && iseq.at(iseq.length() - 1 ) == 'G')
			  {
			  SeqInfo si = SeqInfo( r->name , DONORIN2 );
			  si.splicedist = r_exon.subregion[s].start.position()  - act_bp1 ;
			  si.exin = in_exonsp;
			  si.ofptv = notinframe;
			  si.iseq = iseq;
			  si.eseq = eseq;
			  si.splice_type = "donor";
			  si.alt = getc(*a);
			  if( si.splicedist <= 0 ) si.splicedist--;
			  annot.insert(si);
			}
		      // Donor Intronic +45AG
		      else if((var_iseq.at(var_iseq.length() - 5) != 'G' || var_iseq.at(var_iseq.length() - 4) != 'A' ) && iseq.at(iseq.length() - 5 ) == 'G' && iseq.at(iseq.length() - 4 ) == 'A' ) 
			{
			  SeqInfo si = SeqInfo( r->name , DONORIN45AG );
			  si.splicedist = r_exon.subregion[s].start.position() - act_bp1 ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  si.iseq = iseq;
			  si.eseq = eseq;
			  si.splice_type = "donor";
			  si.alt = getc(*a);
			  if( si.splicedist <= 0 ) si.splicedist--;
			  annot.insert(si);
			}
		      // Donor Exonic AG
		      else if((var_eseq.at(0) != 'G' || var_eseq.at(1) != 'A' ) && eseq.at(0) == 'G' && eseq.at(1) == 'A')
			{
			  SeqInfo si = SeqInfo( r->name , DONOREX2AG ) ;
			  si.splicedist = r_exon.subregion[s].start.position() - act_bp1 ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  si.iseq = iseq;
			  si.eseq = eseq;
			  si.splice_type = "donor";
			  si.alt = getc(*a);
			  if( si.splicedist <= 0 ) si.splicedist--;
			  annot.insert( si );			  
			}
		      else
			{
			  // consider it splice generic
			  SeqInfo si = SeqInfo( r->name , SPLICE );
			  si.splicedist = r_exon.subregion[s].start.position() - act_bp1 ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  si.iseq = iseq;
			  si.eseq = eseq;
			  si.splice_type = "donor";
			  si.alt = getc(*a);
			  if( si.splicedist <= 0 ) si.splicedist--;
			  annot.insert( si );
			}
		    }
		  else // Positive Strand Acceptor
		    {
		      int in_exonsp = s+1;
		      std::string iseq = seqdb->lookup( chr , r_exon.subregion[s].start.position() - ( 5 + act_ref.size() ) , r_exon.subregion[s].start.position() - 1);
		      std::string eseq = seqdb->lookup( chr , r_exon.subregion[s].start.position() , r_exon.subregion[s].start.position() + ( 4 + act_ref.size() ) );
		      Helper::str2upper(iseq);
		      Helper::str2upper(eseq);
		      std::string var_iseq = iseq;
		      std::string var_eseq = eseq;

                      int d = r_exon.subregion[s].start.position() - act_bp1;
		      if ( d > 0 )
                        var_iseq.replace( var_iseq.length() - d  , act_ref.size() ,  act_alt );
                      else
                        var_eseq.replace( abs(d) , act_ref.size() , act_alt );

		      if( ( var_iseq.at(var_iseq.length()-2) != 'A' || var_iseq.at(var_iseq.length()-1) != 'G' ) && iseq.at(iseq.length() - 2) == 'A' && iseq.at(iseq.length() - 1 ) == 'G' )
			{
			  SeqInfo si = SeqInfo( r->name , ACCEPTORIN2 );
			  si.splicedist = act_bp1 - r_exon.subregion[s].start.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  si.iseq = iseq;
			  si.eseq = eseq;
			  si.splice_type = "acceptor";
			  si.alt = *a;
			  if( si.splicedist <= 0 ) si.splicedist--;
			  annot.insert( si );
			}
		      else if ( var_eseq.at(0) != 'G' && eseq.at(0) == 'G' )
			{
			  SeqInfo si = SeqInfo( r->name , ACCEPTOREX1G );
			  si.splicedist = act_bp1 - r_exon.subregion[s].start.position() ;
			  si.ofptv = notinframe;
			  si.exin = in_exonsp;
			  si.iseq = iseq;
                          si.eseq = eseq;
			  si.splice_type = "acceptor";
			  si.alt = *a;
			  if( si.splicedist <= 0 ) si.splicedist--;
			  annot.insert( si );
			}
		      else{
			SeqInfo si = SeqInfo( r->name , SPLICE );
			si.splicedist = act_bp1 - r_exon.subregion[s].start.position() ;
			si.exin = in_exonsp ;
			si.ofptv = notinframe ;
			si.iseq = iseq;
			si.eseq = eseq;
			si.splice_type = "acceptor";
			si.alt = *a;
			if( si.splicedist <= 0 ) si.splicedist--;
			annot.insert( si );
		      }
		    }
		}

	      else if ( s != last_exon && ( r_exon.subregion[s].stop.position() - act_bp1  <= 4 && r_exon.subregion[s].stop.position() - act_bp1  >= -5 ) )
		{
		  int splicedtransc = (transtruncsize - sizeexonint) % 3;
		  int splicedexon = sizeexonint % 3;
		  if ( splicedtransc != splicedexon )
		      notinframe = 1;
		  
		  if ( negative_strand )
		    {
		      int in_exonsp = r_exon.subregion.size() - (s);
		      std::string eseq = seqdb->lookup( chr , r_exon.subregion[s].stop.position() - ( 4 + act_ref.size() ) , r_exon.subregion[s].stop.position() ) ;
		      std::string iseq = seqdb->lookup( chr , r_exon.subregion[s].stop.position() + 1 , r_exon.subregion[s].stop.position() + ( 5 + act_ref.size() ) );
		      Helper::str2upper(iseq);
		      Helper::str2upper(eseq);
		      std::string var_iseq = iseq;
		      std::string var_eseq = eseq;

                      int d = act_bp1 - r_exon.subregion[s].stop.position();
                      if ( d <= 0 ) 
			var_eseq.replace( var_eseq.length() - abs(d) - 1 , act_ref.size() , act_alt );
		      else if ( d == 1 && act_alt.size() > act_ref.size() ) // if exonic or insertion before conserved splice junction 
			var_eseq.replace( var_eseq.length() - 1, act_ref.size() , act_alt );
		      else
                        var_iseq.replace( d - 1 , act_ref.size() ,  act_alt );


		      iseq = getc( iseq ); var_iseq = getc( var_iseq );
		      eseq = getc( eseq ); var_eseq = getc( var_eseq );

		      if ( ( var_iseq.at(0) != 'G' || var_iseq.at(1) != 'A' ) && iseq.at(0) == 'G' && iseq.at(1) == 'A' )
			{
			  SeqInfo si = SeqInfo( r->name , ACCEPTORIN2 );
			  si.splicedist = r_exon.subregion[s].stop.position() - act_bp1;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  si.splice_type = "acceptor";
			  si.eseq = eseq;
			  si.iseq = iseq;
			  si.alt = getrc(*a);
			  if( si.splicedist <= 0 ) si.splicedist--;
			  annot.insert( si );
			}
		      else if( var_eseq.at(var_eseq.length()-1) != 'G' && eseq.at(eseq.length()-1) == 'G' )
			{
			  SeqInfo si = SeqInfo( r->name , ACCEPTOREX1G ) ;
			  si.splicedist = act_bp1 - r_exon.subregion[s].stop.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  si.eseq = eseq;
			  si.iseq = iseq;
			  si.splice_type = "acceptor";
			  si.alt = getrc(*a);
			  if( si.splicedist <= 0 ) si.splicedist--;
			  annot.insert( si );
			}
		      else
			{
			  SeqInfo si = SeqInfo( r->name , SPLICE );
			  si.splicedist = act_bp1 - r_exon.subregion[s].stop.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  si.eseq = eseq;
			  si.iseq = iseq;
			  si.splice_type = "acceptor";
			  si.alt = getrc(*a);
			  if( si.splicedist <= 0 ) si.splicedist--;
			  annot.insert( si );
			}
		    }
		  else
		    {
		      int in_exonsp = s+1;
		      std::string eseq = seqdb->lookup( chr , r_exon.subregion[s].stop.position() -  ( 4 + act_ref.size() ) , r_exon.subregion[s].stop.position() );
		      std::string iseq = seqdb->lookup( chr , r_exon.subregion[s].stop.position() + 1 , r_exon.subregion[s].stop.position() + ( 5 + act_ref.size() ) );
		      Helper::str2upper(iseq);
		      Helper::str2upper(eseq);
		      std::string var_iseq = iseq;
		      std::string var_eseq = eseq;

                      int d = act_bp1 - r_exon.subregion[s].stop.position();
                      if ( d <= 0 )
                        var_eseq.replace( var_eseq.length() - abs(d) - 1 , act_ref.size() , act_alt );
                      else if ( d == 1 && act_alt.size() > act_ref.size() ) // if exonic or insertion before conserved splice junction
			var_eseq.replace( var_eseq.length() - 1, act_ref.size() , act_alt );
                      else
                        var_iseq.replace( d - 1 , act_ref.size() ,  act_alt );

		      if ( ( var_iseq.at(0) != 'G' || var_iseq.at(1) != 'T' ) && iseq.at(0) == 'G' && iseq.at(1) == 'T')
			{
			  SeqInfo si = SeqInfo( r->name , DONORIN2 ) ;
			  si.splicedist = act_bp1 - r_exon.subregion[s].stop.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  si.eseq = eseq;
			  si.iseq = iseq;
			  si.splice_type = "donor";
			  si.alt = getr(*a);
			  if( si.splicedist <= 0 ) si.splicedist--;
			  annot.insert( si );
			}
		      else if ( (var_eseq.at(var_eseq.length()-1) != 'G' || var_eseq.at(var_eseq.length()-2) != 'A' ) && eseq.at(eseq.length()-1) == 'G' && eseq.at(eseq.length()-2) == 'A')
			{
			  SeqInfo si = SeqInfo( r->name , DONOREX2AG ) ;
			  si.splicedist = act_bp1 - r_exon.subregion[s].stop.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  si.eseq = eseq;
			  si.iseq = iseq;
			  si.splice_type = "donor";
			  si.alt = getr(*a);
			  if( si.splicedist <= 0 ) si.splicedist--;
			  annot.insert( si ) ;
			}
		      else if ( (var_iseq.at(3) != 'A' || var_iseq.at(4) != 'G' ) && iseq.at(3) == 'A' && iseq.at(4) == 'G' ) 
			{
			  SeqInfo si = SeqInfo( r->name , DONORIN45AG ) ;
			  si.splicedist = act_bp1 - r_exon.subregion[s].stop.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  si.eseq = eseq;
			  si.iseq = iseq;
			  si.splice_type = "donor";
			  si.alt = getr(*a);
			  if( si.splicedist <= 0 ) si.splicedist--;
			  annot.insert( si );
			}
		      else
			{
			  SeqInfo si = SeqInfo( r->name , SPLICE );
			  si.splicedist = act_bp1 - r_exon.subregion[s].stop.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  si.eseq = eseq;
			  si.iseq = iseq;
			  si.splice_type = "donor";
			  si.alt = getr(*a);
			  if( si.splicedist <= 0 ) si.splicedist--;
			  annot.insert( si );
			}
		    }
		}	      
	    }


	  //	  
          // identify UTR mutations
          //                                                                                                                                                                                                                                                                                                                                        

          if ( inExon > -1 && inCDS == -1 && r_cds.subregion.size() > 0 )
            {
              int cds_start = r_cds.subregion[0].start.position();
              int cds_end = r_cds.subregion[r_cds.subregion.size()-1].stop.position();

	      if ( ( act_bp1 < cds_start && positive_strand ) || ( act_bp2 > cds_end && negative_strand ) )
                annot.insert( SeqInfo( r->name , UTR5 ) );
              if ( ( act_bp1 < cds_start && negative_strand ) || ( act_bp2 > cds_end && positive_strand ) )
		annot.insert( SeqInfo( r->name , UTR3 ) );
	      
	      ++ii;
	      continue;
            }


	  //
	  // Define Non coding elements as having exons but no coding sequence
	  // Allow for splice annotation above but call all mutations here
	  if ( r_cds.subregion.size() == 0 ){
	    if ( inExon > -1 )
	      annot.insert( SeqInfo( r->name , npcRNA ) );
	    else
	      annot.insert( SeqInfo( r->name , INTRON ) );
	    
	    ++ii;
	    continue;
	  }

	  if( inExon == -1 ){
	    annot.insert( SeqInfo( r->name , INTRON ) );
	    
	    ++ii; // next region                                                                                                                                                   
	    continue;	    
	  }
	   
	  /****
	  //
	  // If no exons attached, implies an intronic SNP (or splice site)
	  //
	  if ( CDS_exons.size() == 0 )
	    {
	      // Otherwise
	      annot.insert( SeqInfo( r->name , INTRON ) );
	      
	      ++ii; // next region
	      continue;
	    }
	  *****/

	  //
	  // Get reference sequence
	  //

	  std::string ref_cds;
	  std::set<int>::iterator i = CDS_exons.begin();
	  while ( i != CDS_exons.end() )
	    {
	      if ( negative_strand )
		{
		  if (  *i > first_exon )
		    first_exon = *i;
		}
	      else
		{
		  if ( *i < first_exon )
		    first_exon = *i;
		}
	      
	      ref_cds += seqdb->lookup( chr ,
					r_cds.subregion[ *i ].start.position(),
					r_cds.subregion[ *i ].stop.position() );
	      ++i;
	    }

	  //
	  // Get position of our transcript relative to start of gene
	  //

	  int exon = negative_strand ? r_cds.subregion.size()-1 : 0 ;
	  int pos_extracted_seq = 0;
	  int pos_whole_transcript = 0;

	  while ( exon != inCDS )
	    {
	      // Count all exons
	      pos_whole_transcript += r_cds.subregion[ exon ].stop.position()
		- r_cds.subregion[ exon ].start.position() + 1;
	      
	      // Count only exons extracted from seqdb
	      if ( CDS_exons.find( exon ) != CDS_exons.end() )
		  pos_extracted_seq += r_cds.subregion[ exon ].stop.position() - r_cds.subregion[ exon ].start.position() + 1;
	      exon += negative_strand ? -1 : +1;
	    }

	  // And also add in the distance into the containing exon
	  int pos_in_exon = negative_strand ?
	    r_cds.subregion[ inCDS ].stop.position() - act_bp1 + 1 :
	    act_bp1 - r_cds.subregion[ inCDS ].start.position() + 1 ;
	  
	  pos_extracted_seq += pos_in_exon;
	  pos_whole_transcript += pos_in_exon;
	  
	  // Determine reading frame
	  int frame = r_cds.subregion[ first_exon ].meta.get1_int( PLINKSeq::TRANSCRIPT_FRAME() ) ;
	  
	  // Build variant sequence
	  std::string var_allele = act_alt;
	  
	  // If transcript on -ve strand, get reverse complement
	  if ( negative_strand )
	    {
	      ref_cds    = getrc( ref_cds );
	      var_allele = getrc( var_allele );
	    }
	  
	  // Replace position with alternate allele
	  // Changed this to allow insertions and deletions
	  
	  if( pos_extracted_seq > ref_cds.size() )
	    pos_extracted_seq = ref_cds.size();
	  if( pos_extracted_seq < 1 )
	    pos_extracted_seq = 1;
	  
	  std::string var_cds = ref_cds;
	  var_cds.replace( pos_extracted_seq-1, act_ref.size() , var_allele );

	  // Are reference and variant sequences identical for this gene?
	  if ( ref_cds == var_cds )
	    {
	      annot.insert( SeqInfo( MONO ) );
	      // Next region
	      ++ii;
	      continue;
	    }


	  //
	  // Translate sequence, and populate codon table
	  //

	  std::vector<std::string> ref_codon;
	  std::vector<std::string> alt_codon;
	  std::string trans_ref = translate( ref_cds , frame, ref_codon );
	  std::string trans_var = translate( var_cds , frame, alt_codon );
	  
	  while( alt_codon.size() > ref_codon.size() )
	    ref_codon.push_back("_");
	  while( ref_codon.size() > alt_codon.size() )
            alt_codon.push_back("_");

	  //
	  // Synonymous change?
	  //
	  
	  if ( trans_ref == trans_var )
	    {
	      annot.insert( SeqInfo( r->name , SYN ) );
	      ++ii; // next region
	      continue;
	    }
	  
	  //
	  // For indels, will need a better check of synon than above...
	  //
	  int longest = trans_ref.size() > trans_var.size() ? trans_ref.size() : trans_var.size();
	  int transrefsize = trans_ref.size();
	  int transaltsize = trans_var.size();

	  while ( trans_var.size() < longest ) { trans_var += "_"; }
	  while ( trans_ref.size() < longest ) { trans_ref += "_"; }
	  
	  //
	  // Make calls
	  //

	  std::vector<std::string> difs;
	  
	  // newpos_stop used indicates at which codon position the
	  // new reading frame ends in a stop (*). The position of the
	  // stop is calculated starting at the first changed amino
	  // acid that is created by the frame shift, and ending at
	  // the first stop codon (fs*#), e.g. p.Arg97Glyfs*16

	  int newpos_stop = 0;
	  int posinframe_indel = 0;
	  int firstfs_codon = 0;
	  int firststop_codon = 0;
	  int newpos_start = 0;
	  int origpepsize = 0;
	  int newpepsize = 0;
	  int isnmd = 0;
	  int isofptv = 0;
	  int pposfs = 0;
	  
	  // changed this from trans_var.size() to longest	  
	  for ( unsigned int i = 0; i < longest; i++ ) 
	    {     
	      // for reference -- for substitutions,
	      // ref allele = ref_cds.substr( pos_extracted_seq-1 , 1 )
	      // alt allele = var_allele
	      
	      // Found the new start codon:
	      if ( newpos_start == -1 && trans_var[i] == 'M' )
		newpos_start = i+1;
	      
	      if ( trans_ref[i] != trans_var[i] ) 
		{
		  
		  // Single base substitutions ("SNP"), or multi-nucleotide polymorphism ("MNP")
		  // Either way, the REF and the ALT are of equal lengths here:
		  
		  if ( act_ref.size() == act_alt.size() ) 
		    {
		      /* NOTE: since the REF and ALT are same length, then can correctly use 'i' (i+1, with 1-based offset)
			 as the protein position (and thus properly account for MNPs spanning multiple codons): */
		      
		      if ( i == 0 ) // The first codon has been changed, so the functioning "start" codon will be downstream in the aa sequence:
			newpos_start = -1;
		      
		      if ( trans_var[i] != '*' && trans_ref[i] != '*' ) // MISSENSE
			{
			  seq_annot_t type = MIS;
			  
			  if (trans_var[i] == '?') // an unknown ALT base resulted in an unknown AA
			    type = EXONIC_UNKNOWN;
			  
			  origpepsize = longest;
			  newpepsize = longest;
			  SeqInfo si = SeqInfo( r->name ,
						type ,
						reference ,
						*a ,
						pos_whole_transcript ,
						ref_codon[i] ,
						alt_codon[i] ,
						i+1 ,
						trans_ref.substr(i,1) ,
						trans_var.substr(i,1) ,
						0 ,
						origpepsize ,
						newpepsize );
			  annot.insert( si );
			}
		      else if ( trans_var[i] == '*' ) // NONSENSE
			{
			  seq_annot_t type = NON;
			  origpepsize = trans_ref.size();
			  newpepsize = i+1;
			  SeqInfo si = SeqInfo( r->name ,
						type ,
						reference ,
						*a ,
						pos_whole_transcript ,
						ref_codon[i] ,
						alt_codon[i] ,
						i+1 ,
						trans_ref.substr(i,1) ,
						trans_var.substr(i,1) ,
						0 ,
						origpepsize ,
						newpepsize );
			  si.ofptv = 1;
			  if(newpepsize*3 < sizepenult - 50 )
			    si.nmd = 1;
			  annot.insert( si );
			}
		      else if ( trans_ref[i] == '*' )   // READ-THROUGH
			{
			  seq_annot_t type = RT;
			  origpepsize = trans_ref.size();
			  newpepsize = longest;
			  SeqInfo si =  SeqInfo( r->name ,
						 type ,
						 reference ,
						 *a ,
						 pos_whole_transcript ,
						 ref_codon[i] ,
						 alt_codon[i] ,
						 i+1 ,
						 trans_ref.substr(i,1) ,
						 trans_var.substr(i,1) ,
						 0 ,
						 origpepsize ,
						 newpepsize );
			  annot.insert( si );
			}
		    }
		  else 
		    {

		      //
		      // Indel changes
		      //
		      
		      // Is frameshift?
		      
		      int dif = act_ref.size() > act_alt.size() ? act_ref.size() - act_alt.size() : act_alt.size() - act_ref.size();
		      
		      if ( dif % 3 ) 
			{
			  newpos_stop++;
			  
			  if ( newpos_stop == 1 ) 
			    {
			      firstfs_codon = i;
			      pposfs = i+1;
			    }
			  
			  // Stop of new transcript.
			  if (trans_var[i] == '*' && firststop_codon == 0 ) 
			    {

			      firststop_codon++;

			      seq_annot_t type = FRAMESHIFT;
			      origpepsize = trans_ref.size();
			      newpepsize = i+1;
			      newpos_stop -= 1;
			      
			      // Adding this so we can add ofptv and nmd predictions.
			      
			      SeqInfo si = SeqInfo( r->name ,
						    type ,
						    reference ,
						    *a ,
						    pos_whole_transcript ,
						    ref_codon[firstfs_codon] ,
						    alt_codon[firstfs_codon] ,
						    pposfs ,
						    trans_ref.substr(firstfs_codon,1) ,
						    trans_var.substr(firstfs_codon,1) ,
						    newpos_stop ,
						    origpepsize ,
						    newpepsize );
			      si.ofptv = 1;
			      
			      if(newpepsize*3 < sizepenult - 50 )
				si.nmd = 1;
			      
			      annot.insert( si );
			    }
			  
			  if (i == longest - 1 && firststop_codon == 0 ) 
			    
			    {
			      firststop_codon++;
			      seq_annot_t type = FRAMESHIFT;
			      
			      origpepsize = trans_ref.size();
			      newpepsize = longest;
			      newpos_stop -= 1;
			      
			      SeqInfo si = SeqInfo( r->name ,
						    type ,
						    reference ,
						    *a ,
						    pos_whole_transcript ,
						    ref_codon[firstfs_codon] ,
						    alt_codon[firstfs_codon] ,
						    pposfs ,
						    trans_ref.substr(firstfs_codon,1) ,
						    trans_var.substr(firstfs_codon,1) ,
						    newpos_stop ,
						    origpepsize ,
						    newpepsize );
			      si.ofptv = 1;
			      
			      // Here the protein made is too long and would not undergo NMD.
			      annot.insert( si );
			    }
			}


		      //
		      // Non-FS indels: Codon-deletion
		      //
		    
		      else if ( act_ref.size() > act_alt.size() )
			{
			  
			  posinframe_indel++;
			  
			  if( posinframe_indel == 1) 
			    {
			      seq_annot_t type = CODONDELETION;
			      if( trans_ref[i] == '*' )
				type = STOPDELETION;
			      if( (pos_whole_transcript-1) % 3 != 0 )
				type = OOFCODONDELETION;
			      
			      annot.insert( SeqInfo( r->name ,
						     type ,
						     reference ,
						     *a ,
						     pos_whole_transcript ,
						     ref_codon[i] ,
						     alt_codon[i] ,
						     (int)floor(((pos_whole_transcript-1)/3.0)+1) ,
						     trans_ref.substr(i,1) ,
						     trans_var.substr(i,1)) );
			    }
			}

		      
		      //
		      // Non-FS indels: Codon-insertion
		      //
		      
		      else if ( act_alt.size() > act_ref.size() )  // [nb. should always be true if get to this point]
			{
			  
			  posinframe_indel++;
			  
			  if ( posinframe_indel == 1 ) 
			    {
			      seq_annot_t type = CODONINSERTION;
			      if ( trans_var[i] == '*' )
				type = STOPINSERTION;
			      if ( (pos_whole_transcript-1) % 3 != 0 ){
				type = OOFCODONINSERTION;
			      }
			      
			      annot.insert( SeqInfo( r->name ,
						     type ,
						     reference ,
						     *a ,
						     pos_whole_transcript ,
						     ref_codon[i] ,
						     alt_codon[i] ,
						     (int)floor(((pos_whole_transcript-1)/3.0)+1) ,
						     trans_ref.substr(i,1) ,
						     trans_var.substr(i,1) ) );
			    }
			}
		    }
		  
		}
	    
	    
	      //
	      // Edge case for synon indel (i.e. rare case at end of transcript)
	      //
	    
	      if ( trans_ref[i] == '*' && trans_var[i] == '*' && ( act_ref.size() != act_alt.size() ) )
		{
		  
		  seq_annot_t type = INDEL;
		  
		  annot.insert( SeqInfo( r->name ,
					 type ,
					 reference ,
					 *a ,
					 pos_whole_transcript ,
					 ref_codon[i] ,
					 alt_codon[i] ,
					 (int)floor(((pos_whole_transcript-1)/3.0)+1) ,
					 trans_ref.substr(i,1) ,
					 trans_var.substr(i,1) ) );
		  
		  i = longest;
		}
	    }
	  


	  //
	  // Start-lost
	  //
	  
	  if (newpos_start != 0) 
	    {
	      seq_annot_t type = SL;
	      
		  newpepsize = 0;
		  if (newpos_start != -1)
		    newpepsize = transrefsize - newpos_start + 1;
		  
		  origpepsize = transrefsize;
		  
		  SeqInfo si = SeqInfo( r->name ,
					type ,
					reference ,
					*a ,
					pos_whole_transcript ,
					ref_codon[0] ,
					alt_codon[0] ,
					1,
					trans_ref.substr(0,1) ,
					trans_var.substr(0,1) ,
					newpos_start ,
					origpepsize ,
				  newpepsize );
		  
		  annot.insert( si );
	    }
	  
	  
	  ++ii;
	} // next transcript
      // next alternate allele
      ++a;      
    }

  //
  // other annotations that are usually within locdb regions but not
  // coding, make UNDEFINED at present
  //

  if( annot.size() == 0 )
    annot.insert(SeqInfo( UNDEF ) );

  return annot;  
}
 
std::string SeqInfo::codon() const
{
  if ( intergenic() || intronic() || codondeletion() || codoninsertion() || frameshift() ) return ".";
  return cpos1 == 0 ?
    "." :
    "c." + Helper::int2str( cpos1 ) + ref_seq + ">" + alt_seq ;
}

std::string SeqInfo::genomic() const
{
  if ( intergenic() || intronic() ) return ".";

  return cpos1 == 0 ?
    "." :
    "g." + Helper::int2str( cpos1 ) + genomic_ref + ">" + genomic_alt ;
}


void SeqInfo::details( Variant & var ) const
{
  
  if ( splice() || esplice() || csplice() )
    {
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "DIST=" + Helper::int2str( splicedist ) );
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "NMD=" + Helper::int2str( nmd ) );
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "OFPTV=" + Helper::int2str( ofptv ) );
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "EXIN=" + Helper::int2str( exin ) );
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "ISEQ=" + iseq );
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "ESEQ=" + eseq );
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "SPLICE_TYPE=" + splice_type );
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "ALT=" + alt );
    }
  
  if ( nonsense() )
    {
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "NMD=" + Helper::int2str( nmd ) );
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "OFPTV=" + Helper::int2str( ofptv )  );    	  
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "PEPSIZE=" + Helper::int2str( origpepsize ) + "->" + Helper::int2str(newpepsize) );
    }
  
  if ( frameshift() ) 
    {
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "NMD=" + Helper::int2str( nmd ) );
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "OFPTV=" + Helper::int2str( ofptv )  );    	  
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "FSX=" + Helper::int2str( fs_stop )  );    	  
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "PEPSIZE=" + Helper::int2str( origpepsize ) + "->" + Helper::int2str(newpepsize) );    	  
    }

  if ( startlost() ) 
    {
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "NEWAUG=" + Helper::int2str( fs_stop ) );    	        
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "PEPSIZE=" + Helper::int2str( origpepsize ) + "->" + Helper::int2str(newpepsize)  );    	        
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "OFPTV=" + Helper::int2str( ofptv )  );    	        
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "NMD=" + Helper::int2str(nmd) );    	        
    }

  if ( readthrough() ) 
    {      
      var.meta.add( PLINKSeq::ANNOT_DETAILS() , "PEPSIZE=" + Helper::int2str( origpepsize ) + "->" + Helper::int2str(newpepsize) ) ;
    }
}



std::string SeqInfo::protein() const
{
  
  if ( frameshift() || codondeletion() || codoninsertion() ) return ".";
  
  if ( intergenic() || intronic() ) return ".";
  
  return ppos1 == 0 ?
    "." :
    "p." + Helper::int2str( ppos1 ) + ref_aa + ">" + alt_aa ;
}


void Annotate::init()
{
  if ( GP == NULL ) return ;
  rmap.clear();
  transcript_group_id = 0;  
  seqdb = &(GP->seqdb);
  if ( ! db ) setDB( LOCDB );
  
  
  MetaInformation<VarMeta>::field( PLINKSeq::ANNOT()         , META_TEXT ,  1 , "Highest impact annotation across transcripts" );
  MetaInformation<VarMeta>::field( PLINKSeq::ANNOT_TYPE()    , META_TEXT , -1 , "Transcript-specific list of annotations" );
  MetaInformation<VarMeta>::field( PLINKSeq::ANNOT_GENE()    , META_TEXT , -1 , "List of transcripts used for annotation" );
  MetaInformation<VarMeta>::field( PLINKSeq::ANNOT_CODON()   , META_TEXT , -1 , "Transcript-specific codon changes" );
  MetaInformation<VarMeta>::field( PLINKSeq::ANNOT_PROTEIN() , META_TEXT , -1 , "Transcript-specific amino-acid changes" );
  MetaInformation<VarMeta>::field( PLINKSeq::ANNOT_CODING()  , META_INT  , -1 , "Transcript-specific coding status (0/1)" );
  MetaInformation<VarMeta>::field( PLINKSeq::ANNOT_DETAILS() , META_TEXT , -1 , "Misc. details" );
  MetaInformation<VarMeta>::field( PLINKSeq::ANNOT_SUMMARY() , META_TEXT ,  1 , "Annotation summary" );
 
}



std::string Annotate::translate_reference( const Region & region , bool verbose )
{
  
  if ( ( ! db ) || ( ! db->attached() ) ) return "";
  if ( ( ! seqdb ) || ( ! seqdb->attached() ) ) return "";
  
  // Assume that we have exons as sub-regions
  
  if ( region.subregion.size() == 0 ) return "";
  
  

  // Get strand
  std::map<int,int> cds2ex;
  int n_exons = 0;
  for (int e=0;e<region.subregion.size();e++)
    {
      if ( region.subregion[e].CDS() ) { cds2ex[ n_exons ] = e; ++n_exons ; }
    }

  int s = region.subregion[0].meta.get1_int( PLINKSeq::TRANSCRIPT_STRAND() );
  if ( s == 0 ) Helper::halt( "could not get strand information in translate_reference" );

  bool negative_strand = s < 0 ;
  
  int first_exon = negative_strand ? 0 : n_exons - 1;


  // Get reference coding sequence

  std::string ref_cds;

  std::vector<int> exon_start_idx(n_exons,0);

  for (int i=0; i < region.subregion.size(); i++)
    {
      if ( ! region.subregion[i].CDS() ) continue;
      ref_cds += seqdb->lookup( region.chromosome() ,
				region.subregion[ i ].start.position(),
				region.subregion[ i ].stop.position() );
    }

  // Get scaffold

  if ( negative_strand )
    {
      for (int i=n_exons-1; i >= 0; i--)
	{
	  if ( i == n_exons-1 ) exon_start_idx[n_exons-1] = 0;
	  else exon_start_idx[i] = exon_start_idx[i+1] + ( region.subregion[ cds2ex[i]+1 ].stop.position() - region.subregion[ cds2ex[i]+1 ].start.position() );
	}
    }
  else
    {
      for (int i=0; i < n_exons; i++)
	{
	  if ( i == 0 ) exon_start_idx[0] = 0;
	  else exon_start_idx[i] = exon_start_idx[i-1] + ( region.subregion[ cds2ex[i]-1 ].stop.position() - region.subregion[ cds2ex[i]-1 ].start.position() );
	}
    }

  // Determine reading frame is always 0, as getting whole gene

  int frame = 0;


  // If transcript on -ve strand, get reverse complement

  if ( negative_strand )
    ref_cds = getrc( ref_cds );

  // Translate sequence, and populate codon

  std::vector<std::string> ref_codon;

  std::string trans_ref = translate( ref_cds , frame, ref_codon );
  
  if ( ! verbose ) return trans_ref;


  //
  // Create a verbose representation
  //

  Helper::halt( "needs CDS() changes to work -- internal note to self.." );

  // assume output stream called 'loci' exists
  
  Out & pout = Out::stream( "loci" );
  Out & pbase = Out::stream( "bstats" );
  

  //  std::stringstream verb;

  if ( negative_strand )
    {
      for (int i=n_exons-1; i >= 0; i--)
	{
	  pout << "exon " << i+1 << "\t"
	       << region.chromosome() << "\t"
	       << region.subregion[ i ].start.position() << "\t"
	       << region.subregion[ i ].stop.position() << "\t"
	       << exon_start_idx[i] << "\n";
	}
    }
  else
    {
      for (int i=0; i < n_exons; i++)
	{
	  pout << "exon " << i+1 << "\t"
	       << region.chromosome() << "\t"
	       << region.subregion[ i ].start.position() << "\t"
	       << region.subregion[ i ].stop.position() << "\t"
	       << exon_start_idx[i] << "\n";
	}
    }


  int j = 0;
  std::string chr = Helper::chrCode( region.chromosome() );
  int x = negative_strand ? n_exons-1 : 0 ;
  int exon_start_genomic;

  // i is position in transcript, always goes 0..(n-1) as transcript will
  // have already been put in RC if necessary

  // consider each exon

  int i = 0; // offset in transcript, in bases
  int p = 0; // offset in transcript, in AA

  int bp1 = 0;

  int codon_from_exon = x;


  while ( 1 )
    {

      // pout << ">>>--- exon " << x+1 << " of " << n_exons << " ----\n";

      int offset = 0;

      int start = negative_strand ?
	region.subregion[ x ].stop.position() :
	region.subregion[ x ].start.position() ;

      int len = region.subregion[ x ].stop.position() - region.subregion[ x ].start.position() + 1;

      int dir = negative_strand ? -1 : +1;

      // base-stats (only for sites on a whole exon)
      std::map<std::string,int> bpcnt;
      std::map<std::string,int> aacnt;


      for (int j=0; j<len; j++)
	{
	  int mod = (i+1) % 3 ;

	  if ( mod == 2 ) // start of codon
	    {
	      codon_from_exon = x;
	      bp1 = start + ( j * dir );
	    }
	  else if ( (i+1) % 3 == 0 )
	    {
	      if ( codon_from_exon == x )
		{

		  pout << chr << "\t"
		       << (x+1) << "/" << region.subregion.size() << "\t"
		       << ( negative_strand ? "rev" : "fwd" ) << "\t"
		       << bp1 << ".."
		       << start + ( j * dir ) <<  "\t"
		       << (i-1) << ".." << i+1 << "\t"
		       << p+1 << "\t"
		       << ref_codon[p] << "\t"
		       << trans_ref.substr(p,1) << "\n";
		  
		  bpcnt[ ref_codon[p].substr(0,1) ]++;
		  bpcnt[ ref_codon[p].substr(1,1) ]++;
		  bpcnt[ ref_codon[p].substr(2,1) ]++;
		  aacnt[ trans_ref.substr(p,1) ]++;

		}
	      else
		{
		  // codon is split across two exons
		  pout << chr << "\t"
		       << (x+1) << "&" << (codon_from_exon+1) << "/" << region.subregion.size() << "\t"
		       << ( negative_strand ? "rev" : "fwd" ) << "\t"
		       << bp1 << ".."
		       << start + ( j * dir ) <<  "\t"
		       << (i-1) << ".." << i+1 << "\t"
		       << p+1 << "\t"
		       << ref_codon[p] << "\t"
		       << trans_ref.substr(p,1) << "\n";
		  
		}
	      ++p;
	    }
	  
	  ++i;

	}


      // base-stats output
      
      pbase << "_BSTATS\t"
	    << region.name << "\t"
	    << region.coordinate() << "\t"
	    << x << "\t"
	    << ( negative_strand ? "rev" : "fwd" ) << "\t";
      
      pbase << bpcnt["A"] << "\t"
	    << bpcnt["C"] << "\t"
	    << bpcnt["G"] << "\t"
	    << bpcnt["T"] << "\t"
	    << bpcnt["a"] + bpcnt["c"] + bpcnt["g"] + bpcnt["t"] + bpcnt["N"] ;
      
      std::map<std::string,std::string>::const_iterator pp = Annotate::aa.begin();
      while ( pp != Annotate::aa.end() )
	{
	  pbase << "\t" << aacnt[ pp->first ] ;
	  ++pp;
	}
      pbase << "\n";

      
      // next exon;
      x += dir;

      // done?
      if ( x < 0 || x == n_exons ) break;
    }

  return trans_ref;
}



//
// Helper functions
//

std::string Annotate::translate(std::string & seq, int frame , std::vector<std::string> & codons )
{
  Helper::str2upper(seq);

  // TODO Need to check logic here
  if ( seq.size() - frame == 1 ) seq += "-";
  else if ( seq.size() - frame == 2 ) seq += "--";
  
  std::string trans = "";
  codons.clear();
  
  for (unsigned int i = frame; i<seq.size(); i+=3)
    {
      std::string codon = seq.substr( i,3 );
      codons.push_back( codon );
      
      if ( codon.find("-") != std::string::npos )
	trans += "i";
      else
	{
	  std::string tmp = t[ codon ];
	  if ( tmp == "" ) tmp = "?";
	  trans += tmp;
	}
    }
  return trans;
}


std::string Annotate::getrc(const std::string & s)
{
  int sz = s.size();
  std::string r;
  for ( int i = 0 ; i < sz ; i++ )
    {
      if      ( s[i] == 'a' ) r += "t";
      else if ( s[i] == 'c' ) r += "g";
      else if ( s[i] == 'g' ) r += "c";
      else if ( s[i] == 't' ) r += "a";
      else if ( s[i] == 'A' ) r += "T";
      else if ( s[i] == 'C' ) r += "G";
      else if ( s[i] == 'G' ) r += "C";
      else if ( s[i] == 'T' ) r += "A";
      else r += "N";
    }
  reverse( r.begin(), r.end() );
  return r;
}

std::string Annotate::getr(const std::string & s)
{
  std::string r = s;
  reverse( r.begin(), r.end() );
  return r;
}

std::string Annotate::getc(const std::string & s)
{
  int sz = s.size();
  std::string r;
  for ( int i = 0 ; i < sz ; i++ )
    {
      if      ( s[i] == 'a' ) r += "t";
      else if ( s[i] == 'c' ) r += "g";
      else if ( s[i] == 'g' ) r += "c";
      else if ( s[i] == 't' ) r += "a";
      else if ( s[i] == 'A' ) r += "T";
      else if ( s[i] == 'C' ) r += "G";
      else if ( s[i] == 'G' ) r += "C";
      else if ( s[i] == 'T' ) r += "A";
      else r += "N";
    }
  return r;
}

