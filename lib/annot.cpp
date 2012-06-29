
#include "plinkseq/annot.h"
#include "plinkseq/variant.h"
#include "plinkseq/gstore.h"
#include "plinkseq/filemap.h"

#include <algorithm>
#include <cmath>

extern GStore * GP;

using namespace std;
using namespace Helper;

// Stop codons indicated as "*"

std::map<uint64_t,Region> Annotate::rmap;
uint64_t Annotate::transcript_group_id = 0;

std::map<seq_annot_t,std::string> populate_seqinfo()
{
  std::map<seq_annot_t,std::string> m;
  m[UNDEF]    		= ".";
  m[MONO]     		= "monomorhpic";
  m[IGR]     	    = "intergenic-region";
  m[INTRON]  	    = "intronic";
  m[UTR5]     		= "5-UTR";
  m[UTR3]     		= "3-UTR";
  m[SYN]      		= "silent";
  m[MIS]     		= "missense";
  m[INDEL] 	  		= "indel";
  m[CODONDELETION]  = "codondeletion";
  m[CODONINSERTION] = "codoninsertion";
  m[NON]      		= "nonsense";
  m[SL]	      		= "startlost";
  m[PART]     		= "partial-codon";
  m[DONORIN2] 		= "splice-donor-in2";
  m[DONOREX2AG]     = "splice-donor-ex2ag";
  m[ACCEPTOREX1G]   = "splice-acceptor-ex1g";
  m[ACCEPTORIN2]    = "splice-acceptor-in2";
  m[DONORIN45AG]    = "splice-donor-in45ag";
  m[SPLICE]   		= "splice";
  m[FRAMESHIFT]     = "frameshift";
  m[RT]       		= "readthrough";
  return m;
}

std::map<std::string,std::string> populate_t()
{
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

std::string Annotate::translate(std::string & seq, int frame , std::vector<std::string> & codons )
{

    Helper::str2upper(seq);

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


bool Annotate::load_transcripts( fType t , const std::string & name )
{
  init();
  setDB(t);
  if ( ( ! db ) || ( ! db->attached() ) ) return false;
  uint64_t id = db->lookup_group_id( name );
  if ( id == 0 ) return false;
  return load_transcripts( db->lookup_group_id( name ) );
}



bool Annotate::load_transcripts( uint64_t id )
{
  init();
  transcript_group_id = 0;
  rmap.clear();
  if ( ( ! db ) || ( ! db->attached() ) ) return false ;
  if ( id == 0 ) return false;
  set<Region> regions = db->get_regions( id );

  set<Region>::iterator i = regions.begin();
  while ( i != regions.end() )
    {
      rmap[ i->id ] = *i;
      ++i;
    }
  transcript_group_id = id;
  return true;
}

bool Annotate::load_transcripts( const std::string & grp , const std::set<Region> & regions )
{
  init();
  setDB( LOCDB );
  if ( ( ! db ) || ( ! db->attached() ) ) return false;
  uint64_t id = db->lookup_group_id( grp );
  if ( id == 0 ) return false;
  return load_transcripts( id , regions );
}


bool Annotate::load_transcripts( uint64_t id , const std::set<Region> & regions )
{
  init();
  transcript_group_id = 0;
  rmap.clear();
  if ( ( ! db ) || ( ! db->attached() ) ) return false ;
  if ( id == 0 ) return false;
  std::set<Region>::iterator i = regions.begin();
  while ( i != regions.end() )
    {
      rmap[ i->id ] = *i;
      ++i;
    }
  transcript_group_id = id;
  return true;
}


bool Annotate::annotate(Variant & var , Region * pregion )
{

  std::set<SeqInfo> s;
  if ( ! pregion )
    s = annotate( var.chromosome() ,
		  var.position() ,
		  var.alternate() ,
		  var.reference() ) ;
  else
    s = annotate( var.chromosome() ,
		  var.position() ,
		  var.alternate() ,
		  var.reference() ,
		  *pregion ) ;

  bool annot = false;

  // summary of 'worst' annotation (int, syn, nonsyn)

  int is_silent 		= 0;
  int is_missense 		= 0;
  int is_codondeletion 	= 0;
  int is_codoninsertion = 0;
  int is_frameshift 	= 0;
  int is_indel 			= 0;
  int is_splice 		= 0;
  int is_esplice 		= 0;
  int is_nonsense 		= 0;
  int is_startlost		= 0;
  int is_readthrough 	= 0;
  int is_intergenic 	= 0;
  int is_intronic 		= 0;

  std::set<SeqInfo>::iterator i = s.begin();
  while ( i != s.end() )
    {

      // track whether this is coding, for a 'single' return code

      if ( i->synon() ) ++is_silent;
      if ( i->frameshift() ) ++is_frameshift;
      if ( i->codondeletion() ) ++is_codondeletion;
      if ( i->codoninsertion() ) ++is_codoninsertion;
      if ( i->missense() ) ++is_missense;
      if ( i->nonsense() ) ++is_nonsense;
      if ( i->startlost() ) ++is_startlost;
      if ( i->splice() ) ++is_splice;
      if ( i->esplice() ) ++is_esplice;
      if ( i->readthrough() ) ++is_readthrough;
      if ( i->intergenic() ) ++is_intergenic;
      if ( i->intronic() ) ++is_intronic;
      if ( i->indel() ) ++is_indel;

      // add annotations
      var.meta.add( PLINKSeq::ANNOT_TYPE() , i->status() );
      var.meta.add( PLINKSeq::ANNOT_GENE() , i->gene_name() );
      var.meta.add( PLINKSeq::ANNOT_CODING() , i->coding() );
 // var.meta.add( PLINKSeq::ANNOT_EXONIC() , i->exonic() );

      var.meta.add( PLINKSeq::ANNOT_CHANGE() , i->genomic() );
      var.meta.add( PLINKSeq::ANNOT_CODON() , i->codon() );

      // for splice, use this slot for the details, for now
      if ( i->splice() || i->esplice() )
	{
	  std::string s = "dist=" + Helper::int2str( i->splicedist );
	  std::string s2 = "nmd=" + Helper::int2str( i->nmd );
	  std::string s3 = "ofptv=" + Helper::int2str( i->ofptv );
	  std::string s4 = "exin=" + Helper::int2str( i->exin );
	  var.meta.add( PLINKSeq::ANNOT_PROTEIN() , s );
	  var.meta.add( PLINKSeq::ANNOT_PROTEIN(), s2 );
	  var.meta.add( PLINKSeq::ANNOT_PROTEIN(), s3 );
	  var.meta.add( PLINKSeq::ANNOT_PROTEIN(), s4 );
	}
      if ( i->frameshift() || i->nonsense() )
      {
    	  std::string s = "nmd=" + Helper::int2str( i->nmd );
    	  std::string s2 = "ofptv=" + Helper::int2str( i->ofptv );
    	  var.meta.add( PLINKSeq::ANNOT_PROTEIN() , s );
    	  var.meta.add( PLINKSeq::ANNOT_PROTEIN() , s2 );
    	  var.meta.add( PLINKSeq::ANNOT_PROTEIN(), i->protein() );

      }
      else
	{
	  var.meta.add( PLINKSeq::ANNOT_PROTEIN() , i->protein() );
	}
      ++i;
    }

  // what is the 'worst' annotation?
  std::string aworst = "";
  if ( is_frameshift ) aworst = "frameshift";
  else if ( is_nonsense ) aworst = "nonsense";
  else if ( is_startlost ) aworst = "startlost";
  else if ( is_esplice ) aworst = "esplice";
  else if ( is_readthrough ) aworst = "readthrough";
  else if ( is_codoninsertion ) aworst = "codoninsertion";
  else if ( is_codondeletion ) aworst = "codondeletion";
  else if ( is_missense ) aworst = "missense";
  else if ( is_indel ) aworst = "indel";
  else if ( is_splice ) aworst = "splice";
  else if ( is_silent ) aworst = "silent";
  else if ( is_intronic ) aworst = "intronic";
  else aworst = "intergenic";

  var.meta.set( PLINKSeq::ANNOT() , aworst );

  // what is the 'summary/consensus' annotation?
  int acount = 0;

  if ( is_silent ) ++acount;
  if ( is_missense ) ++acount;
  if ( is_nonsense ) ++acount;
  if ( is_startlost ) ++acount;
  if ( is_esplice ) ++acount;
  if ( is_splice ) ++acount;
  if ( is_readthrough ) ++acount;
  if ( is_intergenic ) ++acount;
  if ( is_frameshift ) ++acount;
  if ( is_codoninsertion ) ++acount;
  if ( is_codondeletion ) ++acount;

  std::string annot_summary = aworst;
  if ( acount > 1 ) annot_summary = "mixed";

  annot_summary += ",NON=" + Helper::int2str( is_nonsense );
  annot_summary += ",SL="  + Helper::int2str( is_startlost );
  annot_summary += ",MIS=" + Helper::int2str( is_missense );
  annot_summary += ",SYN=" + Helper::int2str( is_silent );
  annot_summary += ",ESPL=" + Helper::int2str( is_esplice );
  annot_summary += ",SPL=" + Helper::int2str( is_splice );
  annot_summary += ",RTH=" + Helper::int2str( is_readthrough );
  annot_summary += ",INT=" + Helper::int2str( is_intronic );
  annot_summary += ",IGR=" + Helper::int2str( is_intergenic );
  annot_summary += ",FRAMESHIFT="  + Helper::int2str( is_frameshift );
  annot_summary += ",INDEL="  + Helper::int2str( is_indel );
  annot_summary += ",CODONINSERTION="  + Helper::int2str( is_codoninsertion );
  annot_summary += ",CODONDELETION="  + Helper::int2str( is_codondeletion );

  var.meta.set( PLINKSeq::ANNOT_SUMMARY() , annot_summary );

  // did we receive any annotation?
  return is_silent || is_startlost || is_nonsense || is_missense || is_splice || is_esplice || is_readthrough || is_frameshift || is_codondeletion || is_indel || is_codoninsertion;

}


std::set<SeqInfo> Annotate::lookup(Variant & var)
{
  return annotate( var.chromosome() ,
		   var.position() ,
		   var.alternate() ,
		   var.reference() ) ;

}

std::set<SeqInfo> Annotate::annotate( int chr,
				      int bp1 ,
				      const std::string & alternate ,
				      const std::string & reference ,
				      const Region & r )
{
  std::set<Region> t;
  t.insert(r);
  return annotate(chr,bp1,alternate,reference,&t);
}

std::set<SeqInfo> Annotate::annotate( int chr,
				      int bp1 ,
				      const std::string & alternate ,
				      const std::string & reference ,
				      const std::set<Region> * pregions )
{


  //
  // Store annotations generated in here
  //

  std::set<SeqInfo> annot;

  //
  // Are we attached to the relevant databases?
  //

  if ( ! pregions )
    {
      if ( ( ! db ) || ( ! db->attached() ) ) return annot;
      if ( rmap.size() == 0 ) return annot;
    }

  if ( ( ! seqdb ) || ( ! seqdb->attached() ) ) return annot;


  //
  // Assumptions:
  //

  // At any one position, we assume an individual has either a

  // *Either* a SNP, insertion or deletion:

  // snp .. just flip the base

  //  deletions .. obliterates relevant bases in cds splice range, but does
  //  not otherwise impact cds start, stop, splicing, etc.

  //  insertions .. adds bases within cds splice range, but does not otherwise
  //  impact cds start, stop, splicing, etc. .. that is, the insertion must
  //  fall after the first base of a cds exon and before the last base of that
  //  cds exon


  //
  // Currently, only annot single-base pair polymorphisms
  // Annotation to indels added - prelim code
  //

  // Changes made here by manny commented out
  // if ( reference.size() > 1 )
  //  {
  //   annot.insert( SeqInfo( UNDEF ) ) ;
  //  return annot;
  //  }

  // changes stop

  //
  // Get all transcripts that overlap this position
  //

    int bp2 = alternate.size() > reference.size() ? bp1 + alternate.size() - 1 : bp1 + reference.size() - 1 ;

  // currently bp2 not really used / assume == bp1 (i.e. substitution)
  // uncommented line at top and commented line below
   // int bp2 = bp1 + reference.size() - 1 ;
  // changes stop

  // Pull regions to use for annotation from database, if they weren't supplied directly

  std::set<Region> regions;

  if ( ! pregions )
    {
      std::vector<uint64_t> reg_id = db->get_region_ids( transcript_group_id , chr , bp1 , bp2 );

      for ( unsigned int i=0; i<reg_id.size(); i++)
	{
	  std::map<uint64_t,Region>::const_iterator r = rmap.find( reg_id[i] );
	  if ( r != rmap.end() ) regions.insert( r->second );
	}

      pregions = &regions;
    }


  // Variant does not fall into any transcript

  if ( pregions->size() == 0 )
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

      // std::cout << "Annotating... " << chr << " " << bp1 << " " << *a <<  "\n";

      // If insertion we want to annotate - changes made by manny commented below....
      //      if ( a->size() > 1 )
      //	{
      //	  annot.insert( SeqInfo( UNDEF ) ) ;
      //	  ++a;
      //	  continue;
      //	}
      
      
      std::set<Region>::iterator r = pregions->begin();

      while ( r != pregions->end() )
	{

	  //
	  // Assume that we have exons encoded as sub-regions
	  //

	  if ( r->subregion.size() == 0 )
	    {
	      ++r;
	      continue;
	    }


	  //
	  // Distinguish here between CDS and non-CDS exons
	  //

	  Region r_exon; // only contains CDS (and extra stop-codon)
	  Region r_cds;  // contains all exons (i.e., will include 3UTR , 5UTR, but not 'stop')

	  for ( int ss=0;ss < r->subregion.size();ss++)
	    {
	      
	      if      ( r->subregion[ss].CDS() || r->subregion[ss].stop_codon() )
		r_cds.subregion.push_back( r->subregion[ss] );
	      else if ( r->subregion[ss].exon() )
		r_exon.subregion.push_back( r->subregion[ss] );

	  }

	  // Which exon(s) does this mutation impact?  Pull in
	  // neighbouring exon if needed. Assume all subregions are on
	  // the same chromosome
	  
	  std::set<int> CDS_exons;
	  int in_CDS_exon = 0;
	  
	  int pos = 1;
	  int in_exon = 0;
	  
	  //added this to verify whether or not splice variant is in frame
	  
	  int notinframe = 0;
	  
	  
	  //
	  // Strand and exon status
	  //

	  bool negative_strand = false;
	  bool positive_strand = false; //just to check at least one strand is given
	  
	  for (int ss=0;ss< r_cds.subregion.size(); ss++)
	    {
	      int s = r_cds.subregion[ss].meta.get1_int( PLINKSeq::TRANSCRIPT_STRAND() );
	      if ( s == 0 ) continue;
	      negative_strand = s < 0 ;
	      positive_strand = s > 0 ;
	      break;
	    }
	  
	  // skip if no strand information
	  
	  if ( ! ( negative_strand || positive_strand ) )
	    {
	      ++r;
	      continue;
	    }
	  
	  
	  int first_exon = negative_strand ? 0 : r_cds.subregion.size()-1;
	  
	  int last_exon = negative_strand ? r_cds.subregion.size()-1 : 0;



	  //
	  // Does variant fall within an exon, or near an intron/exon splice-site boundary
	  //
	
	  //size of exons before exon skipping ... Here we assume exon
	  //skipping for splice variants. We want to distinguish
	  //between splice variants that lead to NMD and those that
	  //escape NMD but maintain proper translation frame. MR
	  
	  int transtruncsize = 0;

	  //want to keep track of the size of the exon i am interested
	  //in. MR

	  int sizeexonint = 0;

	  //want to keep track of the size of the transcript where
	  //final splicing occurs. Cases where this will not work:
	  //Exon Number = 1. Have not investigated NMD in genes with
	  //exon size = 1 -- under investigation. MR
	  
	  int sizepenult = 0;

	  for ( unsigned int s  = 0 ; s < r_cds.subregion.size(); s++)
	    {
	      
	      // let us keep all the exons to compute protein
	      // truncated variant size and also to evaluate effect of
	      // splice variants on transcript. Shaun you may have a
	      // better idea of how to deal with this.MR
		
	      bool bb = true;
	      
	      if( bb )
		{
		  // only if my variant is within an exon boundary
		  // will i examine all the exons. This may not be
		  // necessary and may have already been checked. I am
		  // not sure if this is true so being a bit
		  // cautious. MR
			
		  if ( bp1 >= r_cds.subregion[s].start.position() &&
		       bp1 <= r_cds.subregion[s].stop.position() )
		    {
		      // variant is in"_"CDS_ exon s.
		      
		      
		      in_CDS_exon = s;
		      for( unsigned int stmp = 0 ; stmp < r_cds.subregion.size(); stmp++ )
			{
			  CDS_exons.insert(stmp);
			}
		    }
		  
		}
	    }
	  
	  
	  // Is this a SPLICE-SITE?
	  
	  for ( unsigned int s = 0 ; s < r_exon.subregion.size(); s++ )
	    {
	      
	      sizeexonint =     r_exon.subregion[s].stop.position() - r_exon.subregion[s].start.position() + 1;
	      transtruncsize += r_exon.subregion[s].stop.position() - r_exon.subregion[s].start.position() + 1;
	      pos +=            r_exon.subregion[s].stop.position() - r_exon.subregion[s].start.position() + 1;
	      
	      if( s == r_exon.subregion.size() - 2 ){
		sizepenult = transtruncsize;
		
	      }
	      
	      bool splice = false;
	      int splicedist = 0;
	      
	      // Changed splicedist to be defined as - if intronic splice
	      // Changed to + if splice variant is DONOR  EXON 1 | INTRON | EXON2
	      //                                            -  321     321
	      //                                            +     123     123

	     

	      if ( s != first_exon  &&  abs( r_exon.subregion[s].start.position() - bp1 ) <= 5 )
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
		      
		      if ( r_exon.subregion[s].start.position()  - bp1 <= 2 && r_exon.subregion[s].start.position() - bp1  > 0 )
			{
			  
			  SeqInfo si = SeqInfo( r_exon.name , DONORIN2 );
			  si.splicedist = r_exon.subregion[s].start.position()  - bp1 ;
			  si.exin = in_exonsp;
			  si.ofptv = notinframe;
			  if ( si.splicedist > 0 ) annot.insert(si);
			  
			  
			}
			 
		      // Donor Intronic +45AG
			
		      else if( ( r_exon.subregion[s].start.position() - bp1 == 4 || r_exon.subregion[s].start.position() - bp1 == 5 ) 
			       && getrc( seqdb->lookup( chr , r_exon.subregion[s].start.position() - 5 , r_exon.subregion[s].start.position() - 4) ) == "AG" )
			{
			  SeqInfo si = SeqInfo( r_exon.name , DONORIN45AG );
			  si.splicedist = r_exon.subregion[s].start.position() - bp1 ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  if( si.splicedist > 0 ) annot.insert(si);
			}
		      
		      // Donor Exonic AG
		      else if ( bp1 - r_exon.subregion[s].start.position() <= 1 && bp1 - r_exon.subregion[s].start.position() >= 0 
				&& getrc( seqdb->lookup( chr , r_exon.subregion[s].start.position(), r_exon.subregion[s].start.position() + 1 ) ) == "AG" )
			{
			  SeqInfo si = SeqInfo( r_exon.name , DONOREX2AG ) ;
			  si.splicedist = r_exon.subregion[s].start.position() - bp1 ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  if ( si.splicedist <= 0 ) --si.splicedist;
			  if ( si.splicedist <= 5 ) annot.insert( si );			  
			}
		      
		      else
			{
			  // consider it splice generic
			  
			  SeqInfo si = SeqInfo( r_exon.name , SPLICE );
			  si.splicedist = r_exon.subregion[s].start.position() - bp1 ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  if ( si.splicedist <= 0 ) --si.splicedist ;
			  if ( si.splicedist >= -5 ) annot.insert( si );
			  
			}
		      
		    }
		  else
		    {

		      int in_exonsp = s+1;
			
		      if ( r_exon.subregion[s].start.position() - bp1 <= 2 && r_exon.subregion[s].start.position() - bp1 > 0 )
			{
			  SeqInfo si = SeqInfo( r_exon.name , ACCEPTORIN2 );
			  si.splicedist = bp1 - r_exon.subregion[s].start.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  if ( si.splicedist <= 5 ) annot.insert( si );
			  
			}
		      else if ( r_exon.subregion[s].start.position() - bp1 == 0 && seqdb->lookup( chr , bp1 , bp1 ) == "G" )
			{
			  SeqInfo si = SeqInfo( r_exon.name , ACCEPTOREX1G );
			  si.splicedist = bp1 - r_exon.subregion[s].start.position() ;
			  si.ofptv = notinframe;
			  si.exin = in_exonsp;
			  if ( si.splicedist == 0 ) ++si.splicedist;
			  if ( si.splicedist <= 5 ) annot.insert( si );
			}
		      else{
			SeqInfo si = SeqInfo( r_exon.name , SPLICE );
			si.splicedist = bp1 - r_exon.subregion[s].start.position() ;
			si.exin = in_exonsp ;
			si.ofptv = notinframe ;
			if ( si.splicedist >= 0 ) ++si.splicedist;
			if ( si.splicedist <= 5 ) annot.insert( si );
		      }
		      		      
		    }
		}

	      if ( s != last_exon && abs( r_exon.subregion[s].stop.position() - bp1 ) <= 5 )
		{

		  int splicedtransc = (transtruncsize - sizeexonint) % 3;
		  int splicedexon = sizeexonint % 3;
		  if ( splicedtransc != splicedexon )
		    {		      
		      notinframe = 1;
		    }
		  
		  if ( negative_strand )
		    {
		      
		      int in_exonsp = r_exon.subregion.size() - (s);
			
		      if ( bp1 - r_exon.subregion[s].stop.position()  <= 2 && bp1 - r_exon.subregion[s].stop.position() > 0 )
			{
			  SeqInfo si = SeqInfo( r_exon.name , ACCEPTORIN2 );
			  si.splicedist = r_exon.subregion[s].stop.position() - bp1;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  if ( si.splicedist <= 5 ) annot.insert( si );
			  
			}
		      else if( bp1 - r_exon.subregion[s].stop.position() == 0 
			       && getrc( seqdb->lookup( chr , r_exon.subregion[s].stop.position() , r_exon.subregion[s].stop.position() ) ) == "G" )
			{
			  SeqInfo si = SeqInfo( r_exon.name , ACCEPTOREX1G ) ;
			  si.splicedist = bp1 - r_exon.subregion[s].stop.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  if ( si.splicedist <= 0 ) --si.splicedist;
			  if ( si.splicedist <= 5 ) annot.insert( si );
			  
			}
		      else
			{
			  SeqInfo si = SeqInfo( r_exon.name , SPLICE );
			  si.splicedist = bp1 - r_exon.subregion[s].stop.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  if ( si.splicedist <= 0 ) --si.splicedist;
			  if ( si.splicedist <= 5 ) annot.insert( si );
			}
		      
		    }
		  else
		    {
		      int in_exonsp = s+1;
		      
		      if ( bp1 - r_exon.subregion[s].stop.position() <= 2 && bp1 - r_exon.subregion[s].stop.position() > 0 )
			{
			  SeqInfo si = SeqInfo( r_exon.name , DONORIN2 ) ;
			  si.splicedist = bp1 - r_exon.subregion[s].stop.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  if ( si.splicedist >= -5 ) annot.insert( si );
			}
		      else if ( r_exon.subregion[s].stop.position() - bp1 <= 1 && r_exon.subregion[s].stop.position() - bp1 >= 0 
				&& seqdb->lookup( chr , r_exon.subregion[s].stop.position() - 1, r_exon.subregion[s].stop.position() ) == "AG" )
			{
			  SeqInfo si = SeqInfo( r_exon.name , DONOREX2AG ) ;
			  si.splicedist = bp1 - r_exon.subregion[s].stop.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  if ( si.splicedist <= 0 ) --si.splicedist ;
			  if ( si.splicedist >= -5 ) annot.insert( si ) ;
			  
			}
		      else if ( ( bp1 - r_exon.subregion[s].stop.position() == 5 || bp1 - r_exon.subregion[s].stop.position() == 4 ) 
				&& seqdb->lookup( chr , r_exon.subregion[s].stop.position() + 4, r_exon.subregion[s].stop.position() + 5 ) == "AG" )
			{
			  SeqInfo si = SeqInfo( r_exon.name , DONORIN45AG ) ;
			  si.splicedist = bp1 - r_exon.subregion[s].stop.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  if ( si.splicedist >= -5 ) annot.insert( si );
			}
		      else
			{
			  SeqInfo si = SeqInfo( r_exon.name , SPLICE );
			  si.splicedist = bp1 - r_exon.subregion[s].stop.position() ;
			  si.ofptv = notinframe ;
			  si.exin = in_exonsp ;
			  if ( si.splicedist <= 0 ) --si.splicedist;
			  if ( si.splicedist >= -5 ) annot.insert( si );
			}
		    }
		}	      
	    }
	  
	  
	  //
	  // If no exons attached, implies an intronic SNP (or splice site)
	  //

	  if ( CDS_exons.size() == 0 )
	    {
	      // Otherwise
	      annot.insert( SeqInfo( r->name , INTRON ) );
	      ++r; // next region
	      continue;
	    }
	  
	  
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

	  while ( exon != in_CDS_exon )
	    {

	      // Count all exons
	      pos_whole_transcript += r_cds.subregion[ exon ].stop.position()
		- r_cds.subregion[ exon ].start.position() + 1;
	      
	      // Count only exons extracted from seqdb
	      if ( CDS_exons.find( exon ) != CDS_exons.end() )
		{
		  pos_extracted_seq += r_cds.subregion[ exon ].stop.position() - r_cds.subregion[ exon ].start.position() + 1;
		}
	      exon += negative_strand ? -1 : +1;
	    }
	  

	  //
	  // And also add in the distance into the containing exon
	  //

	  int pos_in_exon = negative_strand ?
	    r_cds.subregion[ in_CDS_exon ].stop.position() - bp1 + 1 :
	    bp1 - r_cds.subregion[ in_CDS_exon ].start.position() + 1 ;

	  pos_extracted_seq += pos_in_exon;
	  pos_whole_transcript += pos_in_exon;
	  

	  //
	  // Determine reading frame
	  //

	  int frame = r_cds.subregion[ first_exon ].meta.get1_int( PLINKSeq::TRANSCRIPT_FRAME() ) ;


	  //
	  // Build variant sequence
	  //

	  std::string var_allele = *a;


	  //
	  // If transcript on -ve strand, get reverse complement
	  //

	  if ( negative_strand )
	    {
	      ref_cds    = getrc( ref_cds );
	      var_allele = getrc( var_allele );
	    }


	  //
	  // Replace position with alternate allele
	  //

	  // Changed this to allow insertions and deletions
	  
	  std::string var_cds = ref_cds;

	  if ( reference.size() > 1 )
	    {	      
	      var_cds.replace( pos_extracted_seq-1 , reference.size() , var_allele );
	    }
	  else if ( a->size() > 1 )
	    {
	      var_cds.replace( pos_extracted_seq-1 , 1 , var_allele );	      
	    }
	  else
	    {
	      var_cds.replace( pos_extracted_seq-1 , 1 , var_allele );
	    }
	  
	  std::cout << "(2)so far, so good\n";
	  
	  //
	  // Are reference and variant sequences identical for this gene?
	  //

	  if ( ref_cds == var_cds )
	    {
	      
	      annot.insert( SeqInfo( MONO ) );
	      
	      // Next region
	      ++r;
	      continue;
	    }


	  //
	  // Translate sequence, and populate codon
	  //

	  std::vector<std::string> ref_codon;
	  std::vector<std::string> alt_codon;

	  std::string trans_ref = translate( ref_cds , frame, ref_codon );
	  std::string trans_var = translate( var_cds , frame, alt_codon );
  //
	  // Synomous change?
	  //

	  if ( trans_ref == trans_var )
	    {
	      annot.insert( SeqInfo( r->name , SYN ) );
	      ++r; // next region
	      continue;
	    }


	  //
	  // For indels, will need a better check of synon that above...
	  //

	  int longest = trans_ref.size() > trans_var.size()
	    ? trans_ref.size() : trans_var.size() ;
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

	  for ( unsigned int i=0; i< longest; i++ )
	    {

	      // for reference -- for substitutions,
	      // ref allele = ref_cds.substr( pos_extracted_seq-1 , 1 )
	      // alt allele = var_allele
		
	      if ( newpos_start > 0 )
		{
		  newpos_start++;
		  seq_annot_t type = SL;
		  newpepsize = transrefsize - newpos_start + 1;
		  origpepsize = transrefsize;
		  
		  if ( trans_var[i] == 'M' )
		    {
		      SeqInfo si = SeqInfo( r->name ,
					    type ,
					    reference ,
					    *a ,
					    pos_whole_transcript ,
					    ref_codon[0] ,
					    alt_codon[0] ,
					    (int)floor(((pos_whole_transcript-1)/3.0)+1) ,
					    trans_ref.substr(0,1) ,
					    trans_var.substr(0,1) ,
					    newpos_start ,
					    origpepsize ,
					    newpepsize );
		      
		      annot.insert( si );
		      
		    }
		  
		}

	      if ( trans_ref[i] != trans_var[i] )
		{
	    	  // Single base substitution changes - manny.

		  if( reference.size() == alternate.size() && reference.size() == 1 )
		    {
		      if ( i == 0 )
			{
			  
			  newpos_start++;
			  seq_annot_t type = SL;
			  // Used for debugging to detect start lost sites.

			}
		      
		      else
			{
			  if ( trans_var[i] != '*' && trans_ref[i] != '*'  )
			    {
			      seq_annot_t type = MIS;
			      origpepsize = longest;
			      newpepsize = longest;
			      SeqInfo si = SeqInfo( r->name ,
						    type ,
						    reference ,
						    *a ,
						    pos_whole_transcript ,
						    ref_codon[i] ,
						    alt_codon[i] ,
						    (int)floor(((pos_whole_transcript-1)/3.0)+1) ,
						    trans_ref.substr(i,1) ,
						    trans_var.substr(i,1) ,
						    0 ,
						    origpepsize ,
						    newpepsize );
			      annot.insert( si );
			    }

			  if ( trans_var[i] == '*' )
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
						    (int)floor(((pos_whole_transcript-1)/3.0)+1) ,
						    trans_ref.substr(i,1) ,
						    trans_var.substr(i,1) ,
						    0 ,
						    origpepsize ,
						    newpepsize );
			      si.ofptv = 1;
			      if(newpepsize*3 < sizepenult - 50 ){
				si.nmd = 1;
				
				
				
			      }
			      annot.insert( si );
			      
			  }
			  else if ( trans_ref[i] == '*' )
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
						     (int)floor(((pos_whole_transcript-1)/3.0)+1) ,
						     trans_ref.substr(i,1) ,
						     trans_var.substr(i,1) ,
						     0 ,
						     origpepsize ,
						     newpepsize );
			      annot.insert( si );
			    }
			  
			  
			}
		    }
		  else
 		    {

		      // Indel changes - manny.
		      // Frameshift Indels.

		      // if not modulo 3 consider frameshift -- look
		      // at assumptions made in beg. We will need to
		      // change this for indels that start at splicing
		      // regions? Consider them splice in the meantime

		      int modtmpr = ( reference.size() - 1 ) % 3;
		      int modtmpa = ( alternate.size() - 1 ) % 3;

		      if (  ( reference.size() > 1 && modtmpr != 0 ) || ( alternate.size() > 1 && ( modtmpa != 0 ) ) )
			{
			  
			  newpos_stop++;
			  if(newpos_stop == 1)
			    {
			      firstfs_codon = i;
			      pposfs = i+1;
			    }
			  // Stop of new transcript.
			  if(trans_var[i] == '*' && firststop_codon == 0)
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
			      
			      if(newpepsize*3 < sizepenult - 50 ){
				si.nmd = 1;
				
			      }

			      //uncommented the line below - change back to comment
			      
			      annot.insert( si );
			      			      
			      
			    }
			  
			  if (i == longest - 1 && firststop_codon == 0) 
			    {
			      firststop_codon++;
			      seq_annot_t type = FRAMESHIFT;
			      
			      // uncommented the line below - change back to comment
			      // have to change this for frameshift indels that make elongated transcripts
			      
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

			      //Here the protein made is too long and would not undergo NMD.
			      
			      annot.insert( si );
			      
			    }
			}
		      else if ( reference.size() > 1 && modtmpr % 3 == 0)
			{
			  
			  
			  posinframe_indel++;
			  if( posinframe_indel == 1)
			    {
 			      seq_annot_t type = CODONDELETION;
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
		      else if ( alternate.size() > 1 && modtmpa % 3 == 0){
			
			posinframe_indel++;
			if( posinframe_indel == 1 )
			  {
			    seq_annot_t type = CODONINSERTION;
			    
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
		      else 
			{
			  
			  seq_annot_t type = INDEL;
			
			  
			}
		      
		      
 		    }
		  
		  
		  
		}
	    }
	  
	  ++r;
	} // next transcript
      
      // next alternate allele
      ++a;
      
    }

  return annot;
  
}



std::string SeqInfo::codon() const
{
  if ( intergenic() || intronic() ) return ".";
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

std::string SeqInfo::protein() const
{
 if ( frameshift() ) return "p." + ref_aa + Helper::int2str( ppos1 ) + alt_aa + ",fsX=" + Helper::int2str( fs_stop ) + ",pepsize=" + Helper::int2str( origpepsize ) + "_" + Helper::int2str(newpepsize) + ",ofptv=" + Helper::int2str( ofptv ) + ",nmd=" + Helper::int2str(nmd) ;
 if ( startlost() ) return "p." + ref_aa + Helper::int2str( ppos1 ) + alt_aa + ",newAUG=" + Helper::int2str( fs_stop ) + ",pepsize=" + Helper::int2str( origpepsize ) + "_" + Helper::int2str(newpepsize) + ",ofptv=" + Helper::int2str( ofptv ) + ",nmd=" + Helper::int2str(nmd) ;
 if ( nonsense() ) return "p." + Helper::int2str( ppos1 ) + ref_aa + ">" + alt_aa + ",pepsize=" + Helper::int2str( origpepsize ) + "_" + Helper::int2str(newpepsize) + ",ofptv=" + Helper::int2str( ofptv ) + ",nmd=" + Helper::int2str(nmd)  ;
 if ( readthrough() ) return "p." + Helper::int2str( ppos1 ) + ref_aa + ">" + alt_aa +  ",pepsize=" + Helper::int2str( origpepsize ) + "_" + Helper::int2str(newpepsize) ;
 if ( intergenic() || intronic() ) return ".";
  return ppos1 == 0 ?
    "." :
    "p." + Helper::int2str( ppos1 ) + ref_aa + ">" + alt_aa ;
}


void Annotate::init()
{
  rmap.clear();
  transcript_group_id = 0;

  seqdb = &(GP->seqdb);
  if ( ! db ) setDB( LOCDB );

  // Add slot for additions

  MetaInformation<VarMeta>::field( "_ANNOT" , META_INT , 1 , "Annotation" );

  MetaInformation<VarMeta>::field( "_SYN" , META_TEXT , -1 , "Synonymous allele" );
  MetaInformation<VarMeta>::field( "_MIS" , META_TEXT , -1 , "Missense allele" );
  MetaInformation<VarMeta>::field( "_NON" , META_TEXT , -1 , "Nonsense allele" );
  MetaInformation<VarMeta>::field( "_PART" , META_TEXT , -1 , "Partial codon" );
  MetaInformation<VarMeta>::field( "_SPLICE" , META_TEXT , -1 , "Splice-site" );
  MetaInformation<VarMeta>::field( "_ESPLICE" , META_TEXT , -1 , "Essential splice-site" );
  MetaInformation<VarMeta>::field( "_SL" , META_TEXT , -1 , "Start-lost" );
  MetaInformation<VarMeta>::field( "_INTRON" , META_TEXT , -1 , "Intronic") ;
  MetaInformation<VarMeta>::field( "_FRAMESHIFT" , META_TEXT , -1 , "Frameshift allele");
// made changes here for indels
  MetaInformation<VarMeta>::field( "_CODONDELETION" , META_TEXT , -1 , "Codon-deletion allele");
  MetaInformation<VarMeta>::field( "_CODONINSERTION" , META_TEXT , -1 , "Codon-insertion allele");
 // change stops here
  MetaInformation<VarMeta>::field( "_READTHROUGH" , META_TEXT , -1 , "Read-through allele");
  MetaInformation<VarMeta>::field( "_5UTR" , META_TEXT , -1 , "5' UTR" );
  MetaInformation<VarMeta>::field( "_3UTR" , META_TEXT , -1 , "3' UTR" );
  MetaInformation<VarMeta>::field( "_IGR" , META_TEXT , -1 , "Intergenic region");
  MetaInformation<VarMeta>::field( "_MONO" , META_TEXT , -1 , "Monomorphic");

}



std::string Annotate::translate_reference( const Region & region , bool verbose )
{

  if ( ( ! db ) || ( ! db->attached() ) ) return "";
  if ( ( ! seqdb ) || ( ! seqdb->attached() ) ) return "";

  // Assume that we have exons as sub-regions

  if ( region.subregion.size() == 0 ) return "";

  // Get strand

  int n_exons = region.subregion.size();

  bool negative_strand = region.subregion[0].meta.get1_int( PLINKSeq::TRANSCRIPT_STRAND() ) == -1;

  int first_exon = negative_strand ? 0 : n_exons - 1;


  // Get reference coding sequence

  std::string ref_cds;

  std::vector<int> exon_start_idx(n_exons,0);

  for (int i=0; i < n_exons; i++)
    {
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
	  else exon_start_idx[i] = exon_start_idx[i+1] + ( region.subregion[ i+1 ].stop.position() - region.subregion[ i+1 ].start.position() );
	}
    }
  else
    {
      for (int i=0; i < n_exons; i++)
	{
	  if ( i == 0 ) exon_start_idx[0] = 0;
	  else exon_start_idx[i] = exon_start_idx[i-1] + ( region.subregion[ i-1 ].stop.position() - region.subregion[ i-1 ].start.position() );
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

  //  std::stringstream verb;

  if ( negative_strand )
    {
      for (int i=n_exons-1; i >= 0; i--)
	{
	  plog << "exon " << i+1 << "\t"
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
	  plog << "exon " << i+1 << "\t"
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

      // plog << ">>>--- exon " << x+1 << " of " << n_exons << " ----\n";

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

		  plog << chr << "\t"
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
		  plog << chr << "\t"
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

      plog << "_BSTATS\t"
		<< region.name << "\t"
		<< region.coordinate() << "\t"
		<< x << "\t"
		<< ( negative_strand ? "rev" : "fwd" ) << "\t";

      plog << bpcnt["A"] << "\t"
		<< bpcnt["C"] << "\t"
		<< bpcnt["G"] << "\t"
		<< bpcnt["T"] << "\t"
		<< bpcnt["a"] + bpcnt["c"] + bpcnt["g"] + bpcnt["t"] + bpcnt["N"] ;

      std::map<std::string,std::string>::const_iterator pp = Annotate::aa.begin();
      while ( pp != Annotate::aa.end() )
	{
	  plog << "\t" << aacnt[ pp->first ] ;
	  ++pp;
	}
      plog << "\n";




      // next exon;
      x += dir;

      // done?
      if ( x < 0 || x == n_exons ) break;
    }

  return trans_ref;
}
