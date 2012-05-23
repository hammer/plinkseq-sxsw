
#include "annot.h"
#include "variant.h"
#include "gstore.h"
#include "filemap.h"

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
  m[UNDEF]    = ".";
  m[MONO]     = "monomorhpic";
  m[IGR]      = "intergenic";
  m[INTRON]   = "intronic";
  m[UTR5]     = "UTR-5";
  m[UTR3]     = "UTR-3";
  m[SYN]      = "silent";
  m[MIS]      = "missense";
  m[NON]      = "nonsense";
  m[PART]     = "partial-codon";
  m[SPLICE5]  = "splice-5";
  m[SPLICE3]  = "splice-3";
  m[ESPLICE5] = "esplice-5"; 
  m[ESPLICE3] = "esplice-3";
  m[FS]       = "frameshift";
  m[RT]       = "readthrough";
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
  
  int is_silent      = 0;
  int is_missense    = 0;
  int is_splice      = 0;
  int is_esplice     = 0;
  int is_nonsense    = 0;
  int is_readthrough = 0;
  int is_intergenic  = 0;
  int is_intronic    = 0;
  int is_utr         = 0;

  std::set<SeqInfo>::iterator i = s.begin();
  while ( i != s.end() )
    {
      
      // track whether this is coding, for a 'single' return code      
      
      if ( i->synon() )       ++is_silent;
      if ( i->missense() )    ++is_missense;
      if ( i->nonsense() )    ++is_nonsense;
      if ( i->splice() )      ++is_splice;
      if ( i->esplice() )     ++is_esplice;
      if ( i->utr() )         ++is_utr;
      if ( i->readthrough() ) ++is_readthrough;
      if ( i->intergenic() )  ++is_intergenic;
      if ( i->intronic() )    ++is_intronic;
    
      // add annotations      
      var.meta.add( PLINKSeq::ANNOT_TYPE() ,   i->status() );
      var.meta.add( PLINKSeq::ANNOT_GENE() ,   i->gene_name() );
      var.meta.add( PLINKSeq::ANNOT_CODING() , i->coding() );
      // var.meta.add( PLINKSeq::ANNOT_EXONIC() ,   i->exonic() );
      
      var.meta.add( PLINKSeq::ANNOT_CHANGE() , i->genomic() );
      var.meta.add( PLINKSeq::ANNOT_CODON() ,  i->codon() );

      // for splice, use this slot for the details, for now
      if ( i->splice() )
	{
	  std::string s = "dist=" + Helper::int2str( i->splicedist );
	  var.meta.add( PLINKSeq::ANNOT_PROTEIN() , s );	  
	}
      else
	{
	  var.meta.add( PLINKSeq::ANNOT_PROTEIN() , i->protein() );
	}
      ++i;
    }
  
  // what is the 'worst' annotation?
  std::string aworst = "";

  if ( is_nonsense ) aworst = "nonsense";
  else if ( is_esplice ) aworst = "esplice";
  else if ( is_splice ) aworst = "splice";
  else if ( is_readthrough ) aworst = "readthrough";
  else if ( is_missense ) aworst = "missense";
  else if ( is_silent ) aworst = "silent";
  else if ( is_intronic ) aworst = "intronic";
  else aworst = "intergenic"; 

  var.meta.set( PLINKSeq::ANNOT() , aworst );
  
  // what is the 'summary/consensus' annotation?
  int acount = 0;

  if ( is_silent ) ++acount;
  if ( is_missense ) ++acount;
  if ( is_nonsense ) ++acount;
  if ( is_splice ) ++acount;
  if ( is_esplice ) ++acount;
  if ( is_readthrough ) ++acount;
  if ( is_intergenic ) ++acount;
  
  std::string annot_summary = aworst;
  if ( acount > 1 ) annot_summary = "mixed";  

  annot_summary += ",NON=" + Helper::int2str( is_nonsense );
  annot_summary += ",MIS=" + Helper::int2str( is_missense );
  annot_summary += ",SYN=" + Helper::int2str( is_silent );
  annot_summary += ",SPL=" + Helper::int2str( is_splice );
  annot_summary += ",ESP=" + Helper::int2str( is_esplice );
  annot_summary += ",RTH=" + Helper::int2str( is_readthrough );
  annot_summary += ",INT=" + Helper::int2str( is_intronic );
  annot_summary += ",IGR=" + Helper::int2str( is_intergenic );

  var.meta.set( PLINKSeq::ANNOT_SUMMARY() , annot_summary );
  
  // did we receive any annotation?
  return is_silent || is_nonsense || is_missense || is_splice || is_readthrough || is_esplice;
  
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
  // 

  if ( reference.size() > 1 )  
    {
      annot.insert( SeqInfo( UNDEF ) ) ;
      return annot; 
    }
  

  //
  // Get all transcripts that overlap this position
  //
  
  //   int bp2 = alternate.size() > reference.size() ? bp1 + seq.size() - 1 : bp1 + ref.size() - 1 ; 
  
  // currently bp2 not really used / assume == bp1 (i.e. substitution)
  int bp2 = bp1 + reference.size() - 1 ;


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
      
      if ( a->size() > 1 )
	{
	  annot.insert( SeqInfo( UNDEF ) ) ;
	  ++a; 
	  continue;
	}

      
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
	  
	  Region r_exon;   // only contains CDS (and extra stop-codon)
	  Region r_cds;    // contains all exons (i.e., will include 3UTR, 5UTR, but not 'stop')
	  
	  for (int ss=0;ss< r->subregion.size(); ss++)
	    {
	      if      ( r->subregion[ss].CDS() || r->subregion[ss].stop_codon() ) 
		r_cds.subregion.push_back( r->subregion[ss] );
	      else if ( r->subregion[ss].exon() ) 
		r_exon.subregion.push_back( r->subregion[ss] );
	    }

	  
	  //
	  // Which exon(s) does this mutation impact?  Pull in
	  // neighbouring exon if needed. Assume all subregions are on
	  // the same chromosome
	  //


	  std::set<int> CDS_exons;
	  int in_CDS_exon = 0;
	  

	  //
	  // Strand and exon status
	  //
	  
	  // note that strand encodes type of exon also
	  
	  bool negative_strand = false;
	  bool positive_strand = false; // just to check at least one strand is given
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
	  // Does variant fall within a CDS exon?
	  //
	  
// 	  for ( unsigned int s = 0 ; s < r->subregion.size(); s++ )
// 	    {
// 	      std::cout << "REG " << r->subregion[s].meta.get1_int( PLINKSeq::TRANSCRIPT_FRAME()  ) << " "
// 			<< r->subregion[s].coordinate() << " "
// 			<< r->subregion[s].CDS() << " " 
// 			<< r->subregion[s].exon() << " " 
// 			<< r->subregion[s].start_codon() << " " 
// 			<< r->subregion[s].stop_codon() << " "
// 			<< "\n";
		
// 	    }

	  for ( unsigned int s = 0 ; s < r_cds.subregion.size(); s++ )
	    {	      
	      
// 	      std::cout << "subregion = " << s << " of " << r_cds.subregion.size() << " " 
// 			<< r_cds.subregion[s].CDS() << " " 
// 			<< r_cds.subregion[s].exon() << " " 
// 			<< r_cds.subregion[s].start_codon() << " " 
// 			<< r_cds.subregion[s].stop_codon() << " [" 
// 			<< r_cds.subregion[s].meta.get1_int( PLINKSeq::TRANSCRIPT_FRAME() ) << "] "
// 			<< r_cds.subregion[s].coordinate() << "\n";
		

	      if ( bp1 >= r_cds.subregion[s].start.position() && 
		   bp1 <= r_cds.subregion[s].stop.position() ) 
		{
		  
		  CDS_exons.insert(s);
		  in_CDS_exon = s;
		  
		  if ( bp1 - r_cds.subregion[s].start.position() < 3 && s>0 ) 
		    CDS_exons.insert(s-1);
		  
		  if ( r_cds.subregion[s].stop.position() - bp1 < 3 && s < r_cds.subregion.size()-1 )
		    CDS_exons.insert(s+1);
		  		
		}
	    
	    }

	  
	  //
	  // Is this a SPLICE-SITE? 
	  //


	  for ( unsigned int s = 0 ; s < r_exon.subregion.size(); s++ )
	    {
	      
// 	      std::cout << "s = " << s << " " << first_exon << " " << last_exon << " " << r->subregion.size() << "\n";
// 	      std::cout << "bp1 = " << bp1 << "\t" << r->subregion[s].start.position() << " " << r->subregion[s].stop.position() << "\n";
	      
	      bool splice = false;
	      int splicedist = 0;
	      
	      // Intronic -2 would be the second base before the
	      // splice acceptor site (that is, the A of the invariant
	      // AG that ends an intron) 

	      // Intronic + would be after a 5' donor site Synonymous
	      // +2 would be the second base of an exon,

	      //  5'splice site
	      //    Syn -1 would be last base of exon.
	      //    In  +1

	      // 3'
	      //    Syn +1
	      //    In  -1

	      //  --> + ve strand

	      //  --exon1---|        |--exon2---
	      //             GU    AG
	      // 5' splice site     3' splice site


	      
	      //  ---> +ve strand
	      //
	      //   |---------|          |------------|
	      //             5'         3'
	      // +            123       123
	      // -         321       321
	      

	      //                       <--- -ve strand
	      //
	      //   |---------|          |------------|
	      //             3'         5'
	      // +         321       321 
	      // -            123       123
	      


              if ( s != first_exon  &&  abs( r_exon.subregion[s].start.position() - bp1 ) < 3 )
                {
                  if ( negative_strand )
                    {
                      SeqInfo sie = SeqInfo( r->name , ESPLICE3 );
                      SeqInfo si = SeqInfo( r->name , SPLICE3 );
                      si.splicedist = sie.splicedist = r_exon.subregion[s].start.position() - bp1; 
                      if ( si.splicedist <= 0 ) { --si.splicedist; --sie.splicedist; }
                      if ( si.splicedist > -3 && si.splicedist < 0) annot.insert( si );
                      if ( si.splicedist > 0 && si.splicedist < 3) annot.insert( sie );
                    }
                  else
                    {
                      SeqInfo sie = SeqInfo( r->name , ESPLICE5 );
                      SeqInfo si = SeqInfo( r->name , SPLICE5 );
                      si.splicedist = sie.splicedist = bp1 - r_exon.subregion[s].start.position(); 
                      if ( si.splicedist >= 0 ) { ++sie.splicedist; ++si.splicedist; } 
                      if ( si.splicedist < 3 && si.splicedist > 0) annot.insert( si );
                      if ( si.splicedist > -3 && si.splicedist < 0) annot.insert( sie );

                    }
                }
              
              if ( s != last_exon && abs( r_exon.subregion[s].stop.position() - bp1 ) < 3 )
                {
                  if ( negative_strand )
                    {
                      SeqInfo sie = SeqInfo( r->name , ESPLICE5 );
                      SeqInfo si = SeqInfo( r->name , SPLICE5 );
                      si.splicedist = sie.splicedist = r_exon.subregion[s].stop.position() - bp1;
                      if ( si.splicedist >= 0 ) { ++si.splicedist; ++sie.splicedist; }
                      if ( si.splicedist < 3 && si.splicedist > 0 ) annot.insert( si );
                      if ( si.splicedist > -3 && si.splicedist < 0 ) annot.insert( sie );
                    }
                  else
                    {
                      SeqInfo sie = SeqInfo( r->name , ESPLICE3 );
                      SeqInfo si = SeqInfo( r->name , SPLICE3 );
                      si.splicedist = sie.splicedist = bp1 - r_exon.subregion[s].stop.position() ;
                      if ( si.splicedist <= 0 ) { --si.splicedist; --si.splicedist; }
                      if ( si.splicedist > -3 && si.splicedist < 0 ) annot.insert( si );
                      if ( si.splicedist < 3 && si.splicedist > 0 ) annot.insert( sie );
                    }
                }

            }
	
	  
	  
	  //
	  // If no exons attached (but was pulled into this region)
	  // or (splice site)
	  //

	  if ( CDS_exons.size() == 0 ) 
	    {	      	      
	      annot.insert( SeqInfo( r->name , INTRON ) );
	      ++r; // and skip to next region
	      continue;
	    }
	  
	  
	  //
	  // Get reference sequence (CDS)
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
	      ref_cds = getrc( ref_cds );
	      var_allele = getrc( var_allele );
	    }
	  
	  //
	  // Replace position with alternate allele
	  //
	  
	  std::string var_cds = ref_cds; 
	  
	  var_cds.replace( pos_extracted_seq-1, 1, var_allele );
	  
	  
	  
	  //
	  // Are reference and variant sequences identical for this gene?
	  //
	  
	  if ( ref_cds == var_cds ) 
	    {	      
	      annot.insert( SeqInfo( MONO ) );	      	      
	      ++r; // skip to next transcript
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

	  while ( trans_var.size() < longest ) { trans_var += "_"; }
	  while ( trans_ref.size() < longest ) { trans_ref += "_"; }


	  //
	  // Make calls
	  //
	  
	  std::vector<std::string> difs;
	  
	  for ( unsigned int i=0; i< trans_var.size(); i++ )
	    {
	      
	      if ( trans_ref[i] != trans_var[i] )
		{
		  
		  seq_annot_t type = MIS;

		  if      ( trans_var[i] == '*' ) type = NON; 
		  else if ( trans_ref[i] == '*' ) type = RT;
		  
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
  MetaInformation<VarMeta>::field( "_INTRON" , META_TEXT , -1 , "Intronic") ;
  MetaInformation<VarMeta>::field( "_FRAMESHIFT" , META_TEXT , -1 , "Frameshift allele");
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
