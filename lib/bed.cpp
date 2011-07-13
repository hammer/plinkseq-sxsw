#include "bed.h"

#include "helper.h"
#include "vardb.h"
#include "seqdb.h"
#include "inddb.h"

#include <fstream>
#include <bitset>

int BEDReader::read_fam() 
{

  // Read through once to get IDs and phenotype type.
  
  if ( ! Helper::fileExists( fam_filename ) ) return 0;
  
  std::ifstream FAM( fam_filename.c_str() , std::ios::in );

  std::set<std::string> sfid;
  std::set<std::string> siid;
  std::set<std::string> sjoint;
  
  bool binary = true;

  int ni = 0;
  while ( ! FAM.eof() )
    {
      std::string fid, iid, pat, mat, sex, phe;
      FAM >> fid >> iid >> pat >> mat >> sex >> phe;
      if ( fid == "" ) continue;
      ++ni;
      sfid.insert(fid);
      siid.insert(iid);
      sjoint.insert(fid+"_"+iid);
      if ( phe != "0" && phe != "1" && phe != "2" && phe != "-9" )
	binary = false;
    }



  bool use_fid = ni == sfid.size();
  bool use_iid = ni == siid.size();
  bool use_joint = ni == sjoint.size();
  
  if ( ! ( use_fid || use_iid || use_joint ) ) 
    {
      FAM.close();
      return 0;
    }

  //
  // Storing phenotype information?
  //

  int pid = 0;
  

  vardb->begin();

  if ( inddb )
    {
      inddb->begin();

      pid = inddb->insert_phenotype( phenotype_name , 
				     binary ? "Integer" : "Float" , 
				     ".", 
				     "." );
    }

  
  //
  // Read through second time to store actual data
  //
  
  FAM.clear() ;
  FAM.seekg(0, std::ios::beg) ;
  
  while ( ! FAM.eof() )
    {
      
      std::string fid, iid, pat, mat, sex, phe;
      
      FAM >> fid >> iid >> pat >> mat >> sex >> phe;
      
      if ( fid == "" ) continue;
      
      std::string id;
      if ( use_fid ) id = fid;
      else if ( use_iid ) id = iid;
      else id = fid + "_" + iid;
      
      Individual ind(id);
      
      // Register this individual in VARDB
      vardb->insert( file_id , ind );
      
      // Optionally, store phenotypic information in 
      if ( inddb )
	{	  

	  ind.fid( fid );
	  ind.iid( iid );
	  ind.pat( pat );
	  ind.mat( mat );
	  ind.sex( sex );
	  
	  int idx = inddb->insert(ind);

 	  if ( binary ) 
 	    {
 	      int x;
 	      if ( Helper::str2int( phe,x ) )
 		inddb->insert( idx , pid , x );
 	    }	  
 	  else
 	    {
 	      double x;
 	      if ( Helper::str2dbl( phe,x ) )
 		inddb->insert( idx , pid , x );	     
 	    }
	}
    }

  vardb->commit();

  if ( inddb )
    inddb->commit();
  
  FAM.close();
  
  return ni;
}

int BEDReader::read_bim()
{
  if ( ! Helper::fileExists( bim_filename ) ) return 0;
  std::ifstream BIM( bim_filename.c_str() , std::ios::in );
  int cnt = 0;
  while ( ! BIM.eof() )
    {
      BEDLocus l;
      BIM >> l.chr >> l.name >> l.pos >> l.bp >> l.allele1 >> l.allele2 ;
      if ( BIM.eof() ) continue;
      ++cnt;
      locus.push_back( l );
    }
  BIM.close();
  return locus.size();
}

bool BEDReader::read_header( std::ifstream & B )
{
  
  char ch[1];
  B.read(ch,1);
  std::bitset<8> b;
  b = ch[0];
  
  // If v1.00 file format?
  // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
  
  if (   ( b[2] && b[3] && b[5] && b[6] ) && 
	 ! ( b[0] || b[1] || b[4] || b[7] )    )
    {
      B.read(ch,1);
      b = ch[0];   
      if (   ( b[0] && b[1] && b[3] && b[4] ) && 
	     ! ( b[2] || b[5] || b[6] || b[7] )    )
	{
	  // SNP major coding?
	  B.read(ch,1);
	  b = ch[0];   
	  return ( b[0] );
	}      
      return false;      
    }
  return false;
}


bool BEDReader::read_bed()
{



  //
  // Get slot in VARDB
  //
  
  file_id = vardb->insert( bed_filename , ftag );
  

  //
  // Read FAM and BIM files
  //

  const int ni = read_fam();
  const int nv = read_bim();
  
  plog << "expecting " << nv << " variants on " << ni << " individuals\n";

  const long int ng = (long int)ni * (long int)nv;
  if ( ng == 0 ) return false;
  
  //
  // Open BED file (binary PLINK format file)
  //
  
  if ( ! Helper::fileExists( bed_filename ) ) return false;
  
  std::ifstream B( bed_filename.c_str() , std::ios::in | std::ios::binary );

  
  // Track file source
  vardb->insert_header( file_id , "format", "PLINK/BEDv1.00" );

  
  //
  // This BED file must be v1.00 format, and SNP major  
  //
  
  if ( ! read_header( B ) ) 
    {
      B.close();
      return false;
    }


  //
  // Read genotypes, and save into VARDB
  //

  long int cnt = 0;

  vardb->begin();

  for ( int s = 0 ; s < locus.size(); s++ )
    {
      
      Variant v;
      
      //
      // Assign MAP information
      //
      
      v.chromosome( locus[s].chr );
      v.position( locus[s].bp );
      v.name( locus[s].name );

      // Use SEQDB to look up reference allele. Implies that all
      // alleles are known to be on the positive strand;
      
      bool ref1 = true;  // A1 is reference allele
      
      if ( locus[s].allele1 == "0" || locus[s].allele1 == "X" )
	{
	  std::string a2 = locus[s].allele2;
	  locus[s].allele2 = locus[s].allele1;
	  locus[s].allele1 = a2;
	}
      
      if ( seqdb )
	{
	  std::string ref = seqdb->lookup( locus[s].chr , locus[s].bp );
	  Helper::str2upper(ref);
	  if ( ref != "N" ) 
	    {
	      if ( ref == locus[s].allele2 ) ref1 = false;
	      else if ( ref != locus[s].allele1 ) 
		plog.warn( "mismatching reference allele in BED versus SEQDB" , 
			   v.displaycore() + " " + locus[s].allele1+"/"+locus[s].allele2 + " vs " + ref );
	    }
	}
      
      if ( ref1 ) 
	{
	  v.consensus.reference( locus[s].allele1 );
	  v.consensus.alternate( locus[s].allele2 );	  
	}
      else
	{
	  v.consensus.reference( locus[s].allele2 );
	  v.consensus.alternate( locus[s].allele1 );	  
	}
      
      
      if ( s % 1000 == 0  )
	plog.counter( "inserted " + Helper::int2str( s ) + " variants" );
      
      

      //
      // Get genotypes
      //

      v.resize( ni );
      
      int indx = 0;
      
      while ( indx < ni )
	{
	  
	  std::bitset<8> b;
	  
	  char ch[1];
	  B.read(ch,1);
	  
	  if (!B) 
	    Helper::halt("Problem with the BED file...has the FAM/BIM file been changed?\n");          
	  
	  b = ch[0];      

	  // Parse the four genotypes in this byte

	  int c=0;
	  
	  while (c<7 && indx < ni ) 
	    {
	      
	      bool g1 = b[ c++ ];
	      bool g2 = b[ c++ ];
	      
	      // Missing code? 
	      if ( g1 && !g2 ) 
		v(indx).null( true );
	      else
		{
		  int c = ref1 ? (int)g1 + (int)g2 : 2 - (int)g1 - (int)g2 ;
		  v(indx).set_alternate_allele_count( c );
		}
	      
	      ++cnt;
	      ++indx;
	    } // next individual in byte

	} // next byte/individual     
      

      //
      // Put this variant into the database
      //
      
      vardb->insert_consensus( file_id , v );      
      
      
    } // next SNP
 
  B.close();
  
  
  //
  // Ensure variant index is set on VARDB
  //
  
  vardb->commit();

  vardb->index();

  //
  // Did we read the correct number of genotypes? 
  //

  int2 niv = vardb->make_summary( file_id ) ;


  plog << "inserted " << niv.p1 << " individuals, " << niv.p2 << " variants\n";

  return cnt == ng; 

}

