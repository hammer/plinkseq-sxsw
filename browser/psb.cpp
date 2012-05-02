#include "psb.h"
#include "cgi.h"
#include "lines.h"

#include <iostream>
#include <algorithm>

using namespace ExomeBrowser;


int main()
{

  
  const bool cgi = true;

  
  //
  // Get CGI variables from POST
  //

  char **cgivars ;  int i ;
  
  if ( cgi ) 
    cgivars = getcgivars() ;
  
  
  //
  // Define values
  //
  
  std::string project_path = "";
  std::string pwd = "(if required)";
  std::string loc_set = "refseq";
  std::string genename = "";
  QType q = Q_ERROR;
  std::string var_value = "";
  std::string ind_value = "";
  std::string chr_code = "";
  int chr = 0, bp1 = 0, bp2 = 0;
  std::string pheno = "";
  bool from_top = false;
  bool gview = false;
  bool index = false;

  // Auxilliary information to send to variant printing functions

  Aux a;
  
  
  if ( cgi ) 
    {
      for (i=0; cgivars[i]; i+= 2)
	{
	  
	  std::string str = cgivars[i];
	  
	  if ( str == "q" )
	    {
	      std::string s = cgivars[i+1];
	      
	      if ( s == "v" ) q = Q_VARIANT;
	      else if ( s == "i" ) q = Q_INDIV;
	      else if ( s == "r" ) q = Q_REGION;
	      else if ( s == "gview" ) { gview = true; q = Q_GRAPHICAL_VIEW; } 
	      else if ( s == "varsetlist" ) q = Q_VARSETLIST;
	      else if ( s == "pgrid" )  q = Q_PHENOGRID;
	      else if ( s == "indlist" ) q = Q_INDGRID;
	      else if ( s == "glist" ) q = Q_GENELIST;
	      else if ( s == "mflist" ) q = Q_METALIST;
	      else if ( s == "plist" ) q = Q_PHELIST;
	      else if ( s == "lslist" ) q = Q_LOCSETLIST;
	      else if ( s == "psummary" ) q = Q_PROJSUMMARY;
	      else q = Q_ERROR;	  
	    }
	  
	  if ( str == "addannot" )
	    a.add_annot = true;
	  
	  if ( str == "getgene" ) 
	    from_top = true;
	  
	  if ( str == "gview" ) 
	    gview = true;

	  if ( str == "indgrid" ) 
	    index = true;
	  
	  if ( str == "pgrid" ) 
	    q = Q_PHENOGRID;

	  if ( str == "val" ) 
	    var_value = cgivars[i+1];
	  
	  if ( str == "ind" ) 
	    ind_value = cgivars[i+1];
	  
	  if ( str == "proj" ) 
	    {
	      project_path = cgivars[i+1];
	      a.add_form_value("proj", project_path );
	    }
	  
	  if ( str == "passwd" )
	    {
	      pwd = cgivars[i+1];	      
	    }


	  if ( str == "loc" ) loc_set = cgivars[i+1];

	  
	  if ( str == "indiv_list" ) 
	    {
	      a.indiv_list_url = cgivars[i+1];
	      a.indiv_list_vec = Helper::parse( a.indiv_list_url , " ," );
	      for (int i=0;i<a.indiv_list_vec.size();i++)
		a.indiv_list.insert( a.indiv_list_vec[i] );
	    }
	  

	  if ( str == "varset" ) 
	    {
	      a.varset_url = cgivars[i+1];
	      a.varset = Helper::parse( a.varset_url , ", " );
	      for (int i=0;i<a.varset.size();i++) a.varset_set.insert( a.varset[i] );
	    }
	  
	  
	  if ( str == "regs" )
	    {
	      a.reg_list = a.reg_list_url = cgivars[i+1];
	      while ( 1 ) 
		{
		  int p = a.reg_list.find( "," );
		  if ( p == std::string::npos ) break;
		  a.reg_list.replace( p , 1 , "\n" );
		}
	      std::replace( a.reg_list_url.begin(), a.reg_list_url.end() , '\r' , ' ' );
	      std::replace( a.reg_list_url.begin(), a.reg_list_url.end() , '\n' , ',' );
	    }


	  if ( str == "masks" )
	    a.msk = Helper::parse( cgivars[i+1] , " " );

	  
	  if ( str == "inc_fltr" )
	    a.inc_fltr = cgivars[i+1] ;


	  if ( str == "vinc_fltr" )
	    a.vinc_fltr = cgivars[i+1];


	  if ( str == "meta" ) 
	    {
	      a.mf = Helper::parse( cgivars[i+1] , " ," );
	    }


	  if ( str == "ref_append" ) 
	    {
	      a.ref_append_url = cgivars[i+1];
	      a.ref_append = Helper::parse( cgivars[i+1] , " ," );
	    }


	  if ( str == "loc_append" ) 
	    {
	      a.loc_append_url = cgivars[i+1];
	      a.loc_append = Helper::parse( cgivars[i+1] , " ," );
	    }


	  if ( str == "pheno" )
	    pheno = cgivars[i+1];
	  
	}
      
      
      // ensure correct mode
      if      ( from_top ) q = Q_REGION;
      else if ( gview )    q = Q_GRAPHICAL_VIEW;
      else if ( index )    q = Q_INDGRID;
	
    }
          
  
  
  // ** Free anything that needs to be freed **/
  
  if ( cgi )
    {
      for (i=0; cgivars[i]; i++) free(cgivars[i]) ;
      free(cgivars) ;
    }


  //
  // The 'index' Q_INDGRID uses Q_REGION's code
  //

  if ( q == Q_INDGRID )
    {
      a.indiv_genogrid = true;
      q = Q_REGION;
    }


  //
  // Start HTML form and any output (unless we are deferring this 
  // because we need to populate the functions for a graphical view, 
  // and so this will be called later
  //
  
  if ( q != Q_GRAPHICAL_VIEW ) 
    {
      write_html_header();
      std::cout << "<body>";
    }
  

  //
  // Set up project
  //


  GStore g;
  
  a.g = &g;
  a.loc_set = loc_set;
  g.set_project( project_path ) ;
  
  
  
  //
  // Start page 
  //
  
  if ( q != Q_GRAPHICAL_VIEW ) 
    write_start_page(g,loc_set,q,a,pheno,pwd,project_path);
  


  ////////////////////////////////////////////////////////////////////////////
  

  //
  // Set up basic stuff to perform query
  //


  if ( ! g.pwd(pwd) ) 
    Helper::halt("<b>access denied: password does not match</b>");
  

  //
  // Initial Mask objects
  //

  std::string mstr;
  for (int i = 0 ; i < a.msk.size(); i++)
    mstr += ( i ? " " : "" ) + a.msk[i];


  // Add some additional mask items to speed things up, and 
  // check that the output isn't too large; note: the 'no-geno'
  // assumes that complex mask queries that use genotype data 
  // aren't applied, e.g. 
  // include=" g( DP >= 10 ) > 0.1 "
  
  if ( q == Q_REGION ) 
    {
      // if no mask, safe to add 'no-geno', otherwise, we should keep in case
      // and no dichot phenotype (i.e. for C/C counts)
      
      if ( mstr == "" && pheno == "" && a.inc_fltr == "" && a.vinc_fltr == "" && ! a.indiv_genogrid ) 
	mstr = "limit=1000 no-geno " + mstr;
      else
	mstr = "limit=1000 " + mstr;      
    }
  
  if ( a.vinc_fltr != "" ) 
    { 
      mstr += " v-include=\"" + a.vinc_fltr + "\"";
    }
  
  
  //
  // Create mask
  //
  
  Mask m( mstr , a.inc_fltr , a.inc_fltr != "" );

  
  //
  // Add an individual inclusions/exclusions
  //

  if ( a.indiv_list.size() > 0 ) 
    {
      m.include_indiv( a.indiv_list_vec );
    }
  
  
  //
  // Create/populate individual-map
  //

  g.indmap.populate( g.vardb, g.phmap, m );

  
  //
  // Sanity check
  //

  if ( g.indmap.size() > 50 && a.indiv_genogrid ) 
    {
      std::cout << "Genotype grid mode designed for small (sub)samples (50 or fewer individuals)"
		<< "</BODY></HTML>";
      exit(0);
    }


  //
  // Expand wildcard?
  //

  std::set<std::string> onecopy;
  for (int m=0; m<a.mf.size(); m++)
    {
      if ( a.mf[m] == "*" )
	{
	  a.mf.erase( a.mf.begin() + m );
	  std::map<int,std::string> f = g.vardb.fetch_files();
	  std::map<int,std::string>::iterator i = f.begin();
	  while ( i != f.end() )
	    {
	      std::vector<std::map<std::string,std::string> > m0 = g.vardb.fetch_metatypes( i->first );
	      for (int j=0; j<m0.size(); j++)
		{
		  if ( m0[j]["GRP"] == "Variant" )
		    {
		      std::string mval = m0[j]["NAME"];
		      if ( onecopy.find( mval ) == onecopy.end() )
			{
			  onecopy.insert( mval );
			  a.mf.push_back( mval );
			}
		    }
		}
	      ++i;
	    }
	  break;
	}
    }
  
  
  // check for any non-pp fields
  // i.e. that have '+' suffix
  
  for (int i=0; i<a.mf.size();i++)
    {
      if ( a.mf[i][a.mf[i].size()-1] == '+' )
	{
	  a.mf[i] = a.mf[i].substr(0,a.mf[i].size()-1);
	  a.mfpp[ a.mf[i] ] = 1;
	}
    }




  //
  // Has a phenotype been specified? If so, attach
  //

  // Currently, we assume phenotype corresponds to a 
  // standard (2/1) -> (case/control) variable

  a.show_phenotype = false;

  if ( pheno != "" )
    {
      
      int p_id = g.inddb.fetch_pheno_id( pheno );
      if ( p_id == 0 ) 
	{
	  a.show_phenotype = false;
	}
      else
	{
	  a.show_phenotype = true;
	  a.phenotype_name = pheno;
	  g.phmap.set_phenotype( pheno );
	}
    }
  

  // 
  // Process region list; this only works for initial gene-view
  //
  
  a.genes.clear();
  a.regions.clear();

  // collect all transcript ID here

  std::vector<std::string> tnames; // transcript IDs
  std::vector<std::string> snpids; // SNP IDs
  
  
  std::vector<std::string> tok = Helper::whitespace( a.reg_list );
  
  for (int i=0; i<tok.size(); i++)
    {      
      
      bool okay = true;      
      Region r(tok[i],okay);
      if ( okay ) 
	{
	  a.regions.push_back( r );	  
	}
      else
	{
	  // look for as a gene (assuming upper case for all IDs)
	  std::string tmp = tok[i];
	

	  // Does this look like a gene name? 
	  std::string genename = tmp;
	  Helper::str2upper( genename );
	
	  
	  std::set<std::string> trans_names = g.locdb.targetted_lookup_alias( genename , 
									      ExomeBrowser::symbol, 
									      loc_set );
	  
	  // If no matches, assume the original input was in
	  // transcript form, e.g. NA_12345

	  if ( trans_names.size() == 0 ) 
	    {
	      
	      // this could be a SNP ID or a gene ID. Search for both	      
	      snpids.push_back( tmp );
	      
	      // and search for a gene transcript
	      tnames.push_back( tmp );
	    }
	  else 
	    {
	      
	      // ... otherwise add what we've found from the lookup instead
	      
	      std::set<std::string>::iterator ii = trans_names.begin();
	      while ( ii != trans_names.end() ) 
		{
		  tnames.push_back( *ii ) ; 
		  ++ii;
		}
	    }	  
	}
    }
  

  //
  // Add any rs-ID searchs
  //

  if ( snpids.size() > 0 ) 
    m.include_id( snpids );


  //
  // Variant sets
  //

  if ( a.varset.size() > 0 ) 
    {
      for (int i=0;i<a.varset.size(); i++)
	{
	  if ( a.varset[i][0] == '-' )
	    m.exclude_varset( a.varset[i].substr(1) );
	  else
	    m.require_varset( a.varset[i] );
	}
    }


  // 
  // Reference/Locus filters (if in non GUI mode)
  // 

  //  dbsnp    require dbsnp membership
  //  +dbsnp   append this track
  //  -dbsnp   exclude this track
  
  if ( q == Q_REGION ) 
    {
      for (int i=0; i<a.ref_append.size(); i++) 
	{
	  if ( a.ref_append[i][0] == '-' ) 
	    m.exclude_ref( a.ref_append[i].substr(1) );
	  else if ( a.ref_append[i][0] == '+' ) 
	    m.append_ref( a.ref_append[i].substr(1) );
	  else
	    m.require_ref( a.ref_append[i] );
	}

      for (int i=0; i<a.loc_append.size(); i++) 
	{
	  if ( a.loc_append[i][0] == '-' )
	    m.exclude_loc( a.loc_append[i].substr(1) );
	  else if ( a.loc_append[i][0] == '+' ) 
	    m.append_loc( a.ref_append[i].substr(1) );	  
	  else
	    m.require_loc( a.loc_append[i] );
	}
    }
  

  //
  // Get genomic loci for transcripts
  //
  
  std::vector<Region> trans = g.locdb.fetch( loc_set , tnames );


  //
  // Will we be reporting on 1, or on multi transcripts (e.g. impacts display of exon #)
  //
  
  if ( a.regions.size() == 0 && trans.size() == 1 ) 
    a.single_transcript = true;
  

  
  //
  // For a gene/region-based query, some verbose output regarding transcript info, etc
  //
  
  if ( q == Q_REGION ) 
    {
      
      if ( trans.size() > 0 ) 
	std::cout << "Found " << trans.size() << " matching transcript(s)</b></p>";
      
      std::cout << "<pre><font size=-1>";
      
      for (int r=0; r<trans.size(); r++)
	{
	  
	  //	  std::cout << Helper::sw( trans[r].altname , 15 ) << "  ";
	  
	  std::cout << Helper::sw( a.getURL()->addField("q", "r")
				   ->addField("regs", trans[r].name)
				   ->printLink(trans[r].name) , 15 ) << " ";
	  
	  // attempt to get symbol given gene-name

	  std::set<std::string> gname = g.locdb.targetted_lookup_alias( trans[r].altname , loc_set, ExomeBrowser::symbol );
	  
	  std::set<std::string>::iterator ii = gname.begin();
	  while ( ii != gname.end() )
	    {
	      if ( ii != gname.begin() ) std::cout << ",";
	      std::cout << *ii;	      
	      ++ii;
	    }

	  std::cout << "  ( " << Helper::chrCode(trans[r].start.chromosome()) << ":" 
		    << trans[r].start.position() << ".."
		    << trans[r].stop.position() ;
	  std::cout <<" )<br>";
	}
      
      std::cout << "</font></pre>";
      
      
      //
      // Now, 'genename' will be original search term
      //
      
      if ( a.single_transcript )
	{
	  
	  genename = tnames[0];

	  Region reg = g.locdb.get_region( loc_set , genename ) ;
	  
	  if ( reg.start.chromosome() == 0 )
	    {
	      std::cout << "Could not find gene <b> " << genename << "</b> in <b>" << loc_set << "</b> list "
			<< "</BODY></HTML>";
	      exit(0);
	    }
	  
	  a.region = reg;
	  std::string g_chr_code = Helper::chrCode( reg.start.chromosome() );
	  int g_chr = reg.start.chromosome();
	  int g_bp1 = reg.start.position();
	  int g_bp2 = reg.stop.position();
	  int n_exons = reg.subregion.size();
	  
	  // Query to get all other transcripts in region
	  
	  std::set<Region> others = g.locdb.get_regions( g.locdb.lookup_group_id( loc_set ) , 
							 g_chr , g_bp1 , g_bp2 );
	
	  std::cout << "<b>" << genename << "</b> location : " 
		    << g_chr_code << ":" << g_bp1 << ".." << g_bp2 
		    << "  (view in the UCSC genome browser: ";

	  std::cout << "<a href=\"http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg18&position=" 
		    << g_chr_code << ":" << g_bp1 << ".." << g_bp2     
		    << "\" target=\"_blank\">hg18</a> ";

	    
	  std::cout << "<a href=\"http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg19&position=" 
		    << g_chr_code << ":" << g_bp1 << ".." << g_bp2     
		    << "\" target=\"_blank\">hg19</a>)";

	  std::cout << "</p>";
	  
	  
	  std::cout << "<pre><font size=-1>";
	  
	  double genomic_kb = (g_bp2 - g_bp1 + 1 ) / 1000.0;
	  double coding_kb = 0;
	  for (int s=0; s<reg.subregion.size(); s++)
	    coding_kb += ( reg.subregion[s].stop.position() - reg.subregion[s].start.position() + 1 );
	  coding_kb /= 1000.0;
	  
	  std::cout << "   Genomic kb = " << genomic_kb 
	       << "; coding kb = " << coding_kb 
	       << " in " << n_exons << " exons</p>";
	  
	  // List exons, and other transcripts that overlap this one
      
	  for (int s=0; s<reg.subregion.size(); s++)
	    {
	      
	      std::cout << "   Exon" << Helper::sw(s+1,3) ;
	      std::cout << Helper::sw( Helper::int2str(reg.subregion[s].start.position()) 
				       + " .. " + Helper::int2str(reg.subregion[s].stop.position()) , 16 );
	      std::cout << " | " ;
	      std::cout << Helper::sw( ( reg.subregion[s].stop.position() 
					 - reg.subregion[s].start.position() + 1 )  , 5 );
	      std::cout << "bp | ";
	      
	      // other exons this exon overlaps with? 
	      std::set<Region>::iterator t = others.begin();
	  
	      while ( t != others.end() )
		{
		  if ( t->name != reg.name )
		    for (int u=0; u< t->subregion.size(); u++)
		      {
			if ( reg.subregion[s].overlaps( t->subregion[u] ) )
			  {
			    std::cout << a.getURL()->addField("q", "r")
			      ->addField("regs", t->name)
			      ->printLink(t->name + "(" + t->altname + ")") << " ";
			  }
		      }
		  ++t;
		}
	      std::cout << "<br>";
	    }
	  
	  std::cout << "</font></pre>";
      
	}

      if ( trans.size() > 0 ) 
	std::cout << "<hr>";      


      // report regions
  
      
      if ( a.regions.size() > 0 ) 
	{

	  std::cout << "<pre><b>Regions:</b><br>";
	  for (int i=0; i<a.regions.size(); i++)
	    std::cout << "   " << Helper::chrCode( a.regions[i].chromosome()) << ":"
		      << a.regions[i].start.position() << ".."
		      << a.regions[i].stop.position() << "<br>";
	  std::cout << "</pre>";
	}
      


      if ( a.genes.size() > 0 ) 
	{	  

	  std::vector<std::string> cc = a.genes;
	  a.genes.clear();
	  
	  std::cout << "<pre><b>Genes:</b><br>";
	  for (int i=0; i< cc.size(); i++)
	    {
	      
	      std::set<std::string> trans = 
		g.locdb.targetted_lookup_alias( cc[i] , ExomeBrowser::symbol , loc_set ) ;   
	      
	      // if no aliases match, assume this is a refseq transcript name
	      if ( trans.size() == 0 ) 
		{
		  a.genes.push_back( cc[i] );
		  std::cout << "   " << Helper::sw( cc[i] , 10 ) << "<br>";
		}
	      else // add all matching aliases
		{
		  std::cout << "   " << Helper::sw( cc[i] , 10 ) << " : ";
		  
		  std::set<std::string>::iterator ii = trans.begin();
		  while ( ii != trans.end() ) 
		    {
		      a.genes.push_back( *ii );
		      std::cout << *ii << " ";
		      ++ii;
		    }		 
		  std::cout << "<br>";
		}
	    }
	  std::cout << "</pre>";
	}
      
      if ( a.regions.size() + a.genes.size() > 0 ) 
	std::cout << "<hr>";
    }
  
  
  //
  // In gene, or regional mode
  //
  
  
  if ( q == Q_REGION ) 
    {
      
      
      a.table_row.clear();
      
      a.vcnt = 0;
      
      if ( trans.size() > 0 ) 
	{
	  
	  m.include_loc( loc_set );
	  for (int i=0; i<trans.size(); i++)
	    {
	      m.subset_loc( loc_set , trans[i].name );
	    }
	}


      // If any regions specified, add as a requirement
      
      if ( a.regions.size() > 0 )
	{	  

	  a.region_search = true;

	  for (int r = 0 ; r < a.regions.size() ; r++) 
	    m.include_reg( a.regions[r] ); 

	}


      //
      // Append gene names if one or more region specified
      //
      
      if ( a.region_search )
	m.append_loc( loc_set );
      
      // Other genes specified?

      for (int r = 0 ; r < a.genes.size() ; r++) 
	m.subset_loc( loc_set , a.genes[r] );	    
	
      
      // Get actual information from VARDB
      
      IterationReport irep = g.vardb.iterate( f_display , &a, m );

      if ( irep.reached_limit() ) 
	std::cout << "<b>Reached limit on number of variants that can be returned, " << irep.processed() << "</b><br>";
      else
	std::cout << "Found " << irep.processed() << " variants that match query<br>";

      // Now table should be sorted in order, with duplicates removed, in aux structure;

      // string headers;
      // map<int,std::string> table_row;
      
      std::cout << a.headers << "</p>";
      // LINKS CHR POS NAME ALT/REF FileID QUAL INFO(filter) SAMPLE_CNT VMETA(compressed)
      
      std::cout << "<table border=1>"
		<< "<tr><th>#</th><th>Indiv</th><th>Chr</th><th>Pos</th>";
      
      if ( a.single_transcript )
	std::cout << "<th>Exon</th>";
      
      std::cout << "<th>ID</th><th>Ref/Alt</th>"
		<< "<th># samples</th>"      
		<< "<th>Filter</th>";      
      
      
      // Optional case/control counts?
      
      if ( a.show_phenotype && ( g.phmap.type() == PHE_DICHOT || g.phmap.type() == PHE_FACTOR ) )
	std::cout << "<th>" << a.phenotype_name << "</th>";
      
  
      // Optional meta-fields?
      
      for (int m=0; m< a.mf.size(); m++)
	std::cout << "<th>" << a.mf[m] << "</th>";

      if ( a.region_search ) 
	std::cout << "<th>Locus</th>";
      
      // Optional genotype-grid?
      
      if ( a.indiv_genogrid ) 
	{
	  const int n = g.indmap.size();
	  for (int i=0;i<n;i++) std::cout << "<th align=\"center\">" << g.indmap(i)->id() << "</th>";
	}

     
      // End of header row
      
      std::cout << "</tr>";
      

      // in 'index' mode, show phenotypes

      if ( a.indiv_genogrid )
	{
	  
	  std::cout << "<tr><td>.</td><td>.</td><td>.</td><td>.</td>";

	  if ( a.single_transcript )
	    std::cout << "<td>.</td>";

	  std::cout << "<td>.</td><td>.</td><td>.</td><td>.</td>";
	  
	  if ( a.region_search ) 
	    std::cout << "<th>.</th>";

	  if ( a.show_phenotype && ( g.phmap.type() == PHE_DICHOT || g.phmap.type() == PHE_FACTOR ) )
	    std::cout << "<td>.</td>";

	  for (int m=0; m< a.mf.size(); m++)
	    std::cout << "<td>.</td>";
	  
      	  const int n = g.indmap.size();
	  for (int i=0;i<n;i++) 
	    {
	      std::cout << "<td align=\"center\">";
	      std::stringstream ss;
	      ss << g.indmap.ind(i)->meta;
	      std::vector<std::string> c = Helper::parse( ss.str() , ";" );
	      for (int k=0;k<c.size(); k++) std::cout << c[k] << "<br>";
	      std::cout << "</td>";      
	    } // next indiv
	  std::cout << "</tr>";
	}


      // Rows

      std::map<int,std::string>::iterator k= a.table_row.begin();
      while ( k != a.table_row.end() )
	{
	  std::cout << k->second;
	  ++k;
	}
      std::cout << "</table>";

    }


  //
  // Query for a single variant: no need to call any database iteraton functions; no masks applied here either
  //
  
  if ( q == Q_VARIANT )
    {
      
      std::cout << "Back to "
                << a.getURL()->addField("q", "r")->addField("regs", a.reg_list_url )->printLink("regional variant report");
      std::cout << "<br><hr>";      

       // Expect here we have variant ID
       
       // note == first param is actually redundant      
       
       bool okay = true;
       Region r( var_value , okay  );
       if ( ! okay ) 
	 {
	   std::cout << "Problem processing region code [ " << var_value << " ]</p>";
	   std::cout << "</body></html>";
	   exit(0);
	 }



       //
       // Get variant, based on physical position
       //
       
       Variant var = g.vardb.fetch( r.chromosome() , r.start.position() );

       
       //
       // Display
       //


       std::cout << "<table width=100%><tr><td width=50% valign=top>";

       std::cout << "<h3><font color=\"blue\">Variant information</font></h3>";
         
       std::cout << "<table border=1><tr><th>Field</th><th>Value</th></tr>";
       std::cout << "<tr><td>Name</td><td>" << rs_link( var.name() ) << "</td></tr>";
       std::cout << "<tr><td>Chromosome</td><td>" << Helper::chrCode( var.chromosome() ) << "</td></tr>";
       std::cout << "<tr><td>Position</td><td>" << var.position() << "</td></tr>";
       std::cout << "<tr><td>Reference allele</td><td>" << var.pp_reference() << "</td></tr>";
       std::cout << "<tr><td>Alternate allele(s)</td><td>" << var.pp_alternate() << "</td></tr>";
       std::cout << "<tr><td>Samples</td><td>" << var.n_samples() << "</td></tr>";

       std::vector<std::string> keys = var.meta.keys();
       for (int m=0; m<keys.size(); m++)
	 std::cout << "<tr><td>" << keys[m] << "</td><td>" << var.meta.print( keys[m] ) << "</td></tr>";	  
       std::cout << "</table>";

       std::cout << "</p>";
       

       // Consensus meta-information
       
       SampleVariant & sample = var.consensus;
       
       std::cout << "<p>Consensus</p>";	   
       std::cout << "<table border=1><tr><th>Field</th><th>Value</th></tr>";	   
       std::cout << "<tr><td>Reference allele</td><td>" << sample.pp_reference() << "</td></tr>";
       std::cout << "<tr><td>Alternate allele(s)</td><td>" << sample.pp_alternate() << "</td></tr>";

       if ( sample.quality() < 0 ) 
	 std::cout << "<tr><td>Quality</td><td>" << "NA" << "</td></tr>";
       else
	 std::cout << "<tr><td>Quality</td><td>" << sample.quality() << "</td></tr>";

       std::cout << "<tr><td>Filter</td><td>" << sample.filter() << "</td></tr>";

       std::vector<std::string> keys1 = sample.meta.keys();
       for (int m=0; m< keys1.size(); m++)
	 std::cout << "<tr><td>" << keys1[m] 
	      << "</td><td>" 
	      << sample.meta.print( keys1[m] ) 
	      << "</td></tr>";
       
       std::cout << "</table>";
       
       
       // SampleVariant Tables
       
       const int ns = var.n_samples();
       
       for (int s = 0 ; s < ns ; s++ )
	 {

	   SampleVariant & sample = var.sample( s );
	   
	   std::cout << "<p>Sample " << g.vardb.file_tag( sample.fileset() ) << "</p>";
	   
	   std::cout << "<table border=1><tr><th>Field</th><th>Value</th></tr>";
	   
	   std::cout << "<tr><td>Reference allele</td><td>" << sample.pp_reference() << "</td></tr>";
	   std::cout << "<tr><td>Alternate allele(s)</td><td>" << sample.pp_alternate() << "</td></tr>";

	   if ( sample.quality() < 0 ) 
	     std::cout << "<tr><td>Quality</td><td>" << "NA" << "</td></tr>";
	   else
	     std::cout << "<tr><td>Quality</td><td>" << sample.quality() << "</td></tr>";
	   
	   std::cout << "<tr><td>Filter</td><td>" << sample.filter() << "</td></tr>";
	   
	   std::vector<std::string> keys = sample.meta.keys();
	   for (int m=0; m< keys.size(); m++)
	     std::cout << "<tr><td>" << keys[m] << "</td><td>" << sample.meta.print( keys[m] ) << "</td></tr>";	  
	   
	   std::cout << "</table>";
	   
	 }
       

       //
       // Genotype information
       //
       
       std::cout << "</td><td valign=top>";
       
       std::cout << "<h3><font color=\"blue\">Individual genotypes</font></h3>";
       
       
       // Get list of IDs for consensus set
       std::vector<std::string> id_list = var.ind_id();
       
       std::set<std::string> gmeta;
       std::set<std::string> gRef;
       std::set<std::string> gHet;
       std::set<std::string> gHom;
       std::set<std::string> gMis;

       SampleVariant & con = var.consensus;
       
       if ( id_list.size() != var.size() ) 
	 Helper::halt("internal problem in indmap/var-con sizes");
       
       for (int i=0; i < id_list.size(); i++)
	 {
	   std::vector<std::string> k = var(i).meta.keys();
	   for (int j=0; j<k.size(); j++) gmeta.insert(k[j]);
	   Genotype & g = var(i);
	   if ( g.null() ) gMis.insert( id_list[i] ) ;
	   else 
	     {
	       int ac = g.allele_count( );
	       if ( ac == 0 ) gRef.insert( id_list[i] );
	       else if ( ac == 1 ) gHet.insert( id_list[i] );
	       else if ( ac == 2 ) gHom.insert( id_list[i] );
	     }	  
	 }


       std::cout << "<table border=1><tr><th>Individual ID</th>";
       
       std::cout << "<th nowrap>Sample #</th>";
       
      if ( a.show_phenotype ) 
	std::cout << "<th>Phenotype</th>";
            
      std::cout << "<th>Genotype</th>";
      
      if ( var.flat() )
	{
	  std::set<std::string>::iterator i = gmeta.begin();
	  while ( i != gmeta.end() )
	    {
	      std::cout << "<th>" << *i << "</th>";
	      ++i;
	    }
	}
      else
	std::cout << "<th>Sample/genotype meta-information</th>";
      
      std::cout << "</tr>";


      // Done header, now list individuals, in genotypic order
      
      for (int gt=0; gt<4; gt++)
	{

	  for (int i=0; i<id_list.size(); i++)
	    {
	      
	      
	      if      ( gt == 3 && ! var(i).null() )  continue;
	      else if ( gt != 3 &&   var(i).null() )  continue;
	      else if ( gt != 3 && var(i).allele_count( ) != 2-gt ) continue;
	      
	
	      // Include link to individual-report

	      std::cout << "<tr><td>" 
	                << a.getURL()->addField("q", "i")
		->addField("regs", a.reg_list_url )
		->addField("ind", id_list[i])
		->printLink(id_list[i])
	                << "</td>";
	      
	      
	      // File-set for this sample variant
	      
	      std::cout << "<td>" << var.sample_label(i) << "</td>";
	      

	      // Phenotype (optionally)

	      if ( a.show_phenotype ) 
		{

		  Individual * person = var.ind( i );
		  		  
		  if ( person ) 
		    {
		      
		      if ( g.phmap.type() == PHE_DICHOT )
			{
			  if ( person->affected() == CASE )
			    std::cout << "<td>CASE</td>";
			  else if ( person->affected() == CONTROL )
			    std::cout << "<td>CONTROL</td>";
			  else
			    std::cout << "<td>MISSING</td>";
			}
		      else if ( g.phmap.type() == PHE_QT ) 
			{
			  if ( ! person->missing() )			  
			    std::cout << "<td>" << person->qt() << "</td>";			  
			  else
			    std::cout << "<td>MISSING</td>";			  
			}
		      else if ( g.phmap.type() == PHE_FACTOR ) 
			{
			  if ( ! person->missing() )			  
			    std::cout << "<td>" << person->group_label() << "</td>";			  
			  else
			    std::cout << "<td>MISSING</td>";			  
			}
		    }		  
		  else
		    std::cout << "<td>n/a</td>";
		}

	      if ( gt == 0 ) 
		std::cout << "<td bgcolor=\"red\">";
	      else if ( gt == 1 ) 
		std::cout << "<td bgcolor=\"yellow\">";
	      else if ( gt == 2 )
		std::cout << "<td bgcolor=\"lightgreen\">";
	      else 
		std::cout << "<td bgcolor=\"gray\">";
	      
	      if ( var.flat() )
		{		  
		  std::cout << var.geno_label( var(i) ) ;		  
		}
	      else
		std::cout << var.label(i) ;
	      
	      std::cout << "</td>";
	      
	      if ( var.flat() )
		{
		  std::set<std::string>::iterator k = gmeta.begin();
		  while ( k != gmeta.end() )
		    {
		      if ( var(i).meta.has_field( *k ) )
			{
			  meta_name_t t = *k;
			  std::string value = var(i).meta.print( t );
			  std::cout << "<td>" << ( ( value == "" || value == " " ) ? "&nbsp;" : value ) << "</td>";
			}
		      else
			std::cout << "<td>.</td>";			
		      ++k;
		    }
		}
	      else
		std::cout << "<td>" << var.gmeta_label(i) << "</td>";
	      
	      std::cout << "</tr>";	  
	    }
	}
      
      std::cout << "</table>";
             


      //
      // End of table
      //
      

      std::cout << "</td></tr></table>";
      
    }



  //
  // Individua report
  //

  if ( q == Q_INDIV )
    {

      a.genename = genename;

      std::cout << "Back to "
                << a.getURL()->addField("q", "r")->addField("regs", a.reg_list_url )->printLink("regional variant report");
      
      std::cout << "<hr>";
      

      // 1) Report all phenotypic information for this individual
      
      if ( a.show_phenotype ) 
	{

	  // A single individual ID
	  
	  Individual * person = g.indmap.ind( ind_value );
	  
	  if ( ! person )
	    {
	      std::cout << "Individual [ " << ind_value << " ] not found...</p>";
	      std::cout << "</BODY></HTML>";
	      exit(0);
	    }
	  
	  std::cout << "Phenotype for variable <em>" << pheno << "</em> is [ " ;
	  if ( g.phmap.type() == PHE_DICHOT )
	    {
	      if ( person->affected() == CASE )
		std::cout << "CASE";
	      else if ( person->affected() == CONTROL )
		std::cout << "CONTROL";
	      else
		std::cout << "MISSING";
	    }
	  else if ( g.phmap.type() == PHE_QT ) 
	    {
	      if ( ! person->missing() ) std::cout << person->qt() ;
	      else std::cout << "MISSING";
	    }
	  else if ( g.phmap.type() == PHE_FACTOR )
	    {
	      if ( ! person->missing() ) std::cout << person->group_label() ;
	      else std::cout << "MISSING";
	    }
	  else 
	    std::cout << ".";

	  std::cout << " ]</p>";
	}

      
      // 2) Report all variants for gene for ths individual
      
      a.indiv_id = ind_value;
       

      // Add any gene transcripts      

      if ( trans.size() > 0 )
	{
	  m.include_loc( loc_set );      
	  for (int i=0; i<trans.size(); i++)
	    m.subset_loc( loc_set , trans[i].name );	      
	}

      
      // Add any included regions 

      for (int r = 0 ; r < a.regions.size() ; r++) 
	m.include_reg( a.regions[r] );

      
            
      VariantGroup vars(m);


      // Accumulate in group, but from a single function call

      g.vardb.iterate( f_display_indiv , &vars, m );


      // Apply existing group function

      g_display_indiv( vars , &a );

    }



  //
  // Misc reports
  //

  // 1) Enumerate all variant sets/supersets

  if ( q == Q_VARSETLIST )
    {
      show_varsets( a ) ;
    }

  
  // 2) Enumerate all genes 

  
  // 3) Tablulate all individuals and phenotypes

  if ( q == Q_PHENOGRID )
    {
      show_indiv(a);      
    }


  //
  // 'Index' (individual-grid) view
  //
  


  //
  // Graphical view
  //
  
  if ( q == Q_GRAPHICAL_VIEW ) 
    {
      show_graphical_view(g,loc_set,q,a,pheno,pwd,project_path,m);
    }

  
  // 
  // All done
  //

  std::cout << "</body></html>";
  

  exit(0);
  
}


void ExomeBrowser::f_display(Variant & var, void *p)
{

  // Aux. data structure
  
  Aux * a = (Aux*)p;

  std::stringstream o1;
  
  o1 << "<tr>";
      
  // Row count
  o1 << "<td>" << ++(a->vcnt) << "</td>";

  // Links per variant
  o1 << "<td>";

  // call pbrowse.cgi, but with all arguments explicitly formed.      

  o1 << a->getURL()->addField("q", "v")
    ->addField("regs", a->reg_list_url )
    ->addField("inc_fltr", a->inc_fltr )
    ->addField("vinc_fltr", a->vinc_fltr )
    ->addField("val", var.coordinate())->printLink("view");
  
  o1 << "</td>";
  
  o1 << "<td>" << Helper::chrCode( var.chromosome() ) << "</td>";

  if ( var.stop() == 0 || var.stop() == var.position() )
    o1 << "<td>" << var.position() << "</td>";
  else
    o1 << "<td>" << var.position() << ".." << var.stop() << "</td>";

  if ( a->single_transcript )
    {
      
      int exon = exon_overlap( a->region , var.position() );
      if ( exon ) 
	o1 << "<td>" << exon << "</td>";
      else
	o1 << "<td><em>intron</em></td>";
    }

  
  // Print name, with link to dbSNP if appropriate
  o1 << "<td>" << rs_link( var.name() ) << "</td>";
  
  o1 << "<td>" << var.pp_reference() << "/" << var.pp_alternate() << "</td>";

  o1 << "<td>" << var.n_samples() << "</td>";
  
  o1 << "<td>" << var.print_meta_filter( "<br>" ) << "</td>";

  SampleVariant & con = var.consensus;


  // Allele count?
  
  if ( a->show_phenotype ) 
    {

      // case/control count for dichotomous traits
      if ( a->g->phmap.type() == PHE_DICHOT )
	{
	  int case_n = 0 , control_n = 0;
	  for (int j=0; j< var.size(); j++)
	    {
	      affType aff = var.ind( j )->affected();
	      
	      if ( var(j).nonreference() )
		{
		  if ( aff == CASE ) ++case_n;
		  else if ( aff == CONTROL ) ++control_n;
		}
	    }
	  o1 << "<td>A=" << case_n << ";U=" << control_n << "</td>";
	}

      // factor count 
      if ( a->g->phmap.type() == PHE_FACTOR )
	{
	  std::map<std::string,int> cnt;
	  for (int j=0; j< var.size(); j++)
	    {
	      Individual * person = var.ind(j);
	      std::string l = ( person->missing() || person->group_label() == "." ) ? 
		"&lt;?&gt;" : person->group_label() ;	      
	      if ( var(j).nonreference() ) cnt[l]++;
	    }
	  o1 << "<td>";
	  std::map<std::string,int>::iterator ii = cnt.begin();
	  while ( ii != cnt.end() )
	    {
	      if ( ii != cnt.begin() ) o1 << "<br>";
	      o1 << ii->first << "=" << ii->second ;
	      ++ii;
	    }
	  o1<< "</td>";	  
	}
      
    }

  // Optional meta-information?
  
  for (int m=0; m < a->mf.size(); m++)
    {
      std::string s = var.print_meta( a->mf[m] , "\t" );
      if ( s == "" ) s = ".";
      std::vector<std::string> sv = Helper::char_split( s , '\t' , true );
      o1 << "<td nowrap>"; 
      for (int i=0;i<sv.size(); i++)
	{ 
	  if (i) o1 << "<br>"; 
	  o1 << ( a->mfpp[ a->mf[m] ] ? sv[i] : pp( sv[i] ) ) ;	  
	}
      o1 << "</td>";
    } 
     
  // Appends? 
  
  if ( a->region_search )
    {
      if ( var.meta.has_field( PLINKSeq::META_LSET() ) )
	{
	  std::string s = PLINKSeq::META_LSET();
	  o1 << "<td nowrap>" << var.meta.print( s ) << "</td>";
	}
      else
	o1 << "<td>" << "n/a" << "</td>";
    }

  // In 'index' / geno-grid mode, add extra colums for individual genotypes

  if ( a->indiv_genogrid ) 
    {
      
      const int n = var.size();

      for (int i=0; i<n; i++)
	{

	  Genotype & g = var(i);
	  
	  // color-coding
	  if ( g.null() ) o1 << "<td align=\"center\" bgcolor=\"gray\">";
	  else if ( g.reference() ) o1 << "<td align=\"center\" bgcolor=\"lightgreen\">";
	  else if ( g.heterozygote() ) o1 << "<td align=\"center\" bgcolor=\"red\">";
	  else if ( g.alternate_homozygote() ) o1 << "<td align=\"center\" bgcolor=\"yellow\">";
	  
	  // actual genotype
	  o1 << "<b>" << ( var.flat() ? var.geno_label( var(i) ) : var.label(i) ) << "</b><br>";
	  
	  // genotype meta-information?

	  if ( var.flat() )
	    {	      
	      std::stringstream ss;
	      ss << g.meta;
	      std::vector<std::string> c = Helper::parse( ss.str() , ";" );
	      for (int k=0;k<c.size(); k++) o1 << c[k] << "<br>";
	    }
	  else
	    {
	      std::cout << "pbrowse does not currently handle projects in which individuals feature more than once"
			<< "</p></body></html>";
	      exit(1);
	      
	    }

	  // done
	  o1 << "</td>";
	}

    }

  o1 << "</tr>";


  // Save table row
  a->table_row[ a->vcnt ] = o1.str() ;

  
}


void ExomeBrowser::f_display_indiv(Variant & vars, void *p)
{
  VariantGroup * d = (VariantGroup*)p; 
  d->force_add(vars);
}

void ExomeBrowser::g_display_indiv(VariantGroup & vars, void *p)
{
  
  // Aux. data structure
  
  Aux * a = (Aux*)p;
  
  const int ni = a->g->indmap.ind_n( a->indiv_id );
  
  // Display variant group as a HTML table
  
  // Overall group details
  
  std::cout << "Individual report for <b>" << a->indiv_id << "</b>";
  if ( vars.name() != "" ) std::cout << "for gene " << vars.name() << "</b>";
  std::cout << "</p>";
  
  std::cout << vars.size() << " variants in total<br>";
  
  if ( vars.size() == 0 ) return;
  
  //vars.prepare();

  // LINKS CHR POS NAME ALT/REF  GENO  GENE_META...
  
  std::cout << "<table border=1>"
       << "<tr><th>#</th><th>Indiv</th><th>Chr</th><th>Pos</th><th>ID</th><th>Ref/Alt</th>";

  // Optional meta-fields?
  
  for (int m=0; m< a->mf.size(); m++)
    std::cout << "<th>" << a->mf[m] << "</th>";
  
  std::cout << "<th>Genotype</th>";
  
  //
  // Get genotype meta-fields
  //

  std::set<std::string> gmeta;
  
  for (int i=0; i < vars.size(); i++)
    {            
      Variant & var = vars(i);
      std::vector<std::string> k = var( ni ).meta.keys();
      for (int z=0; z<k.size(); z++) 
	gmeta.insert(k[z]);      
    }
  

  // We now have a list of all relevant genotypic meta fields, and we
  // know which position in each file this individual is sitting in
  
  std::set<std::string>::iterator i = gmeta.begin();
  while ( i != gmeta.end() )
    {
      std::cout << "<th>" << *i << "</th>";
      ++i;
    }
  
  // End of header row
  
  std::cout << "</tr>";
  


  //
  // Per-variant summary
  //

  for (int i=0; i<vars.size(); i++)
    {
      
      Variant & var = vars(i);

      std::cout << "<tr>";

      // counter
      std::cout << "<td>" << i+1 << "</td>";

      // Links per variant
      std::cout << "<td>";

      // call pbrowse.cgi, but with all arguments explicitly formed.      

      std::cout << a->getURL()->addField("q", "v")
	->addField("val", vars.var(i).coordinate())
	->addField("regs", a->reg_list_url)
	->printLink("view");

      std::cout << "</td>";

      std::cout << "<td>" << Helper::chrCode( vars.var(i).chromosome() ) << "</td>" 
	   << "<td>" << vars.var(i).position() << "</td>";

      // Print name, with link to dbSNP if appropriate
      std::cout << "<td>" << rs_link( vars.var(i).name() ) << "</td>";
           
      std::cout << "<td>" << vars.var(i).pp_reference() << "/" << vars.var(i).pp_alternate() << "</td>" ;


      // Optional meta-information?
      
      for (int m=0; m< a->mf.size(); m++)
	{
	  std::string ms = vars.var(i).print_meta( a->mf[m] , "<br>" );
	  std::cout << "<td>" << ( ms == "" ? "." : ms ) << "</td>";
	}

      // Genotype; genotype meta-fields;

      int gt = var(ni).null() ? 3 : 
	2 - var(ni).allele_count( );
      
      if ( gt == 0 ) 
	std::cout << "<td bgcolor=\"red\">";
      else if ( gt == 1 ) 
	std::cout << "<td bgcolor=\"yellow\">";
      else if ( gt == 2 )
	std::cout << "<td bgcolor=\"lightgreen\">";
      else 
	std::cout << "<td bgcolor=\"gray\">";
      
      std::cout << var.geno_label( var( ni ) )
		<< "</td>";
      
      std::set<std::string>::iterator k = gmeta.begin();
      while ( k != gmeta.end() )
	{
	  if ( (vars.var(i))(ni).meta.hasField( *k ) )
	    {
	      meta_name_t t = *k;
	      std::cout << "<td>" << (vars.var(i))(ni).meta.print( t ) << "</td>";
	    }
	  else
	    std::cout << "<td>.</td>";
	  ++k;
	}
      std::cout << "</tr>";	  

    }
  
  std::cout << "</table>";

} 



//
// Graphical view
//


void ExomeBrowser::show_graphical_view( GStore & g , 
					const std::string & loc_set,
					const QType & q, 
					Aux & a, 
					const std::string & pheno,
					const std::string & pwd , 
					const std::string & project_path , 
					Mask & m )

{



  //
  // Some constants
  //
  
  const int N_RARE_MAC       = 50;
  const int N_TOO_MANY_INDIV = 100;
  const int N_TOO_MANY_TRANS = 30;
  const int N_TOO_MANY_VAR   = 100;


  //
  // Get list of regions to query (either corresponding to a region or transcript
  //
  
  // original input is a.reg_list (regions or genes)

  std::vector<std::vector<Region> > regions;
  std::vector<int2> minmax;
  std::vector<std::string> tok = Helper::whitespace( a.reg_list );

  for (int i=0; i<tok.size(); i++)
    {
      
      std::vector<Region> dummy;
      regions.push_back(dummy);
      const int lst = regions.size() - 1 ;

      bool okay = true;
      Region r(tok[i],okay);
      if ( okay ) regions[lst].push_back( r );
      else
	{
	  std::string tmp = tok[i];
	  // Does this look like a gene name? 
	  std::string genename = tmp;
	  Helper::str2upper( genename );
	  
	  std::set<std::string> trans_names = g.locdb.targetted_lookup_alias( genename , 
									      ExomeBrowser::symbol, 
									      loc_set );
	  
	  // If no matches, assume the original input was in
	  // transcript form, e.g. NA_12345

	  std::vector<std::string> tnames;

	  if ( trans_names.size() == 0 ) 
	    {
	      tnames.push_back( tmp );
	    }
	  else 
	    {
	      
	      // ... otherwise add what we've found from the lookup instead
	      
	      std::set<std::string>::iterator ii = trans_names.begin();
	      while ( ii != trans_names.end() ) 
		{
		  tnames.push_back( *ii ) ; 
		  ++ii;
		}
	    }
	  
	  // get genomic loci for these transcripts
	  
	  std::vector<Region> trans = g.locdb.fetch( loc_set , tnames );
	  for (int i=0;i<trans.size();i++) regions[lst].push_back( trans[i] );
	  
	}
      

      // now find min/max for what will be region, and look up *all* transcripts
      // but do not expand
	 
      int min = 1;
      int max = 300000000;
      
      if ( regions[lst].size() > 0 ) 
	{
	  min = regions[lst][0].start.position();
	  max = regions[lst][0].stop.position();
	  
	  for (int i=0;i<regions[lst].size();i++)
	    {
	      if ( min > regions[lst][i].start.position() ) min = regions[lst][i].start.position() ;
	      if ( max < regions[lst][i].stop.position() ) max = regions[lst][i].stop.position() ;
	    }
	}
      
      minmax.push_back(int2(min,max));
      
      if ( g.locdb.attached() )
	{
	  std::set<Region> r = g.locdb.get_regions( loc_set , Region( regions[lst][0].start.chromosome() , min , max ) );
	  regions[lst].clear();
	  std::set<Region>::iterator ii = r.begin();
	  while ( ii != r.end() ) { regions[lst].push_back( *ii ) ; ++ii; } 
	}
      
    }
  



  //
  // Now we have a list of regions for this query line
  //
  
  const int nreg = regions.size();

  
  // Fixed x-axes
  
  const int xsize = 1200;
  const int border = 150;
  const int rborder = 25;
  const int plotsize = xsize - border - rborder;
  

  // Use one canvas per transcript 
  
  std::vector<Pseq::Helper::Lines> canvas(nreg);

  //
  // 1) loci from core refseq group
  // 2) ref-dbs (e.g. could be SNP intervals)
  // 3) loc-d bs
  // 3) variants and connect individuals with rare variants
  //

  //
  // Iterate over each region / plot 
  //

  for ( int r = 0 ; r < nreg ; r++ ) 
    {
      
      if ( regions[r].size() == 0 ) continue;

      Pseq::Helper::Lines & cnv = canvas[r];
      
      //
      // Determine span of region
      //
      
      int c = regions[r].size();
      
      // find min/max  (note, which is based on original query/transcripts --
      //  new transcripts may have been added, but if these go off the edge of the 
      // plot, we do not automatically expand to include those in the definition
      // of the span )
      
      int min = minmax[r].p1;
      int max = minmax[r].p2;      
      int span = max - min;      

      // a little border, and some constraints

      min -= span*0.02;
      max += span*0.02;
      if ( min < 1 ) min = 1;
      if ( max > 300000000 ) max = 300000000;
      span = max - min + 1 ;

      int chr = regions[r][0].start.chromosome();

      // Set mask to include this region
      // Question: what if other includes/requires have been specified? 
      //           hmm.. should blank those first
      //           excludes should be fine

      Mask m2 = m;
      
      m2.include_reg( Region( chr , min , max ) );      
      
      // ensure ref-db and loc-db appends

      bool add_ref = g.refdb.attached() && a.ref_append.size() ;
      bool add_loc = g.locdb.attached() && a.loc_append.size() ;
      
      // < 100bp  1bp == 10pixels   10:1
      // < 1kb    1bp == 1pixel      1:1
      // 1-10kb   1bp == 0.1 pixel   1:10
      // 10-100kb 1bp == 0.01 pixel  1:100
      // 100-1M          0.0001      1:1000
 
      
      double scale = 1;
      if      ( span > 1e6 ) scale = 0.0001;
      else if ( span > 1e5 ) scale = 0.001;
      else if ( span > 1e4 ) scale = 0.01;
      else if ( span > 1e3 ) scale = 0.1;
      else if ( span < 1e2 ) scale = 10;
      
      // include genotypes? meta-information

      bool load_genotypes =  span < 1e6 ;

      if ( ! load_genotypes ) 
	{
	  m2.load_genotype_data( false );      
	  m2.load_genotype_meta( false );
	  m2.load_variant_meta( false );
	}
  
      
      // show genome sequence? 
      bool show_sequence = scale == 10 && g.seqdb.attached();
            
      // Obtain variants
      VariantGroup vars(m);
      g.vardb.iterate( f_display_indiv , &vars, m2 );      
      const int nv = vars.n_variants();
      const int n  = g.indmap.size();
      
      // Using case/control status?
      bool case_control = g.phmap.type() == PHE_DICHOT;

      // get expected size of plot (y-axis)

      // y-axis will be a certain number of units, where each unit is 20 pixels
      //  1 * genomic bp element labels
      //  1 * axis
      //  1 * spacer
      //  (at small scales) 1 * SEQUENCE TRACK + 1 * spacer
      //  N * transcripts
      //  1 * spacer
      //  1 * variant 
      //  1 * spacer
      //  R * ref-groups 
      //  L * loc-vgroups
      //  1 * spacer
      //  P * people 
      
      // get numbers for N, M and P now


      // TODO
      //  ref-details
      //  varianr meta-info
      //  family data?
      //  annotation
      //  mac/maf representation

      
      //
      // Annoation
      //
      
      bool use_annotation = load_genotypes && g.seqdb.attached() && g.locdb.attached() && ! N_TOO_MANY_TRANS ;
      
      if ( use_annotation ) 
	{

	  std::set<Region> aregs;
	  for (int j=0;j<regions[r].size();j++) aregs.insert( regions[r][j] );
	  if ( ! Annotate::load_transcripts( loc_set , aregs ) ) use_annotation = false;
   
	  std::vector<std::string> annot;  // worst annot per trans.
	  std::map<int,std::map<std::string,int> > change; // variant-N -->  transcript : 0/1/2/3/4/5 int/sil/mis/spl/non/oth
	  std::map<int,std::map<std::string,std::string> > protchange; // variant-N -->  transcript : change
	  
	  for (int v=0;v<nv;v++)
	    {
	      
	      bool exonic = Annotate::annotate( vars(v) );
	      
	      annot.push_back( vars(v).meta.get1_string( PLINKSeq::ANNOT() ) );
	      
	      // 	  plog << var.coordinate() << "\t"
	      // 	       << "transcript" << "\t"
	      // 	       << var.meta.as_string( PLINKSeq::ANNOT_GENE() , "," ) << "\n";
	      
	      // 	  plog << var.coordinate() << "\t"
	      // 	       << "protein" << "\t"
	      // 	       << var.meta.as_string( PLINKSeq::ANNOT_PROTEIN() , "," ) << "\n";
	      
	      
	    }
	}
      
      
      
      
      //
      // Determine size of plot
      //

      // each line is exactly 20 pixels

      const int ystep = 15; 
      
      int ylines = 7;  // ( defaults from above)

      if ( show_sequence ) ylines += 2;  // sequence line?

      // number of transcripts (or collapsed if too many)
      ylines += regions[r].size() < N_TOO_MANY_TRANS
       	? regions[r].size() 
	: 1 ;  

      // ref and loc tracks (and a group spacer, if any)
      ylines += a.ref_append.size() ? a.ref_append.size() + 1 : 0 ;   
      ylines += a.loc_append.size() ? a.loc_append.size() + 1 : 0 ;  
      
      // Get list of individuals with 1+ *rare* variant
      std::vector<int> mac(nv);
      std::vector<double> maf(nv);
      std::map<int,std::set<int> > sites; 
      
      if ( load_genotypes ) 
	{
	  for (int v=0;v<nv;v++)
	    {
	      int c     = 0; // minor allele
	      int c_tot = 0; // total counts	  
	      bool altmin = vars(v).n_minor_allele( &c , &c_tot );      
	      if ( altmin ) { mac[v] = c; maf[v] = (double)c/(double)c_tot; }
	      else { mac[v] = c_tot - c ; maf[v] = 1 - (double)c/(double)c_tot; }
	      if ( mac[v] <= N_RARE_MAC )
		for (int i=0;i<n;i++)
		  if ( vars(v,i).nonreference() ) sites[i].insert(v);
	    }
	}
      
      
      if ( sites.size() > N_TOO_MANY_INDIV )
	{
	  ylines += 2 ; // only show message saying that there are too many people to show
	}
      else if ( sites.size() > 0 )
	ylines += 1 + sites.size();

      
      // add some more for luck (not everthing fully counted)
      ylines += 20;


      // Now we know what size the plot should be 
      const int ysize = ystep * ylines;


      //
      // Create plot
      //
      
      cnv.setobj( "c" + Helper::int2str( r ) , xsize , ysize );
      cnv.box( 1, 1 , xsize , ysize , "white" , "black" , 1 );      
      
      if ( scale < 0.001 ) 
	cnv.text( 8, 12 , 
		  Helper::chrCode(chr) + ":" +  Helper::int2str( min ) + ".." + Helper::int2str( max ) 
		  + "  (" + Helper::dbl2str( (max-min+1)/(double)1e6) + "Mb)", "blue" );
      else if ( scale < 1 ) 
	cnv.text( 8, 12 , 
		  Helper::chrCode(chr) + ":" +  Helper::int2str( min ) + ".." + Helper::int2str( max ) 
		  + "  (" + Helper::dbl2str( (max-min+1)/(double)1e3) + "kb)", "blue" );	
      else 
	cnv.text( 8, 12 , 
		  Helper::chrCode(chr) + ":" +  Helper::int2str( min ) + ".." + Helper::int2str( max ) 
		  + "  (" + Helper::int2str( max-min+1 ) + "bp)", "blue" );

      
      int yoff = ystep * 1.5;


      //
      // Genomic positions
      //
      //

      if ( scale < 0.001 ) 
	for (int t=0;t<10;t++) // do not show last position
	  cnv.text( border + plotsize * ( t/10.0 ) - 20 , yoff , 
		    Helper::dbl2str( 1e-6 * ( min + (t/10.0)*span ) ) + "Mb" , "black" );
      else if ( scale < 1 ) 
	for (int t=0;t<10;t++) // do not show last position
	  cnv.text( border + plotsize * ( t/10.0 ) - 20 , yoff , 
		    Helper::dbl2str( 1e-3 * ( min + (t/10.0)*span ) ) + "kb" , "black" );
      else 
	for (int t=0;t<10;t++) // do not show last position
	  cnv.text( border + plotsize * ( t/10.0 ) - 20 , yoff , 
		    Helper::int2str( ( min + (t/10.0)*span ) ) + "bp" , "black" );
      
      yoff += ystep;


      //
      // Axis and navigation links
      //

      cnv.line( border , yoff , xsize - rborder , yoff , "gray" , 1 ) ;
      for (int t=0;t < 11 ; t++)
        cnv.line(  border + plotsize * ( t/10.0 )  , yoff - 5 , border + plotsize * ( t/10.0 ) , yoff + 5 , "gray" , 1 );      
      
      // Left (shift 20%, 50%)
      
      std::string nreg = Helper::chrCode(chr) + ":" 
	+ Helper::int2str(  min - span*0.2 ) + ".." 
	+ Helper::int2str( max - span*0.2 ) ; 
      
      cnv.text(  border - 19 , yoff - 7 , 
		 "[<]" , 
		 "blue" ,
		 "navlink_left1" ,
		 a.getURL()->addField("q", "gview")
		 ->addField("regs", nreg )
		 ->addField("inc_fltr", a.inc_fltr )
		 ->addField("vinc_fltr", a.vinc_fltr )->printLink() );
      
      nreg = Helper::chrCode(chr) + ":" 
	+ Helper::int2str(  min - span*0.5 ) + ".." 
	+ Helper::int2str( max - span*0.5 ) ; 
      
      cnv.text(  border - 22 , yoff + 9 , 
		 "[<<]" , 
		 "blue" ,
		 "navlink_left2" ,
		 a.getURL()->addField("q", "gview")
		 ->addField("regs", nreg )
		 ->addField("inc_fltr", a.inc_fltr )
		 ->addField("vinc_fltr", a.vinc_fltr )->printLink() );


      // Right (shift 20%, 50%)

      nreg = Helper::chrCode(chr) + ":" 
	+ Helper::int2str(  min + span*0.2 ) + ".." 
	+ Helper::int2str( max + span*0.2 ) ; 
      
      cnv.text(  xsize-19 , yoff - 7 , 
		 "[>]" ,
		 "blue" ,
		 "navlink_right1" ,
		 a.getURL()->addField("q", "gview")
		 ->addField("regs", nreg )
		 ->addField("inc_fltr", a.inc_fltr )
		 ->addField("vinc_fltr", a.vinc_fltr )->printLink() );
	      
      nreg = Helper::chrCode(chr) + ":" 
	+ Helper::int2str(  min + span*0.5 ) + ".." 
	+ Helper::int2str( max + span*0.5 ) ; 
      
      cnv.text(  xsize - 23 , yoff + 9 , 
		 "[>>]" , 
		 "blue" ,
		 "navlink_right2" ,
		 a.getURL()->addField("q", "gview")
		 ->addField("regs", nreg )
		 ->addField("inc_fltr", a.inc_fltr )
		 ->addField("vinc_fltr", a.vinc_fltr )->printLink() );


      // zoom links
      
      for (int t=0; t < 10 ; t++)
	{

	  // zoom in
	  if ( scale < 10 )
	    {
	      
	      std::string nreg = Helper::chrCode(chr) + ":" 
		+ Helper::int2str( ( min + (t/10.0)*span ) ) + ".." 
		+ Helper::int2str( ( min + ((t+1)/10.0)*span ) ) ; 
	      
	      cnv.text(  border + plotsize * ( t/10.0 ) + ( plotsize/22 ) , yoff + 12 , 
			 "[+]" , 
			 "gray" ,
			 "navlink_plus_" + Helper::int2str( t ) , // some unique name
			 a.getURL()->addField("q", "gview")
			 ->addField("regs", nreg )
			 ->addField("inc_fltr", a.inc_fltr )
			 ->addField("vinc_fltr", a.vinc_fltr )->printLink() );
	    }
	  

	  // zoom out (allow 20Mb max)

	  if ( span < 20e6 )
	    { 
	      // check chr bounds
	      int chr_min = 1 ; 
	      int chr_max = 300000000 ; 
	      if  ( g.seqdb.attached() ) 
		{
		  std::map<int,int2> cmm = g.seqdb.getMinMax();
		  if ( cmm.find( chr ) != cmm.end() ) 
		    {
		      chr_min = cmm[chr].p1;
		      chr_max = cmm[chr].p2;
		    }
		}
	      
	      int center = min + ((t+0.5)/10.0)*span ;
	      if ( chr_min < center - 10 * span ) chr_min = center - 10 * span ;
	      if ( chr_max > center + 10 * span ) chr_max = center + 10 * span ;
	      
	      std::string nreg = Helper::chrCode(chr) + ":" 
		+ Helper::int2str( chr_min ) + ".." 
		+ Helper::int2str( chr_max ) ; 
	      
	      cnv.text(  border + plotsize * ( t/10.0 ) + ( plotsize/22 )+1 , yoff - 5 , 
			 "[-]" , 
			 "gray" ,
			 "navlink_minus_" + Helper::int2str( t ) , // some unique name
			 a.getURL()->addField("q", "gview")
			 ->addField("regs", nreg )
			 ->addField("inc_fltr", a.inc_fltr )
			 ->addField("vinc_fltr", a.vinc_fltr )->printLink() );
	    }

	}

      yoff += ystep;

      
      //
      // Track 3 : spacer
      //
	  
      yoff += ystep;

      //
      // Track 4 : Sequence
      //

      if ( show_sequence ) 
	{
	  std::string s = g.seqdb.lookup( chr , min , max ); 

	  if ( s.size() != span ) 
	    cnv.text( 10 , yoff , "problem reading sequence data...");
	  else
	    {
	      // will be at least 10 pixels per base
	      for (int i=0;i<s.size();i++)
		{		  
		  if      ( s[i] == 'A' || s[i] == 'a' ) cnv.text( border + ( i/(double)span ) * plotsize -5 , yoff , "A" , "red" );
		  else if ( s[i] == 'C' || s[i] == 'c' ) cnv.text( border + ( i/(double)span ) * plotsize -5, yoff , "C" , "blue" );
		  else if ( s[i] == 'G' || s[i] == 'g' ) cnv.text( border + ( i/(double)span ) * plotsize -5, yoff , "G" , "green" );
		  else if ( s[i] == 'T' || s[i] == 't' ) cnv.text( border + ( i/(double)span ) * plotsize -5, yoff , "T" , "yellow" ); 
		  else cnv.text( border + ( i / (double)span ) * plotsize -5, yoff , s.substr(i,1) , "gray" );		  
		}
	    }
	}

      // and add another spacer
      yoff += ystep;


      //
      // Transcripts
      //
      
      cnv.text( 5 , yoff , "Transcripts" , "blue" ) ;   

      yoff += ystep;

      bool show_full_transcript = regions[r].size() < N_TOO_MANY_TRANS ; 

      for (int tr = 0 ; tr < regions[r].size(); tr++ ) 
	{
	  
	  Region & reg = regions[r][tr];
	  
	  int p1 = border + ( ( reg.start.position() - min ) / (double)span ) * plotsize ;
	  int p2 = border + ( ( reg.stop.position()  - min ) / (double)span ) * plotsize ;

	  // either link to make full transcript
	  if ( show_full_transcript )
	    cnv.box( p1 , yoff-5 , p2 , yoff+5 , "white","white",1,
		     "trans_link_" + Helper::int2str(tr) , 		
		     a.getURL()->addField("q", "gview")
		     ->addField("regs", reg.name )
		     ->addField("inc_fltr", a.inc_fltr )
		     ->addField("vinc_fltr", a.vinc_fltr )->printLink() );
   
	  cnv.line( p1 , yoff , p2 , yoff , "gray" , 1 ) ;

	  
	  const int s = reg.subregion.size();

	  // exons first, then CDS
	  for (int ss=0;ss<s;ss++)
	    {
	      
	      int p1 = border + ( ( reg.subregion[ss].start.position() - min) / (double)span ) * plotsize;
	      int p2 = border + ( ( reg.subregion[ss].stop.position() - min) / (double)span ) * plotsize;
	      
	      int w = 1;
	      std::string bgcol = "gray";
	      std::string fgcol = "gray";
	      
	      if ( ! reg.subregion[ss].exon() ) continue;
	      
	      if ( show_full_transcript ) 
		cnv.box( p1, yoff-w , p2 , yoff+w , bgcol , fgcol , 1, 
			 "el_" + Helper::int2str(tr) + "_" + Helper::int2str(ss) , 
			 a.getURL()->addField("q", "gview")
			 ->addField("regs", reg.subregion[ss].coordinate() )
			 ->addField("inc_fltr", a.inc_fltr )
			 ->addField("vinc_fltr", a.vinc_fltr )->printLink() );
	      else
		cnv.box( p1, yoff - w , p2 , yoff + w , bgcol , fgcol );
	    }
	  
	  // CDS 
	  for (int ss=0;ss<s;ss++)
	    {
	      
	      int p1 = border + ( ( reg.subregion[ss].start.position() - min) / (double)span ) * plotsize;
	      int p2 = border + ( ( reg.subregion[ss].stop.position() - min) / (double)span ) * plotsize;
	      
	      int w = 4;
	      std::string bgcol = "lightblue";
	      std::string fgcol = "black";

	      if      ( reg.subregion[ss].start_codon() ) { bgcol = fgcol = "red";w=7; } 
	      else if ( ! reg.subregion[ss].CDS() ) continue;
		   
	      if ( show_full_transcript ) 
		cnv.box( p1, yoff-w , p2 , yoff+w , bgcol , fgcol , 1, 
			 "el_" + Helper::int2str(tr) + "_" + Helper::int2str(ss) , 
			 a.getURL()->addField("q", "gview")
			 ->addField("regs", reg.subregion[ss].coordinate() )
			 ->addField("inc_fltr", a.inc_fltr )
			 ->addField("vinc_fltr", a.vinc_fltr )->printLink() );
	      else
		cnv.box( p1, yoff - w , p2 , yoff + w , bgcol , fgcol );
	    }


	  // transcript name in left 

	  if ( show_full_transcript ) 
	    {
	      cnv.box(2 , yoff-5 , border , yoff+5 , "white" , "white" );
	      cnv.box(xsize - rborder + 1 , yoff-5 , xsize - 1 , yoff+5 , "white" , "white" );
	      cnv.text( 5,yoff+3 , 
			reg.altname != reg.name && reg.altname != "" && reg.altname != "." 
			? reg.altname + "(" + reg.name + ")" 
			: reg.name  );	  
	      yoff += ystep;
	    }
	  

	}
      
      yoff += ystep;

      //
      // Reference DBs
      //
      
      if ( add_ref ) 
	{

	  // get reference variants in this region	  
	  cnv.text( 5 , yoff , "Reference variants" , "blue" ) ;   
	  yoff += ystep;

	  for (int i=0;i<a.ref_append.size();i++)
	    {
	      std::set<RefVariant> rvar = g.refdb.lookup( Region(chr,min,max) , a.ref_append[i] , 10000 );
	      // did we hit limit?
	      if ( rvar.size() >= 10000 ) 
		{
		  cnv.text(5,yoff,a.ref_append[i] + " (more items in this interval than can be displayed)" , "gray" );
		  rvar.clear();
		}
	      else
		cnv.text(5,yoff,a.ref_append[i] + " (" + Helper::int2str(rvar.size()) + " items)" , "gray" );
	      
	      bool show_meta =  rvar.size() < 20 ;
	      std::set<RefVariant>::iterator ii = rvar.begin();
	      while ( ii != rvar.end() )
		{
		  int p1 = border + ( ( ii->start() - min ) / (double)span ) * plotsize; 
		  int p2 = border + ( ( ii->stop() - min ) / (double)span ) * plotsize; 
		  
 		  cnv.line( p1 , yoff-5 , 
			    p1 , yoff+5 , 
			    "black" ); 
		  

		  if ( show_meta ) 
		    {
		      std::stringstream ss ;
		      ss << *ii << "; " << ii->value() ;
		      cnv.text( p1 + 3 , yoff+5 , ss.str() , "gray");
		      yoff += ystep;
		    }
		  
		  //cnv.box(2 , yoff-1 , border , yoff+11 , "white" , "white" );
		  
		  ++ii;
		}
	      yoff += ystep;
	    }

	  yoff += ystep;
	}

      
      //
      // Locus DBs
      //

      if ( add_loc ) 
	{
	  
	  cnv.text( 5 , yoff , "Reference intervals/loci" , "blue" ) ;   
	  
	  yoff += ystep;

	  for (int i=0;i<a.loc_append.size();i++)
	    {

	      std::set<Region> aloc = g.locdb.get_regions( a.loc_append[i] , Region(chr,min,max) );
	      
	      cnv.box( 2 , yoff , border , yoff+20 , "white" , "white" );
	      cnv.box(xsize - rborder + 1 , yoff , xsize - 1 , yoff+20 , "white" , "white" );
	      cnv.text( 10 ,yoff,a.loc_append[i] + " (" + Helper::int2str( aloc.size() ) + " intervals)", "gray" );
	      
	      std::set<Region>::iterator ii = aloc.begin();
	      while ( ii != aloc.end() ) 
		{
		  int p1 = border + ( ( ii->start.position() - min ) / (double)span ) * plotsize ;
		  int p2 = border + ( ( ii->stop.position()  - min ) / (double)span ) * plotsize ;
		  cnv.line( p1 , yoff + 5 , p1 , yoff + 15 , "black" , 1 ); 
		  cnv.line( p1 , yoff + 10 , p2 , yoff + 10 , "black" , 1 ); 
		  cnv.line( p2 , yoff + 5 , p2 , yoff + 15 , "black" , 1 ); 
		  std::stringstream ss;
		  ss << ii->name;
		  cnv.text( p2 + 2 , yoff + 15 , ss.str() , "gray" );

		  yoff += ystep;
		  ++ii;
		}

	      yoff += ystep;
	    }
	}
      
      yoff += ystep;
      
      
      
      //
      // Variant ticks
      //
      
      cnv.text( 5 , yoff , "Variants" , "blue" ) ;   
      yoff += ystep;
      cnv.text( 5 , yoff , " (" + Helper::int2str( nv ) + " in region)" , "gray" ) ;   
      for (int v=0;v<nv;v++)
	{
	  int p1 = border + ( ( vars(v).position() - min ) / (double)span ) * plotsize; 
	  
	  if ( nv < N_TOO_MANY_VAR ) 
	    {
	      if ( show_sequence )
		{
		  std::string alt = vars(v).alternate();
		  std::string col = "gray";
		  if ( alt.size() == 1 ) // not multiallelic, indel
		    {
		      if      ( alt == "A" || alt == "a" ) col = "red";
		      else if ( alt == "C" || alt == "c" ) col = "blue";
		      else if ( alt == "G" || alt == "g" ) col = "green";
		      else if ( alt == "T" || alt == "t" ) col = "yellow";
		    }
		 
		  
		  cnv.line( p1 , yoff , p1 , ystep*5 , "gray" );
		  
		  cnv.box( p1-8 , yoff-8 , 
			   p1+8 , yoff+8 , 
			   "white" , 
			   "black",
			   1,
			   "var_" + Helper::int2str(v) , 
			   a.getURL()->addField("q", "v")
			   ->addField("regs", a.reg_list_url )
			   ->addField("inc_fltr", a.inc_fltr )
			   ->addField("vinc_fltr", a.vinc_fltr )
			   ->addField("val", vars(v).coordinate() )->printLink() );
		  
		  cnv.text( p1 - 4 , yoff + 4  , 
			    vars(v).alternate() , 
			    col );
		  
		}		  
	      else
		cnv.box( p1-5 , yoff-5 , 
			 p1+5 , yoff+5 , 
			 "white" , 
			 "black",
			 1,
			 "var_" + Helper::int2str(v) , 
			 a.getURL()->addField("q", "v")
			 ->addField("regs", a.reg_list_url )
			 ->addField("inc_fltr", a.inc_fltr )
			 ->addField("vinc_fltr", a.vinc_fltr )
			 ->addField("val", vars(v).coordinate() )->printLink() );
	    }
	  else
	    cnv.line( p1 , yoff -5, 
		      p1 , yoff +5 , 
		      "black" , 
		      1 );
	  
	}

      yoff += ystep;

    
      //
      // Annotation
      //

      
  
      //
      // Phenotype/frequency (MAF/MAC indicators)
      //

      
      if ( nv < 30 ) 
	{

	  yoff += ystep;

	  cnv.text( 5 , yoff , "Minor allele count/frequency" , "blue" );


	  if ( case_control ) 
	    {
	    
  	      for (int v=0;v<nv;v++)
		{
		  yoff += ystep;
		  
		  int p1 = border + ( ( vars(v).position() - min ) / (double)span ) * plotsize;

		  cnv.line( p1, yoff - 5 , p1 , yoff + 5 , "black" ) ;

		  // need to recalculate C/C frequencies

		  // cases
		  int c1 = 0 , c1_tot = 0;
		  double maf1 = 0;
		  bool altmin1 = vars(v).n_minor_allele( &c1 , &c1_tot , &maf1 , CASE );

		  int c2 = 0 , c2_tot = 0;
		  double maf2 = 0;
		  bool altmin2 = vars(v).n_minor_allele( &c2 , &c2_tot , &maf2 , CONTROL );
		  
		  if ( !altmin1 ) { c1 = c1_tot - c1; maf1 = 1 - maf1; }
		  if ( !altmin2 ) { c2 = c2_tot - c2; maf2 = 1 - maf2; }

		  if ( maf1 > 0.01 || maf2 > 0.01 ) 
		    cnv.text( p1 + 4, yoff +5 , 
			      "MAF=" + Helper::dbl2str( maf1 , 3 ) + "/" + Helper::dbl2str( maf2 , 3) , "gray" );
		  else
		    cnv.text( p1 + 4 , yoff +5 , 
			      "AAC=" + Helper::int2str( c1 ) + "/" + Helper::int2str( c2 ) , "gray" );
		  
		}	      
	    
	    }
	  else
	    {
	      for (int v=0;v<nv;v++)
		{
		  yoff += ystep;

		  int p1 = border + ( ( vars(v).position() - min ) / (double)span ) * plotsize;

		  cnv.line( p1, yoff - 5 , p1 , yoff + 5 , "black" ) ;
		  
		  if ( maf[v]>0.01 ) 
		    cnv.text( p1 , yoff +5 , 
			      "MAF=" + Helper::dbl2str( maf[v] , 3 ) , "gray" );
		  else
		    cnv.text( p1 , yoff +5 , 
			      "AAC=" + Helper::int2str( mac[v] ) , "gray" );
		 		 
		  //   cnv.text( p1+1 , yoff +5 , Helper::int2str( mac[v] ) + " alt-alleles", "green" );
		}	      
	      
	    }
	}
      
      yoff += ystep;

      

      //
      // Variant ticks
      //
      
      yoff += ystep;
      cnv.text( 5 , yoff , "Individuals (rare variant genotypes)" , "blue" ) ;   
      yoff += ystep;

      
      //
      // Calculate summary counts
      //

      if ( sites.size() < N_TOO_MANY_INDIV ) 
	{

	  std::string nreg = Helper::chrCode(chr) + ":" +  Helper::int2str( min ) + ".." + Helper::int2str( max ) ; 

	  if ( case_control ) // do caes first )
	    {

	      std::map<int,std::set<int> >::iterator ii = sites.begin();	  
	      while ( ii != sites.end() )
		{	      
		  if ( g.indmap(ii->first)->affected() != CASE ) { ++ii; continue; }

		  // link to table-browser for individual
		  cnv.box( 5 , yoff , xsize - rborder , yoff + 20 , "white" , "white" , 1 , 
			   "indiv_" + Helper::int2str( ii->first ) , 
			   a.getURL()->addField("q", "i")
			   ->addField("regs", nreg )
			   ->addField("ind", g.indmap.ind( ii->first )->id() ) 
			   ->printLink() );

		  cnv.text( 5 , yoff+10 , g.indmap.ind( ii->first )->id() + " (A)", "red" );
		  
		  std::set<int>::iterator jj = ii->second.begin();	      
		  int tmin = *jj;
		  int tmax = tmin;
		  while ( jj != ii->second.end() ) 
		    {		  
		      int p1 = border + ( ( vars( *jj ).position() - min ) / (double)span  ) * plotsize ;		  
		      cnv.box( p1 , yoff-5 , p1 , yoff+5 , "black" );
		      tmax = *jj;
		      ++jj;
		    }
		  
		  int p1 = border + ( ( vars(tmin).position() - min  ) / (double)span ) * plotsize;
		  int p2 = border + ( ( vars(tmax).position() - min  ) / (double)span ) * plotsize;		  
		  cnv.line( p1 , yoff , p2 , yoff , "gray" );		  
		  yoff += ystep;	      
		  ++ii;	      
		}

	      yoff += ystep;

	      // controls
	      ii = sites.begin();	  
	      while ( ii != sites.end() )
		{	      
		  if ( g.indmap(ii->first)->affected() != CONTROL ) { ++ii; continue; }

		  cnv.box( 5 , yoff , xsize - rborder , yoff + 20 , "white" , "white" , 1 , 
			   "indiv_" + Helper::int2str( ii->first ) , 
			   a.getURL()->addField("q", "i")
			   ->addField("regs", nreg )
			   ->addField("ind", g.indmap.ind( ii->first )->id() ) 
			   ->printLink() );

		  cnv.text( 5 , yoff , g.indmap.ind( ii->first )->id() + " (U)", "blue" );
		  
		  std::set<int>::iterator jj = ii->second.begin();	      
		  int tmin = *jj;
		  int tmax = tmin;
		  while ( jj != ii->second.end() ) 
		    {		  
		      int p1 = border + ( ( vars( *jj ).position() - min ) / (double)span  ) * plotsize ;		  
		      cnv.box( p1 , yoff-5 , p1 , yoff+5 , "black" );		  
		      tmax = *jj;
		      ++jj;
		    }
		  
		  int p1 = border + ( ( vars(tmin).position() - min  ) / (double)span ) * plotsize;
		  int p2 = border + ( ( vars(tmax).position() - min  ) / (double)span ) * plotsize;		  
		  cnv.line( p1 , yoff , p2 , yoff , "gray" );		  
		  yoff += ystep;	      
		  ++ii;	      
		}
	      
	      yoff += ystep;

	      // missing data
	      ii = sites.begin();	  
	      while ( ii != sites.end() )
		{	      
		  if ( g.indmap(ii->first)->affected() != UNKNOWN_PHE ) { ++ii; continue; }

		  cnv.box( 5 , yoff , xsize - rborder , yoff + 20 , "white" , "white" , 1 , 
			   "indiv_" + Helper::int2str( ii->first ) , 
			   a.getURL()->addField("q", "i")
			   ->addField("regs", nreg )
			   ->addField("ind", g.indmap.ind( ii->first )->id() ) 
			   ->printLink() );
		  
		  cnv.text( 5 , yoff , g.indmap.ind( ii->first )->id() + " (?)", "gray" );
		  
		  std::set<int>::iterator jj = ii->second.begin();	      
		  int tmin = *jj;
		  int tmax = tmin;
		  while ( jj != ii->second.end() ) 
		    {		  
		      int p1 = border + ( ( vars( *jj ).position() - min ) / (double)span  ) * plotsize ;		  
		      cnv.box( p1 , yoff-5 , p1 , yoff+5 , "black" );		  
		      tmax = *jj;
		      ++jj;
		    }
		  
		  int p1 = border + ( ( vars(tmin).position() - min  ) / (double)span ) * plotsize;
		  int p2 = border + ( ( vars(tmax).position() - min  ) / (double)span ) * plotsize;		  
		  cnv.line( p1 , yoff , p2 , yoff , "gray" );		  
		  yoff += ystep;	      
		  ++ii;	      
		}
	      
	      yoff += ystep;
	      
	    }
	  else
	    {
	      std::map<int,std::set<int> >::iterator ii = sites.begin();	  
	      while ( ii != sites.end() )
		{	      
		  cnv.text( 5 , yoff , g.indmap.ind( ii->first )->id() );
		  
		  std::set<int>::iterator jj = ii->second.begin();	      
		  int tmin = *jj;
		  int tmax = tmin;
		  while ( jj != ii->second.end() ) 
		    {		  
		      int p1 = border + ( ( vars( *jj ).position() - min ) / (double)span  ) * plotsize ;		  
		      cnv.box( p1 , yoff , p1 , yoff+10 , "black" );		  
		      tmax = *jj;
		      ++jj;
		    }
		  
		  int p1 = border + ( ( vars(tmin).position() - min  ) / (double)span ) * plotsize;
		  int p2 = border + ( ( vars(tmax).position() - min  ) / (double)span ) * plotsize;
		  
		  // using a case/control phenotype?
		  int ph = 0;
		  if ( case_control && g.indmap( ii->first)->affected() == CASE ) ph = 2;
		  else if ( case_control && g.indmap( ii->first)->affected() == CONTROL ) ph = 1;		  
		  cnv.line( p1 , yoff+5 , p2 , yoff+5 , ph == 0 ? "gray" : ph == 2 ? "red" : "blue"  );
		  
		  yoff += ystep;	      
		  ++ii;	      
		}
	    }

	}

      
      // re-write box around edge
      cnv.line( 1, 1 , xsize , 1 );
      cnv.line( xsize , 1 , xsize , ysize );
      cnv.line( xsize , ysize , 1 , ysize );
      cnv.line( 1  , ysize , 1 , 1 );

    }

  

  //
  // Write to page 
  //

  std::string preamble = "";
  for (int r=0;r<nreg;r++)
    preamble += canvas[r].preamble();
  preamble += Pseq::Helper::Lines::loader() ;
  
  write_html_header( preamble );

  // special header to enable the graphical functions
  std::cout << "<body onload=\"canvasdraw()\">";

  // write top of page
  write_start_page(g,loc_set,q,a,pheno,pwd,project_path);
  
  // now place each transcript
  for (int r=0;r<nreg;r++)
    std::cout << "<center>" << canvas[r].place() << "</center><br><hr>";
    
  return;
} 




//
// Helper functions
//


int ExomeBrowser::exon_overlap( const Region & reg , int pos )
{
  int chr = reg.start.chromosome();
  for (int s=0; s < reg.subregion.size(); s++)
    if ( reg.subregion[s].overlaps( Region(chr,pos,pos) ) ) return s+1;  
  return 0; // code for intronic variant
}



std::string ExomeBrowser::rs_link(const std::string & label )
{  
  if ( label == "0" ) return "n/a";  // fix old missing = 0 code      
  if ( label.size() > 2 && label.substr(0,2) == "rs" ) 
    return "<a href=\"http://www.ncbi.nlm.nih.gov/sites/entrez?db=snp&cmd=search&term=" 
      + label + "\" target=\"_blank\">" + label + "</a>"; 
  return label;
}

std::string ExomeBrowser::pp(const std::string & str , const int len ) 
{
  if ( str.size() < len ) return str;
  return "<abbr title=\"" + str + "\">" + str.substr(0,len) + "...</abbr>";
}

void ExomeBrowser::make_gene_list(Aux * a)
{

  // true param indicates to fetch alternate names
  std::set<std::string> genes = a->g->locdb.fetch_names( a->loc_set );
  
  std::set<std::string>::iterator ii = genes.begin();
  while ( ii != genes.end() )
    {
      std::cout << a->getURL()->addField("q", "r")->addField("regs", *ii )->printLink( *ii ) << "<br/>";
      ++ii;
    }
  std::cout << "</body></html>";
  exit(0);
}


void ExomeBrowser::make_phe_list(Aux * a)
{
  std::map<std::string,std::vector<std::string> > m = a->g->inddb.fetch_phenotype_info();
  std::map<std::string,std::vector<std::string> >::iterator i = m.begin();
  std::cout << "<h3>Available phenotypes</h3>";
  std::cout << "<table border=1><tr>"
       << "<th align=left>Phenotype</th>"
       << "<th align=left>Type</th>"
       << "<th align=left>Description</th></tr>";
  while ( i != m.end() )
    {
      std::cout << "<tr>"
	   << "<td>" << i->first << "</td>"
	   << "<td>" << i->second[0] << "</td>"
	   << "<td>" << i->second[1] << "</td>"
	   << "</tr>";
      ++i; 
    }
  std::cout << "</table></body></html>";
  exit(0);
}



void ExomeBrowser::make_mf_list(Aux * a)
{
  
  std::map<int,std::string> f = a->g->vardb.fetch_files();
  std::map<int,std::string>::iterator i = f.begin();
  while ( i != f.end() )
    {
      std::cout << "<h3><font color=blue>" << i->second << "(" << i->first << ")</font></h3>";
      // Meta-types
      std::cout << "<h4>Variant meta-information fields</h4>";
      std::cout << "<table border=1><th align=left>Meta-field</th>"
 	   << "<th align=left>Type</th>"
 	   << "<th align=left>Length</th>"
 	   << "<th align=left>Group</th>"
 	   << "<th align=left>Description</th>"
 	   << "</tr>";
      
      std::vector<std::map<std::string,std::string> > m = a->g->vardb.fetch_metatypes( i->first );
      for (int j=0; j<m.size(); j++)
	{
	  std::cout << "<tr><td>";
	  
	  if ( MetaMeta::static_variant( m[j]["NAME"] ) )
	    std::cout << "<font color=\"green\"> " << m[j]["NAME"] << "</font></td><td>";
	  else
	    std::cout << m[j]["NAME"] << "</td><td>";
	  
	  std::cout << m[j]["TYPE"] << "</td><td>"
	       << m[j]["NUM"] << "</td><td>"
	       << m[j]["GRP"] << "</td><td>"
	       << m[j]["DESC"] << "</td></tr>";
	}
      std::cout << "</table>";

      // headers
      std::cout << "<h4>Header rows in VCF</h4>";
      std::cout << "<table border=1>"
 	   << "<th align=left>Key</th>"
 	   << "<th align=left>Value</th></tr>";
      
     std::vector<std::map<std::string,std::string> > h = a->g->vardb.fetch_headers( i->first );
      for (int j=0; j<h.size(); j++)
	{
	  // do not print out first line of VCF #CHROM...
	  std::string ky =  h[j]["KEY"];
	  if ( ky.size() > 4 && ky.substr(0,4) == "#CHR" )
	    continue;
	  std::cout << "<tr><td>" 
	       << ky << "</td><td>"
	       << h[j]["VALUE"] << "</td></tr>";
	}
      std::cout << "</table>";      
      std::cout <<"<hr>";
      ++i;
    }
  std::cout << "</table></body></html>";
  exit(0);
}


void ExomeBrowser::make_locset_list(Aux * a)
{
  
  std::set<GroupInfo> g = a->g->locdb.group_information();
  std::set<GroupInfo>::iterator i = g.begin();

  std::cout << "<h3>Available locus sets</h3>";
  std::cout << "<table border=1><tr>"
       << "<th align=left>Name</th>"
       << "<th align=left>Number of entries</th>"
       << "<th align=left>Decscription</th>"
       << "</tr>";
  
  while ( i != g.end() )
    {
      std::cout <<"<tr>"
	   << "<td>" << i->name << "</td>"
	   << "<td>" << a->g->locdb.count( i->idx )<< "</td>"
	   << "<td>" << i->description << "</td>"
	   << "</tr>";
      ++i; 
    }
  std::cout << "</table></body></html>";
  exit(0);
}

void ExomeBrowser::make_proj_summary(Aux * a)
{
  std::cout << "<h3>Project Summary</h3>"
	    << "<p><pre>"
	    << a->g->summary( false ) 
	    << "</pre></p>";
  std::cout << "</table></body></html>";
  exit(0);
}


void ExomeBrowser::write_html_header( const std::string & head_preamble )
{

  std::cout << "Content-type: text/html\n\n"
	    << "<html><head>";
 

  std::cout << "<title>PLINK/SEQ browser</title>"
	    << "<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />";
  
  // CSS

  std::cout << "<style media=\"screen\" type=\"text/css\">"
	    << "table { border-collapse: collapse; border: 1px solid #666666; font: normal 11px verdana, arial, helvetica, sans-serif;"
	    << "color: #363636; background: #f6f6f6; text-align:left;}"
	    << "caption {"
	    << "text-align: center;"
	    << "font: bold 16px arial, helvetica, sans-serif;"
	    << "background: transparent;"
	    << "padding:6px 4px 8px 0px;"
	    << "color: #CC00FF;"
	    << "   text-transform: uppercase;"
	    << "}"
	    << "thead, tfoot {"
	    << "  text-align:left;"
	    << "height:30px;"
	    << "}"
	    << "thead th, tfoot th {"
	    << "padding:5px;"
	    << "}"
	    << "table a {"
	    << "color: #333388;"
	    << "  text-decoration:underline;"
	    << "}"
	    << "table a:hover {"
	    << "  text-decoration:underline;"
	    << "}"
	    << "tr.odd {"
	    << "background: #f171f1;"
	    << "}"
	    << "tbody th, tbody td {"
	    << "padding:5px;"
	    << "}"
	    << "</style>";

  // some basic JavaScript functions

  std::cout << "<SCRIPT LANGUAGE=\"JavaScript\">"
	    << " function checkAll(field) {"
	    << " for (i = 0; i < field.length; i++)"
	    << " field[i].checked = true ; }"
	    << " function uncheckAll(field) { for (i = 0; i < field.length; i++) field[i].checked = false ; }"
	    << "</script>";


  // any on-the-fly-defined functions related to the graphical 
  // views

  std::cout << head_preamble ;  

  // end of HEAD
  
  std::cout << "</head>";

}


void ExomeBrowser::write_start_page( const GStore & g , 
				     const std::string & loc_set,
				     const QType & q, 
				     Aux & a, 
				     const std::string & pheno ,
				     const std::string & pwd , 
				     const std::string & project_path )
{

  std::cout << "<form name=\"myform\" action=\"pbrowse.cgi\" method=\"GET\"> ";
  
  std::cout << "<table width=100% CELLPADDING=0 CELLSPACING=0>"
	    << "<tr><td width=50% valign=center align=left>"
	    << "<h1><a style=\"color:black;text-decoration:none;\" href=\""; 
  
  // std::cout << a.getURL()->addField("q", "r")
  //   ->addField("passwd",pwd)
  //   ->printURL();

  std::cout << "\">p<font color=\"darkred\">browse</font></h1>"
	    << "</td><td width=50% valign=center align=right>";
  
 
  // project/pwd specification
  
  if ( g.pwd( pwd ) )
    std::cout << "(" << a.getURL()->addField("q", "psummary")->addField("passwd",pwd)->printLink("show project summary") << ")" << "<br>";

  std::cout << "Project: <input type=\"text\" size=\"50\" name=\"proj\" value=\"" 
	    << Helper::html_encode( project_path ) << "\">"
	    << "<br>"
	    << "Password: <input type=\"text\" size=\"50\" name=\"passwd\" value=\"" 
	    << Helper::html_encode( pwd ) << "\">"
	    

    // end of header
	    << "</td></tr></table>"

	    << "<hr>";
  

  if ( q == Q_GENELIST 
       || q == Q_PHELIST 
       || q == Q_METALIST 
       || q == Q_LOCSETLIST 
       || q == Q_PROJSUMMARY )
    {
      
      if ( ! Helper::fileExists( project_path ) )
	{
	  std::cout << "File [ " 
		    << project_path 
		    << " ] could not be found "
		    << "</BODY></HTML>";
	  exit(0);
	}

      
      if ( ! g.pwd( pwd ) ) 
	{
	  Helper::halt( "<b>access denied: password does not match</b>" );
	  exit(0);
	}

      if ( q == Q_GENELIST ) 
	make_gene_list(&a);
      else if ( q == Q_PHELIST )
	make_phe_list(&a);
      else if ( q == Q_METALIST )
	make_mf_list(&a);
      else if ( q == Q_LOCSETLIST )
	make_locset_list(&a);
      else if ( q == Q_PROJSUMMARY )
	make_proj_summary(&a);

      exit(0);

    }


  
  //
  // Draw main query box, with saved defaults
  //
  
  
  std::cout << "<input type=\"hidden\" NAME=\"proj\" "
	    << " value=\"" << project_path << "\"> ";
  
  std::cout << "<table width=100%><tr>";
  
  // Hidden value to indicate query type (needed?)

  std::cout << "<input type=\"hidden\" name=\"q\" value=\"";
  if ( q == Q_REGION || q == Q_ERROR ) std::cout << "r";
  else if ( q == Q_VARIANT ) std::cout << "v";
  else if ( q == Q_INDIV ) std::cout << "i";
  std::cout << "\">";

  
  //
  // Panels
  //


  std::cout << "<td width=22% valign=top>";


  //
  // Query
  //

  std::cout << "<p><b>Variants</b>: gene/transcript ID (" << loc_set << " " 
	    << a.getURL()->addField("q", "lslist")->printLink("change") << "|"
	    << a.getURL()->addField("q", "glist")->addField("pheno", pheno)->printLink("list") << "),<br>region or variant ID";

  // add separate page to change this (on lslist page)
//   std::cout << "<br><input type=\"text\" size=\"45\" ";  
//   if ( loc_set != "" ) 
//     std::cout << " value=\"" << Helper::html_encode(loc_set) << "\"";
//   std::cout << " name=\"loc\"></p>";

  
  //
  // Region list
  // 

  //  std::cout << "<p>Additional regions and genes<br>";
  std::cout << "<textarea ";
  std::cout << "name=\"regs\" rows=\"7\" cols=\"30\">"
	    << a.reg_list 
    //	    << Helper::html_encode( a.reg_list )
	    << "</textarea></p>";

  
  //
  // Second table column
  //

  std::cout << "</td><td width=4% valign=top> &nbsp ";
  std::cout << "</td><td width=22% valign=top>";



  //
  // Optional meta-fields
  //
  
  std::cout << "<p>Optional variant meta-fields ";
  std::cout << "(" << a.getURL()->addField("q", "mflist")->addField("pheno", pheno)->printLink("list") << ")";

  std::cout << "<br><input type=\"text\" size=\"45\" name=\"meta\"";
  
  if ( a.mf.size() != 0 ) 
    {
      std::cout << " value=\"";
      for (int m=0; m< a.mf.size(); m++) 
	{
	  if ( m>0 ) std::cout << " ";
	  std::cout << Helper::html_encode(a.mf[m]);
	}
      std::cout << "\"";
    }
  std::cout << ">";
  std::cout << "</p>";


  //
  // Phenotype
  //

  std::cout << "<p>Optional case/control phenotype";

  std::cout << "(" << a.getURL()->addField("q", "plist")->addField("pheno", pheno)->printLink("list") << ")";

  
  std::cout << "<br><input type=\"text\" size=\"45\" name=\"pheno\"";
  
  if ( pheno != "" ) 
    {
      std::cout << " value=\"";
      std::cout << Helper::html_encode(pheno) << "\"";
    }
  
  std::cout << "</p>";


  //
  // Individual list/mask
  //

  std::cout << "<p>Individuals include "
	    << "(" << a.getURL()->addField("q", "pgrid")->printLink("list") << ")";	
  std::cout << "<br><input type=\"text\" size=\"45\" name=\"indiv_list\"";
  std::cout << " value=\""<< Helper::html_encode( a.indiv_list_url ) << "\"";  
  std::cout << ">";
  std::cout << "</p>";
  


  //
  // Third column
  //

  std::cout << "</td><td width=4% valign=top> &nbsp ";
  std::cout << "</td><td width=22% valign=top>";
  
  //
  // Variant set filters
  //
  
  std::cout << "<p>Variant sets/superset filters "
	    << "(" << a.getURL()->addField("q", "varsetlist")->printLink("list") << ")";	
  std::cout << "<br><input type=\"text\" size=\"45\" name=\"varset\"";
  std::cout << " value=\""<< Helper::html_encode( a.varset_url ) << "\"";  
  std::cout << ">";
  std::cout << "</p>";

  //
  // Reference appends
  //
  
  std::cout << "<p>Reference DB filters/appends";
  std::cout << "<br><input type=\"text\" size=\"45\" name=\"ref_append\"";
  std::cout << " value=\""<< Helper::html_encode( a.ref_append_url ) << "\"";  
  std::cout << ">";
  std::cout << "</p>";
  

  //
  // Locus DB appends
  //

  std::cout << "<p>Locus DB filters/appends";
  std::cout << "<br><input type=\"text\" size=\"45\" name=\"loc_append\"";
  std::cout << " value=\""<< Helper::html_encode( a.loc_append_url ) << "\"";  
  std::cout << ">";
  std::cout << "</p>";




  //
  // Fourth table column
  //

  std::cout << "</td><td width=4% valign=top> &nbsp ";
  std::cout << "</td><td width=22% valign=top>";


  //
  // Mask specification
  //
  
  std::cout << "<p>Optional mask specification";
  std::cout << "<br><input type=\"text\" size=\"45\" name=\"masks\"";
  
  if ( a.msk.size() != 0 ) 
    {
      std::cout << " value=\"";
      for (int m=0; m< a.msk.size(); m++) 
	{
	  if ( m>0 ) std::cout << " ";
	  std::cout << Helper::html_encode(a.msk[m]);
	}
      std::cout << "\"";
    }
  std::cout << ">";
  std::cout << "</p>";
  

  //
  // Mask specification
  //
  
  std::cout << "<p>Include filter (<tt>include</tt> in mask)";
  std::cout << "<br><input type=\"text\" size=\"45\" name=\"inc_fltr\"";
  std::cout << " value=\""<< Helper::html_encode( a.inc_fltr ) << "\"";  
  std::cout << ">";
  std::cout << "</p>";
  

  //
  // V-include mask
  //

  std::cout << "<p>Variant include filter (<tt>v-include</tt> in mask)";
  std::cout << "<br><input type=\"text\" size=\"45\" name=\"vinc_fltr\"";
  std::cout << " value=\""<< Helper::html_encode( a.vinc_fltr ) << "\"";  
  std::cout << ">";
  std::cout << "</p>";



  //
  // End of table
  //

  std::cout << "</td></tr></table>";

  //
  // Submit buttons
  //
  
  std::cout << " <br> <input type=\"submit\" name=\"getgene\" value=\"Variant table\"> "
	    << " <input type=\"submit\" name=\"gview\" value=\"Regional Figure\"> "
	    << " <input type=\"submit\" name=\"indgrid\" value=\"Genotype Grid\"> "
	    << " <input type=\"submit\" name=\"pgrid\" value=\"Individual table\"> ";
    
  std::cout << " </form> ";
  std::cout << "<hr> ";


  //
  // Draw query 
  //
  

  if ( ! Helper::fileExists( project_path ) )
    {
      std::cout << "File [ " << project_path << " ] could not be found "
		<< "</BODY></HTML>";
      exit(0);
    }
  
    if ( q == Q_ERROR )
    {
      //      std::cout << "Problem processing input...";
      
      if ( ! g.pwd(pwd) ) 
	std::cout << "Please enter a password to access this project<br>";

      std::cout << "</body></html>";
      exit(0);
    }  

}


void ExomeBrowser::index_grid( Aux & a )
{
  // show multiple individuals (presumably a smallish number) for selected genotypes in a table



}


void ExomeBrowser::show_indiv( Aux & a )
{


  // Get *all* individuals from the VARDB
  
  std::map<std::string,std::set<std::string> > ind2file;

  int cnt = 0;
  if ( ! a.g->vardb.attached() ) return;

  std::map<int,std::string> files = a.g->vardb.fetch_files();
  std::map<int,std::string>::iterator i = files.begin();
  while ( i != files.end() )
    {
      std::vector<std::string> inds = a.g->vardb.fetch_individuals( i->first );
      for (int j = 0 ; j < inds.size(); j++)
	ind2file[ inds[j] ].insert( a.g->vardb.file_tag( i->first ) + " : " + i->second );
      ++i;
    }
  
  const int n = ind2file.size();
  
  if ( a.indiv_list.size() == 0 ) 
    std::cout << "All " << n << " individuals selected<br>";
  else 
    std::cout << a.indiv_list.size() << " of " << n << " individuals selected ("
	      << a.getURL()->addField("q", "pgrid")
      ->addField("indiv_list", "" )
      ->printLink( "clear list" )
	      << ")<br>";
      
  // make table
  
  std::cout << "<table border=1><tr>"
	    << "<th>Selected</th>"
	    << "<th>Primary ID</th>"
	    << "<th>FID</th>"
	    << "<th>IID</th>"
	    << "<th>Paternal ID</th>"
	    << "<th>Maternal ID</th>"
	    << "<th>Sex</th>"
	    << "<th>Phenotype(s)</th>"
	    << "</tr>";
  

  std::map<std::string,std::set<std::string> >::iterator ii = ind2file.begin();

  while ( ii != ind2file.end() )
    {

      std::string id = ii->first;

      // get sex of a particular individual
      sType s = a.g->inddb.sex( id );

      Individual person = a.g->inddb.fetch( id );
      
      // checked?

      if ( a.indiv_list.size() == 0 || a.indiv_list.find(id) != a.indiv_list.end() )
	std::cout << "<td><input type=\"checkbox\" disabled=\"disabled\" name=\"indsel\" value=\""<<id<<"\" checked></td>";
      else
	std::cout << "<td><input type=\"checkbox\" disabled=\"disabled\" name=\"indsel\" value=\""<<id<<"\"></td>";

      // ID
      
      std::cout << "<td>" 
		<< a.getURL()->addField("q", "pgrid")
	->addField("indiv_list", a.indiv_list_url == "" ? id : a.indiv_list_url + " " + id )
	->printLink( id ) 	
		<< "</td>";
      
       
      // IID/FID
      std::cout << "<td>" << person.fid() << "</td>"
		<< "<td>" << person.iid() << "</td>";

      // Paternal/maternal ID
      std::cout << "<td>" << person.father() << "</td>"
		<< "<td>" << person.mother() << "</td>";
      
      // Sex

      if ( s == MALE ) std::cout << "<td>Male</td>";
      else if ( s == FEMALE ) std::cout << "<td>Female</td>";
      else std::cout << "<td>.</td>";

      // Phenotypes
      std::cout << "<td>";
      std::stringstream ss;
      ss << person.meta;
      std::vector<std::string> c = Helper::parse( ss.str() , ";" );
      for (int k=0;k<c.size(); k++) std::cout << c[k] << "<br>";
      std::cout << "</td>";      
      
      // Done
      std::cout << "</tr>";
      ++ii;
    }
              
  std::cout << "</table>";
    
  
}


void ExomeBrowser::show_varsets( Aux & a )
{

  std::cout << "<h3>Variant sets</h3>";

  std::vector<std::string> sets = a.g->vardb.get_sets();
  
  std::cout << "<table border=1><tr><th>Selected</th><th>Set</th><th>Number of variants</th><th>Description</th></tr>";

  for (int s = 0 ; s < sets.size(); s++)
    {
      std::cout << ( a.varset_set.find( sets[s] ) != a.varset_set.end() ? 
		     "<td><input type=\"checkbox\" disabled=\"disabled\" name=\"setsel\" value=\""+sets[s]+"\" checked></td>" :
		     "<td><input type=\"checkbox\" disabled=\"disabled\" name=\"setsel\" value=\""+sets[s]+"\"></td>" );
      
      std::cout << "<td>" 
		<< a.getURL()->addField("q", "varsetlist")
	->addField("varset", a.varset_url == "" ? sets[s] : a.varset_url + " " + sets[s] )
	->printLink( sets[s] ) 	
		<< "</td>";
      
      std::cout << "<td>" << a.g->vardb.get_set_size( sets[s] ) << "</td><td>" 
		<< a.g->vardb.get_set_description( sets[s] ) << "</td></tr>";
    }

  std::cout << "</table>";

  //
  // Variant super-sets
  //

  std::cout << "<h3>Variant super-sets</h3>";

  std::vector<std::string> ss = a.g->vardb.get_supersets();
  if ( ss.size() == 0 ) std::cout << "<br><em>none</em><br>";
  else
    {
      std::cout << "<table border=1><tr><th>Selected</th><th>Set</th><th>Number of sets</th><th>Sets</th><th>Description</th></tr>";
            
      for (int s = 0 ; s < ss.size(); s++)
	{
	  std::vector<std::string> sets = a.g->vardb.get_sets( ss[s] );
	  
	  std::cout << "<tr><td>" ;

	  std::cout << "<td>" 
		    << a.getURL()->addField("q", "varsetlist")
	    ->addField("varset", a.varset_url == "" ? ss[s] : a.varset_url + " " + ss[s] )
	    ->printLink( ss[s] ) 	
		    << "</td>";
     	
	  std::cout << "<td>" << sets.size() << "</td>";

	  for (int t = 0 ; t < sets.size(); t++)
	    std::cout << sets[t] << "<br>";
	  std::cout << "</td>" << a.g->vardb.get_superset_description( ss[s] ) << "</td>";
	  std::cout << "</tr>";
	}
      std::cout << "</table>";
    }
}
