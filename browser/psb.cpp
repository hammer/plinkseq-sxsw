#include "psb.h"
#include "cgi.h"

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
  std::string loc_set = "refseq";
  std::string genename = "";
  QType q = Q_ERROR;
  std::string var_value = "";
  std::string ind_value = "";
  std::string chr_code = "";
  int chr = 0, bp1 = 0, bp2 = 0;
  std::string pheno = "";
  bool from_top = false;
  

  // Auxilliary information to send to variant printing functions

  Aux a;
  
  
  //
  // Start HTML form and any output
  //

  std::cout << "Content-type: text/html\n\n"
	    << "<html><head><title>PLINK/SEQ Browser</title>"
	    << "<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />"
	    << "<title>PLINK/SEQ genetics library</title>";
  
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
    //	    << "background: #f1f1f1;"
	    << "background: #f171f1;"
	    << "}"
	    << "tbody th, tbody td {"
	    << "padding:5px;"
	    << "}"
	    << "</style>";
  
  // end of HEAD
	    std::cout << "</head><body>";


  
  if ( cgi ) 
    {
      for (i=0; cgivars[i]; i+= 2)
	{
	  
	  std::string str = cgivars[i];
	  
	  if ( str == "q" )
	    {
	      std::string s = cgivars[i+1];
	      
	      if ( s == "v" ) q = Q_VARIANT;
	      else if ( s == "g" ) q = Q_GENE;
	      else if ( s == "i" ) q = Q_INDIV;
	      else if ( s == "r" ) q = Q_REGION;
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
	  
	  if ( str == "val" ) 
	    var_value = cgivars[i+1];
	  
// 	  if ( str == "gene" ) 
// 	    {
// 	      genename = cgivars[i+1];
// 	      Helper::str2upper( genename );
// 	    }
	  
	  if ( str == "ind" ) 
	    ind_value = cgivars[i+1];
	  
	  if ( str == "proj" ) 
	    {
	      project_path = cgivars[i+1];
	      a.add_form_value("proj", project_path );
	    }
	  
	  if ( str == "loc" ) loc_set = cgivars[i+1];
	  
	  if ( str == "regs" )
	    {
	      a.reg_list = a.reg_list_url = cgivars[i+1];
	      while ( 1 ) 
		{
		  int p = a.reg_list.find( "," );
		  if ( p == std::string::npos ) break;
		  a.reg_list.replace( p , 1 , "\n" );
		}
	      std::replace( a.reg_list_url.begin(), a.reg_list_url.end() , '\n' , ',' );
	    }

	  if ( str == "masks" )
	    a.msk = Helper::parse( cgivars[i+1] , " " );
	  
	  if ( str == "meta" ) 
	    a.mf = Helper::parse( cgivars[i+1] , " ," );
	  
	  if ( str == "pheno" )
	    pheno = cgivars[i+1];
	  
	}
      
      
      if ( from_top ) 
	{
	  q = Q_GENE;
	}

    }
  
  
  // ** Free anything that needs to be freed **/
  
  if ( cgi )
    {
      for (i=0; cgivars[i]; i++) free(cgivars[i]) ;
      free(cgivars) ;
    }

  
  
  //
  // Start page 
  //
  
  std::cout << "<table width=100% CELLPADDING=0 CELLSPACING=0>"
	    << "<tr><td width=50% valign=top align=left>"
	    << "<h1>PLINK<font color=\"darkred\">SEQ</font> exome browser</h1>"
	    << "</td><td width=50% valign=top align=right>"
	    << "<font size=-1>"
	    << "Project: <em>" << project_path << "</em></font>"
	    << "<br>(<a href=\"pbrowse.cgi?q=psummary&" << a.print_form_value("proj")
	    << "\">show summary</a>)"
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

      GStore g;

      a.g = &g;
      a.loc_set = loc_set;
      g.set_project( project_path ) ;

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
  
  std::cout << "<form name=\"myform\" action=\"pbrowse.cgi\" method=\"GET\"> ";
  
  std::cout << "<input type=\"hidden\" NAME=\"proj\" "
	    << " value=\"" << project_path << "\"> ";
  
  std::cout << "<table width=100%><tr>";
  
  // Hidden value to indicate query type (needed?)

  std::cout << "<input type=\"hidden\" name=\"q\" value=\"";
  if ( q == Q_GENE || q == Q_ERROR ) std::cout << "g\">";
  else if ( q == Q_VARIANT ) std::cout << "v\">";
  else if ( q == Q_INDIV ) std::cout << "i\">";
  

  
  //
  // Panels
  //


  std::cout << "<td width=30% valign=top>";


  //
  // Query 
  //

  std::cout << "<p><b>Gene ID</b> (symbol or NM_012345) ";
  std::cout << "(<a href=\"pbrowse.cgi?q=glist&" << a.print_form_value("proj")
	    << "&meta=" << a.mf_print()
	    << "&masks=" << a.msk_print() 
	    << "&pheno=" << pheno
	    << "&regs=" << a.reg_list_url
	    << "\">list</a>)";
  
//   std::cout << "<br><input type=\"text\" size=\"45\" name=\"gene\"";
//   if ( genename != "" )
//     std::cout << " value=\"" << genename << "\"";
//   std::cout << ">";
//   std::cout << "</p> ";
  

  //
  // Region list
  // 

  //  std::cout << "<p>Additional regions and genes<br>";
  std::cout << "<textarea ";
  std::cout << "name=\"regs\" rows=\"7\" cols=\"30\">"
	    << a.reg_list 
	    << "</textarea></p>";

  //
  // Submit buttons
  //
  
  std::cout << " <br> <input type=\"submit\" name=\"getgene\" value=\"Fetch\"> "
	    << " <input type=\"reset\" value=\"Reset\"> ";

  
  //
  // Second table column
  //

  std::cout << "</td><td width=5% valign=top> &nbsp ";
  std::cout << "</td><td width=30% valign=top>";



  //
  // Optional meta-fields
  //
  
  std::cout << "<p>Optional variant meta-fields ";
  std::cout << " (<a href=\"pbrowse.cgi?q=mflist&" << a.print_form_value("proj")
	    << "&meta=" << a.mf_print()
	    << "&masks=" << a.msk_print()
	    << "&pheno=" << pheno
	    << "&regs=" << a.reg_list_url
	    << "\">list</a>) ";			  

  std::cout << "<br><input type=\"text\" size=\"45\" name=\"meta\"";
  
  if ( a.mf.size() != 0 ) 
    {
      std::cout << " value=\"";
      for (int m=0; m< a.mf.size(); m++) 
	{
	  if ( m>0 ) std::cout << " ";
	  std::cout << a.mf[m];
	}
      std::cout << "\"";
    }
  std::cout << ">";
  std::cout << "</p>";


  //
  // Phenotype
  //

  std::cout << "<p>Optional case/control phenotype";

  std::cout << " (<a href=\"pbrowse.cgi?q=plist&" << a.print_form_value("proj")
	    << "&meta=" << a.mf_print()
	    << "&masks=" << a.msk_print() 
	    << "&pheno=" << pheno
	    << "&regs=" << a.reg_list_url
	    << "\">list</a>) ";			  
  
  std::cout << "<br><input type=\"text\" size=\"45\" name=\"pheno\"";
  
  if ( pheno != "" ) 
    {
      std::cout << " value=\"";
      std::cout << pheno << "\"";
    }
  
  std::cout << "</p>";

  // Gene-set
  std::cout << "<p>Gene set ";
  std::cout << "(<a href=\"pbrowse.cgi?q=lslist&" << a.print_form_value("proj")
       << "\">list</a>)";			  

  std::cout << "<br><input type=\"text\" size=\"45\" ";  
  if ( loc_set != "" ) 
    std::cout << " value=\"" << loc_set << "\"";
  std::cout << " name=\"loc\"></p>";
    



  //
  // Third column
  //

  std::cout << "</td><td width=5% valign=top> &nbsp ";
  std::cout << "</td><td width=30% valign=top>";
  

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
	  std::cout << a.msk[m];
	}
      std::cout << "\"";
    }
  std::cout << ">";
  std::cout << "</p>";
  

  //
  // include mask
  //


  //
  // V-include mask
  //



  
  //
  // End of table
  //

  std::cout << "</td></tr></table>";


  std::cout << " </form> ";
  std::cout << "<hr> ";


  //
  // Draw query 
  //

  std::vector<std::string> toks = Helper::parse( a.reg_list , " \n" );
  if ( toks.size() == 1 ) 
    {
      genename = toks[0];
      Helper::str2upper( genename );
    }
  else genename = "";

  if ( a.reg_list == "" )
    {
      std::cout << "No genes or regions specified...</BODY></HTML>";
      exit(0);
    } 
  
  if ( ! Helper::fileExists( project_path ) )
    {
      std::cout << "File [ " << project_path << " ] could not be found "
		<< "</BODY></HTML>";
      exit(0);
    }
  
    if ( q == Q_ERROR )
    {
      std::cout << "Problem processing input...";
      std::cout << "</body></html>";
      exit(0);
    }
  


  ////////////////////////////////////////////////////////////////////////////
  

  //
  // Set up basic stuff to perform query
  //


  GStore g;

  a.g = &g;

  a.loc_set = loc_set;

  g.set_project( project_path );

  // Initial Mask objects

  std::string mstr;
  for (int i = 0 ; i < a.msk.size(); i++)
    mstr += ( i ? " " : "" ) + a.msk[i];


  // Add some additional mask items to speed things up, and 
  // check that the output isn't too large; note: the 'no-geno'
  // assumes that complex mask queries that use genotype data 
  // aren't applied, e.g. 
  // include=" g( DP >= 10 ) > 0.1 "
  
  if ( q == Q_GENE || q == Q_REGION ) 
    {
      // if no mask, safe to add 'no-geno', otherwise, we should keep in case
      if ( mstr == "" ) 
	mstr = "limit=5000 no-geno " + mstr;
      else
	mstr = "limit=5000 " + mstr;      
    }
  
  //  std::cout << "m = ["<<mstr<<"]\n";

  Mask m( mstr );
  
  g.indmap.populate( g.vardb, g.phmap, m );


  //
  // Expand wildcard?
  //

  std::set<std::string> onecopy;
  for (int m=0; m<a.mf.size(); m++)
    {
      if ( a.mf[m] == "*" )
	{
	  a.mf.clear();
	  std::map<int,std::string> f = g.vardb.fetch_files();
	  std::map<int,std::string>::iterator i = f.begin();
	  while ( i != f.end() )
	    {
	      std::vector<std::map<std::string,std::string> > m = g.vardb.fetch_metatypes( i->first );
	      for (int j=0; j<m.size(); j++)
		{
		  if ( m[j]["GRP"] == "Variant" )
		    {
		      std::string mval = m[j]["NAME"];
		      if ( onecopy.find( mval ) == onecopy.end() )
			{
			  onecopy.insert(mval);
			  a.mf.push_back(mval);
			}
		    }
		}
	      ++i;
	    }
	  break;
	}
    }



  //
  // If a gene-name has been specified, then find the actual transcript(s) that match
  //

  // Translate symbol into 1+ transcripts, and pick the first

  std::set<std::string> trans_names = g.locdb.targetted_lookup_alias( genename , "symbol" , "refseq" ) ;   

  std::vector<std::string> tnames;
  std::set<std::string>::iterator ii = trans_names.begin();
  while ( ii != trans_names.end() ) 
    {
      tnames.push_back( *ii ) ; 
      ++ii;
    }
  
  

  //
  // Get transcript information from database
  //
  
  std::vector<Region> trans = g.locdb.fetch( loc_set , tnames );


  //
  // If this name matches an alternate name, assume it is a gene
  // symbol; find alias
  //

  a.multi_transcripts = trans.size() > 0;
  
  std::string main_gene = genename;

  if ( ! a.multi_transcripts ) 
    {
      Region reg = g.locdb.get_region( loc_set , main_gene ) ;
      if ( reg.name == main_gene ) 
	{
	  main_gene = reg.altname; 
	}     
    }
  
  // If only a single transcript, change to this now
  
  if ( trans.size() == 1 ) 
    {
      a.multi_transcripts = false;
      genename = trans[0].name;
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

  a.other_genes.clear();
  a.regions.clear();
  a.extended_search = false;

  if ( a.reg_list != "" )
    {

      std::vector<std::string> tok = Helper::whitespace( a.reg_list );
      
      for (int i=0; i<tok.size(); i++)
	{

	  bool okay = true;

	  Region r(tok[i],okay);
	  if ( okay ) 
	    {
	      a.regions.push_back( r );
	      a.extended_search = true;
	    }
	  else
	    {
	      // look for as a gene (assuming upper case for all IDs)
	      std::string tmp = tok[i];
	      Helper::str2upper( tmp );
	      a.other_genes.push_back( tmp );
	      a.extended_search = true;
	    }
	}
     
    }

  
  
  if ( q == Q_GENE 
       && a.regions.size() > 0 
       && genename == "" ) 
    {      
      q = Q_REGION;
    }
  



  //
  // For a gene-based query, some versbose output regarding transcript info, etc
  //

  if ( q == Q_GENE ) 
    {

	if ( a.multi_transcripts )
	    std::cout << "Found " << trans.size() 
		      << " transcript(s) matching gene name <b>" 
		      << main_gene << "</b></p>";
	
	std::cout << "<pre><font size=-1>";
	
	for (int r=0; r<trans.size(); r++)
	{
	    
	    std::cout << trans[r].altname << "   " ;
	    std::cout << "<a href=\"pbrowse.cgi?q=g&" << a.print_form_value("proj")
		      << "&gene="<< trans[r].name	       
		      << "&meta=" << a.mf_print()
		      << "&masks=" << a.msk_print() 
		      << "&pheno=" << pheno
		      << "&regs=" << a.reg_list_url
		      << "\">" << trans[r].name <<"</a>  "
		      << Helper::chrCode(trans[r].start.chromosome()) << ":" 
		      << trans[r].start.position() << ".."
		      << trans[r].stop.position() ;
	    std::cout <<"<br>";
	    if ( !cgi ) std::cout << "\n";
	}
	
      std::cout << "</font></pre>";
            

      //
      // Now, 'genename' will be original search term
      //
      
      if ( ! a.multi_transcripts )
	{
	  
	    Region reg = g.locdb.get_region( loc_set , genename ) ;
	    
	    if ( reg.start.chromosome() == 0 && ! a.extended_search )
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
			    std::cout << "<a href=\"pbrowse.cgi?q=g&" << a.print_form_value("proj")
				      << "&gene="<< t->name
				      << "&meta=" << a.mf_print()
				      << "&masks=" << a.msk_print() 			
				      << "&pheno=" << pheno
				      << "&regs=" << a.reg_list_url
				      << "\">"<<t->name<<"("<<t->altname<<")</a> ";			  
			  }
		      }
		  ++t;
		}
	      std::cout << "<br>";
	    }
	  
	  std::cout << "</font></pre>";
      
	}

      std::cout << "<hr>";      
    }


  if ( a.extended_search ) 
    {
      std::cout << "<p>Additionally ";
      if ( genename != "" ) std::cout << "requiring ";
      else std::cout << "including ";
      std::cout << "the following regions, and including the following genes</p>";
      
      if ( a.regions.size() > 0 ) 
	{
	  std::cout << "<pre><b>Regions:</b><br>";
	  for (int i=0; i<a.regions.size(); i++)
	    std::cout << "   " << Helper::chrCode( a.regions[i].chromosome()) << ":"
		      << a.regions[i].start.position() << ".."
		      << a.regions[i].stop.position() << "<br>";
	  std::cout << "</pre>";
	}
      
      if ( a.other_genes.size() > 0 ) 
	{
	  std::vector<std::string> cc = a.other_genes;
	  a.other_genes.clear();
	  
	  std::cout << "<pre><b>Additional genes:</b><br>";
	  for (int i=0; i< cc.size(); i++)
	    {
	      
	      std::set<std::string> trans = 
		g.locdb.targetted_lookup_alias( cc[i] , "symbol" , "refseq" ) ;   
	      
	      // if no aliases match, assume this is a refseq transcript name
	      if ( trans.size() == 0 ) 
		{
		  a.other_genes.push_back( cc[i] );
		  std::cout << "   " << Helper::sw( cc[i] , 10 ) << "<br>";
		}
	      else // add all matching aliases
		{
		  std::cout << "   " << Helper::sw( cc[i] , 10 ) << " : ";
		  
		  std::set<std::string>::iterator ii = trans.begin();
		  while ( ii != trans.end() ) 
		    {
		      a.other_genes.push_back( *ii );
		      std::cout << *ii << " ";
		      ++ii;
		    }		 
		  std::cout << "<br>";
		}
	    }
	  std::cout << "</pre>";
	}

      std::cout << "<hr>";
    }


  //
  // In gene, or regional mode
  //
  

  if ( q == Q_GENE || q == Q_REGION ) 
    {
      
      
      a.table_row.clear();
      
      a.vcnt = 0;
        
      if ( q == Q_GENE ) 
	{
	  
	  m.include_loc( loc_set );
      
	  if ( a.multi_transcripts )
	    {
	      for (int i=0; i<trans.size(); i++)
		m.subset_loc( loc_set , trans[i].name );
	      a.genename = main_gene;
	    }
	  else
	    {
	      m.subset_loc( loc_set , genename );      
	      a.genename = genename;
	    }
	}
      

      // If any regions specified, add as a requirement
      
      if ( a.regions.size() > 0 )
	{
	  
	  if ( q == Q_GENE ) 
	    {
	      for (int r = 0 ; r < a.regions.size() ; r++) 
		m.require_reg( a.regions[r] );
	    }
	  else
	    {
	      for (int r = 0 ; r < a.regions.size() ; r++) 
		m.include_reg( a.regions[r] );	    
	    }
	}


      //
      // Append gene names in extended search mode
      //

      if ( a.extended_search )
	m.append_loc( loc_set );

      // Other genes specified?

      for (int r = 0 ; r < a.other_genes.size() ; r++) 
	{
	  m.subset_loc( loc_set , a.other_genes[r] );	    
	}
      

      g.vardb.iterate( f_display , &a, m );

      // Now table should be sorted in order, with duplicates removed, in aux structure;

      // string headers;
      // map<int,std::string> table_row;
      
      std::cout << a.headers << "</p>";
      // LINKS CHR POS NAME ALT/REF FileID QUAL INFO(filter) SAMPLE_CNT VMETA(compressed)
      
      std::cout << "<table border=1>"
	   << "<tr><th>#</th><th>Indiv</th><th>Chr</th><th>Pos</th>";
      
      if (  ! ( a.extended_search || a.multi_transcripts ) )
	std::cout << "<th>Exon</th>";
      
      std::cout << "<th>ID</th><th>Ref/Alt</th>"
	   << "<th># samples</th>"      
	   << "<th>Filter</th>";      
      
      
      // Optional case/control counts?
      
      if ( a.show_phenotype && g.phmap.type() == PHE_DICHOT )
	std::cout << "<th>C/C count</th>";
      
  
      // Optional meta-fields?
      
      for (int m=0; m< a.mf.size(); m++)
	std::cout << "<th>" << a.mf[m] << "</th>";

      if ( a.extended_search ) 
	std::cout << "<th>Locus</th>";
      
     
      // End of header row
      
      std::cout << "</tr>";
      
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
		<< "<a href=\"pbrowse.cgi?q=g&" << a.print_form_value("proj")
		<< "&gene="<< genename  
		<< "&meta=" << a.mf_print()
		<< "&masks=" << a.msk_print() 
		<< "&pheno=" << a.phenotype_name 
		<< "&regs=" << a.reg_list_url
		<< "\">gene report</a>";
      if ( genename != "" ) 
	std::cout << " for <b>" << genename << "</b>";
      
      std::cout << "<hr>";
      


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
       std::cout << "<tr><td>Chromosome</td><td>" << var.chromosome() << "</td></tr>";
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
	      
	     
	      if ( gt == 3 && ! var(i).null() ) continue;
	      else if ( var(i).allele_count( ) != 2-gt ) continue;
	      
	
	      // Include link to individual-report

	      std::cout << "<tr><td>" 
			<< "<a href=\"pbrowse.cgi?q=i&" << a.print_form_value("proj")
			<< "&gene="<< genename  
			<< "&meta=" << a.mf_print()
			<< "&masks=" << a.msk_print() 
			<< "&ind=" << id_list[i] 
			<< "&pheno=" << a.phenotype_name 
			<< "&regs=" << a.reg_list_url
			<< "\">"
			<< id_list[i] 
			<< "</a></td>";
	      
	      
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
			  std::cout << "<td>" << var(i).meta.print( t ) << "</td>";
			}
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
  //
  //

  if ( q == Q_INDIV )
    {

      a.genename = genename;

      std::cout << "Back to "
		<< "<a href=\"pbrowse.cgi?q=g&" << a.print_form_value("proj")
		<< "&gene="<< genename  
		<< "&meta=" << a.mf_print()
		<< "&masks=" << a.msk_print() 
		<< "&pheno=" << a.phenotype_name 
		<< "&regs=" << a.reg_list_url
		<< "\">gene report</a>";
      if ( genename != "" ) 
	std::cout << " for <b>" << genename << "</b>";
      
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
       
      m.include_loc( loc_set );

      if ( a.multi_transcripts )
	{
	  for (int i=0; i<trans.size(); i++)
	    {
	      m.subset_loc( loc_set , trans[i].name );	      
	    }
	  a.genename = main_gene;
	}
      else
	{
	  m.subset_loc( loc_set , genename );      
	  a.genename = genename;	  
	}
      

      // If any regions specified, add as a requirement
      for (int r = 0 ; r < a.regions.size() ; r++) 
	{
	  m.include_reg( a.regions[r] );
	}
      
            
      VariantGroup vars(m);

      // Accumulate in group, but from a single function call
      g.vardb.iterate( f_display_indiv , &vars, m );

      // Apply existing group function
      g_display_indiv( vars , &a );

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

  o1 << "<a href=\"pbrowse.cgi?q=v&" << a->print_form_value("proj")
     << "&gene=" << a->genename
     << "&meta=" << a->mf_print()
     << "&masks=" << a->msk_print() 
     << "&pheno=" << a->phenotype_name  
     << "&val=" << var.coordinate()
     << "&regs=" << a->reg_list_url 
     << "\">view</a>";

  o1 << "</td>";
  
  o1 << "<td>" << var.chromosome() << "</td>" 
     << "<td>" << var.position() << "</td>";
  

  if ( ! ( a->extended_search || a->multi_transcripts ) )
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
  
  if ( a->show_phenotype && a->g->phmap.type() == PHE_DICHOT )
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
      o1 << "<td>" << case_n << "/" << control_n << "</td>";
    }
  
  // Optional meta-information?
  
  for (int m=0; m < a->mf.size(); m++)
    {
      std::string s = var.print_meta( a->mf[m] , "\t" );
      std::vector<std::string> sv = Helper::char_split( s , '\t' , true );      
      o1 << "<td nowrap>"; 
      for (int i=0;i<sv.size(); i++)
	{ 
	  if (i) o1 << "<br>"; 
	  o1 << pp( sv[i] );	  
	}
      o1 << "</td>";
    } 
     
  // Appends? 
  
  if ( a->extended_search ) 
    {
      if ( var.meta.has_field( PLINKSeq::META_LSET() ) )
	{
	  std::string s = PLINKSeq::META_LSET();
	  o1 << "<td nowrap>" << var.meta.print( s ) << "</td>";
	}
      else
	o1 << "<td>" << "n/a" << "</td>";
    }


//   // Get normal appends from meta-fields
//   for ( int i = 0; i < a->app.size() ; i++)
//     {
//       RefVariant rv = a->g->refdb.lookup( var , a->app[i] );
//       o1 << "<td nowrap>" << rv << "</td>";
//     }
  
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
  
  std::cout << "Individual report for individual <b>" << a->indiv_id << "</b>";
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

      std::cout << "<a href=\"pbrowse.cgi?q=v&" << a->print_form_value("proj")
		<< "&gene=" << a->genename
		<< "&meta=" << a->mf_print()
		<< "&masks=" << a->msk_print() 
		<< "&pheno=" << a->phenotype_name 
		<< "&val=" << vars.var(i).coordinate()	
		<< "&regs=" << a->reg_list_url
		<< "\">view</a>";

      std::cout << "</td>";

      std::cout << "<td>" << vars.var(i).chromosome() << "</td>" 
	   << "<td>" << vars.var(i).position() << "</td>";

      // Print name, with link to dbSNP if appropriate
      std::cout << "<td>" << rs_link( vars.var(i).name() ) << "</td>";
           
      std::cout << "<td>" << vars.var(i).pp_reference() << "/" << vars.var(i).pp_alternate() << "</td>" ;


      // Optional meta-information?
      
      for (int m=0; m< a->mf.size(); m++)
	std::cout << "<td>" << vars.var(i).print_meta( a->mf[m] , "<br>" ) << "</td>";
      

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
	  ++k;
	}
      std::cout << "</tr>";	  

    }
  
  std::cout << "</table>";

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
  std::vector<std::string> genes = a->g->locdb.fetch_names( a->loc_set , true );
  
  // these are sorted, so only display unique ones
  std::string lastgene = "";
  for (int g=0; g<genes.size(); g++)
    {      
      if ( genes[g] != lastgene )
	{
	  std::cout << "<a href=\"pbrowse.cgi?q=g&" << a->print_form_value("proj")
		    << "&gene="<< genes[g]
		    << "&meta=" << a->mf_print()
		    << "&masks=" << a->msk_print() 
		    << "&pheno=" << a->phenotype_name 
		    << "&regs=" << a->reg_list_url	 
		    << "\">" << genes[g] << "</a><br>";
	  lastgene = genes[g];
	}
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
	    << a->g->summary() 
	    << "</pre></p>";
  std::cout << "</table></body></html>";
  exit(0);
}
