#include "psb.h"

#include "cgi.h"

#include <iostream>

using namespace std;
using namespace Helper;
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
  
  string project_path = "";
  string loc_set = "refseq";
  string genename = "";
  QType q = Q_ERROR;
  string var_value = "";
  string ind_value = "";
  string chr_code = "";
  int chr = 0, bp1 = 0, bp2 = 0;
  string pheno = "";
  bool from_top = false;
  string reg_list = "";


  // Auxilliary information to send to variant printing functions

  Aux a;
  
  
  //
  // Start HTML form and any output
  //

  cout << "Content-type: text/html\n\n"
       << "<html><head><title>PLINK/SEQ Browser</title>"
       << "</head><body>";
  
  if ( cgi ) 
    {
      for (i=0; cgivars[i]; i+= 2)
	{
	  
	  string str = cgivars[i];
	  
	  if ( str == "q" )
	    {
	      string s = cgivars[i+1];
	      
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
	  
	  if ( str == "gene" ) 
	    {
	      genename = cgivars[i+1];
	      Helper::str2upper( genename );
	    }
	  
	  if ( str == "ind" ) 
	    ind_value = cgivars[i+1];
	  
	  if ( str == "proj" ) 
	    {
	      project_path = cgivars[i+1];
	      a.add_form_value("proj", project_path );
	    }
	  
	  if ( str == "loc" ) loc_set = cgivars[i+1];
	  
	  if ( str == "regs" )
	    reg_list = cgivars[i+1];
	  
	  if ( str == "masks" )
	    a.msk = parse( cgivars[i+1] , " " );
	  
	  if ( str == "meta" ) 
	    a.mf = parse( cgivars[i+1] , " ," );
	  
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


  if ( !cgi ) 
    {
      q = Q_GENE;
      genename = "NLRP11";
      project_path = "../plinkseq/browser/ceu.exome";
      loc_set = "refseq";
    }
  
  
  //
  // Start page 
  //
  
  cout << "<table width=100% CELLPADDING=0 CELLSPACING=0>"
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
      
      if ( ! fileExists( project_path ) )
	{
	  cout << "File [ " 
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
  
  cout << "<form name=\"myform\" action=\"pbrowse.cgi\" method=\"GET\"> ";
  
  cout << "<input type=\"hidden\" NAME=\"proj\" "
       << " value=\"" << project_path << "\"> ";
  
  cout << "<table width=100%><tr>";

  // Hidden value to indicate query type (needed?)

  cout << "<input type=\"hidden\" name=\"q\" value=\"";
  if ( q == Q_GENE || q == Q_ERROR ) cout << "g\">";
  else if ( q == Q_VARIANT ) cout << "v\">";
  else if ( q == Q_INDIV ) cout << "i\">";
  
  
  //
  // Panels
  //


  cout << "<td width=20% valign=top>";


  //
  // Query 
  //

  cout << "<p><b>Gene ID</b> (symbol or NM_012345) ";
  cout << "(<a href=\"pbrowse.cgi?q=glist&" << a.print_form_value("proj")
       << "&meta=" << a.mf_print()
       << "&masks=" << a.msk_print() 
       << "&pheno=" << pheno
       << "\">list</a>)";			  
  
  cout << "<br><input type=\"text\" size=\"15\" name=\"gene\"";
  if ( genename != "" )
    cout << " value=\"" << genename << "\"";
  cout << ">";
  cout << "</p> ";
  
  // Gene-set
  cout << "<p>Gene set ";
  cout << "(<a href=\"pbrowse.cgi?q=lslist&" << a.print_form_value("proj")
       << "\">list</a>)";			  

  cout << "<br><input type=\"text\" size=\"15\" ";  
  if ( loc_set != "" ) 
    cout << " value=\"" << loc_set << "\"";
  cout << " name=\"loc\"></p>";
    

  //
  // Submit buttons
  //
  
  cout << " <br> <input type=\"submit\" name=\"getgene\" value=\"Fetch\"> "
       << " <input type=\"reset\" value=\"Reset\"> ";


  cout << "</td><td width=5% valign=top> &nbsp ";
  cout << "</td><td width=30% valign=top>";


  //
  // Optional meta-fields
  //
  
  cout << "<p>Optional variant meta-fields ";
  cout << " (<a href=\"pbrowse.cgi?q=mflist&" << a.print_form_value("proj")
       << "&meta=" << a.mf_print()
       << "&masks=" << a.msk_print()
       << "&pheno=" << pheno
       << "\">list</a>) ";			  

  cout << "<br><input type=\"text\" size=\"30\" name=\"meta\"";
  
  if ( a.mf.size() != 0 ) 
    {
      cout << " value=\"";
      for (int m=0; m< a.mf.size(); m++) 
	{
	  if ( m>0 ) cout << " ";
	  cout << a.mf[m];
	}
      cout << "\"";
    }
  cout << ">";
  cout << "</p>";


  //
  // Phenotype
  //

  cout << "<p>Optional case/control phenotype";

  cout << " (<a href=\"pbrowse.cgi?q=plist&" << a.print_form_value("proj")
       << "&meta=" << a.mf_print()
       << "&masks=" << a.msk_print() 
       << "&pheno=" << pheno
       << "\">list</a>) ";			  

  cout << "<br><input type=\"text\" size=\"30\" name=\"pheno\"";
  
  if ( pheno != "" ) 
    {
      cout << " value=\"";
      cout << pheno << "\"";
    }
  
  cout << "</p>";


  //
  // Mask specification
  //
  
  std::cout << "<p>Optional mask specification";
  std::cout << "<br><input type=\"text\" size=\"30\" name=\"masks\"";
  
  if ( a.msk.size() != 0 ) 
    {
      std::cout << " value=\"";
      for (int m=0; m< a.msk.size(); m++) 
	{
	  if ( m>0 ) cout << " ";
	  std::cout << a.msk[m];
	}
      std::cout << "\"";
    }
  std::cout << ">";
  std::cout << "</p>";
  

  //
  // Third column
  //

  cout << "</td><td width=5% valign=top> &nbsp ";
  cout << "</td><td width=30% valign=top>";
  

  //
  // Region list
  // 

  cout << "<p>Additional regions and genes<br>";
  cout << "<textarea ";
  cout << "name=\"regs\" rows=\"6\" cols=\"30\">"
       << "</textarea></p>";


  cout << "</td></tr></table>";


  cout << " </form> ";
  cout << "<hr> ";


  //
  // Draw query 
  //

  if ( q == Q_ERROR )
    {
      cout << "Problem processing input...";
      cout << "</body></html>";
      exit(0);
    }
 
  
  if ( ! fileExists( project_path ) )
    {
      cout << "File [ " << project_path << " ] could not be found "
           << "</BODY></HTML>";
      exit(0);
    }
  
  
  if ( q == Q_GENE && genename == "" && reg_list == "" )
    {
      cout << "No gene specified...</BODY></HTML>";
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
  
  string main_gene = genename;

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

  if ( reg_list != "" )
    {

      vector<string> tok = whitespace( reg_list );
      
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
	cout << "Found " << trans.size() 
	     << " transcript(s) matching gene name <b>" 
	     << main_gene << "</b></p>";
      
      cout << "<pre><font size=-1>";
      
      for (int r=0; r<trans.size(); r++)
	{
	  
	  cout << trans[r].altname << "   " ;
	  cout << "<a href=\"pbrowse.cgi?q=g&" << a.print_form_value("proj")
	       << "&gene="<< trans[r].name	       
	       << "&meta=" << a.mf_print()
	       << "&masks=" << a.msk_print() 
	       << "&pheno=" << pheno
	       << "\">" << trans[r].name <<"</a>  "
	       << chrCode(trans[r].start.chromosome()) << ":" 
	       << trans[r].start.position() << ".."
	       << trans[r].stop.position() ;
	  cout <<"<br>";
	  if ( !cgi ) cout << "\n";
	}
      
      cout << "</font></pre>";
            

      //
      // Now, 'genename' will be original search term
      //
      
      if ( ! a.multi_transcripts )
	{
	  
	  Region reg = g.locdb.get_region( loc_set , genename ) ;
	  
	  if ( reg.start.chromosome() == 0 && ! a.extended_search )
	    {
	      cout << "Could not find gene <b> " << genename << "</b> in <b>" << loc_set << "</b> list "
		   << "</BODY></HTML>";
	      exit(0);      
	    }
	  
	  a.region = reg;
	  string g_chr_code = chrCode( reg.start.chromosome() );
	  int g_chr = reg.start.chromosome();
	  int g_bp1 = reg.start.position();
	  int g_bp2 = reg.stop.position();
	  int n_exons = reg.subregion.size();
	  
	  // Query to get all other transcripts in region
	  
	  set<Region> others = g.locdb.get_regions( g.locdb.lookup_group_id( loc_set ) , 
						    g_chr , g_bp1 , g_bp2 );
	  
	  cout << "<b>" << genename << "</b> location : " 
	       << g_chr_code << ":" << g_bp1 << ".." << g_bp2 
	       << "  (view in the "
	       << "<a href=\"http://genome.ucsc.edu/cgi-bin/hgTracks?org=Human&db=hg18&position=" 
	       << g_chr_code << ":" << g_bp1 << ".." << g_bp2     
	       << "\" target=\"_blank\">UCSC genome browser</a>)</p>";
	  
	  
	  cout << "<pre><font size=-1>";
	  
	  double genomic_kb = (g_bp2 - g_bp1 + 1 ) / 1000.0;
	  double coding_kb = 0;
	  for (int s=0; s<reg.subregion.size(); s++)
	    coding_kb += ( reg.subregion[s].stop.position() - reg.subregion[s].start.position() + 1 );
	  coding_kb /= 1000.0;
	  
	  cout << "   Genomic kb = " << genomic_kb 
	       << "; coding kb = " << coding_kb 
	       << " in " << n_exons << " exons</p>";
	  
	  // List exons, and other transcripts that overlap this one
      
	  for (int s=0; s<reg.subregion.size(); s++)
	    {
	      
	      cout << "   Exon" << sw(s+1,3) ;
	      cout << sw( int2str(reg.subregion[s].start.position()) 
			  + " .. " + int2str(reg.subregion[s].stop.position()) , 16 );
	      cout << " | " ;
	      cout << sw( ( reg.subregion[s].stop.position() 
			    - reg.subregion[s].start.position() + 1 )  , 5 );
	      cout << "bp | ";
	      
	      // other exons this exon overlaps with? 
	      set<Region>::iterator t = others.begin();
	  
	      while ( t != others.end() )
		{
		  if ( t->name != reg.name )
		    for (int u=0; u< t->subregion.size(); u++)
		      {
			if ( reg.subregion[s].overlaps( t->subregion[u] ) )
			  {
			    cout << "<a href=\"pbrowse.cgi?q=g&" << a.print_form_value("proj")
				 << "&gene="<< t->name
				 << "&meta=" << a.mf_print()
				 << "&masks=" << a.msk_print() 			
				 << "&pheno=" << pheno
				 << "\">"<<t->name<<"("<<t->altname<<")</a> ";			  
			  }
		      }
		  ++t;
		}
	      cout << "<br>";
	    }
	  
	  cout << "</font></pre>";
      
	}

      cout << "<hr>";      
    }


  if ( a.extended_search ) 
    {
      cout << "<p>Additionally ";
      if ( genename != "" ) cout << "requiring ";
      else cout << "including ";
      cout << "the following regions, and including the following genes</p>";
      
      if ( a.regions.size() > 0 ) 
	{
	  cout << "<pre><b>Regions:</b><br>";
	  for (int i=0; i<a.regions.size(); i++)
	    cout << "   " << chrCode( a.regions[i].chromosome()) << ":"
		 << a.regions[i].start.position() << ".."
		 << a.regions[i].stop.position() << "<br>";
	  cout << "</pre>";
	}
      
      if ( a.other_genes.size() > 0 ) 
	{
	  vector<string> cc = a.other_genes;
	  a.other_genes.clear();
	  
	  cout << "<pre><b>Additional genes:</b><br>";
	  for (int i=0; i< cc.size(); i++)
	    {
	      
	      std::set<std::string> trans = 
		g.locdb.targetted_lookup_alias( cc[i] , "symbol" , "refseq" ) ;   
	      
	      // if no aliases match, assume this is a refseq transcript name
	      if ( trans.size() == 0 ) 
		{
		  a.other_genes.push_back( cc[i] );
		  cout << "   " << sw( cc[i] , 10 ) << "<br>";
		}
	      else // add all matching aliases
		{
		  cout << "   " << sw( cc[i] , 10 ) << " : ";
		  
		  std::set<std::string>::iterator ii = trans.begin();
		  while ( ii != trans.end() ) 
		    {
		      a.other_genes.push_back( *ii );
		      cout << *ii << " ";
		      ++ii;
		    }		 
		  cout << "<br>";
		}
	    }
	  cout << "</pre>";
	}

      cout << "<hr>";
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
      // map<int,string> table_row;
      
      cout << a.headers << "</p>";
      // LINKS CHR POS NAME ALT/REF FileID QUAL INFO(filter) SAMPLE_CNT VMETA(compressed)
      
      cout << "<table border=1>"
	   << "<tr><th>#</th><th>Indiv</th><th>Chr</th><th>Pos</th>";
      
      if (  ! ( a.extended_search || a.multi_transcripts ) )
	cout << "<th>Exon</th>";
      
      cout << "<th>ID</th><th>Ref/Alt</th>"
	   << "<th># samples</th>"      
	   << "<th>Filter</th>";      
      
      
      // Optional case/control counts?
      
      if ( a.show_phenotype )
	cout << "<th>C/C count</th>";
      
  
      // Optional meta-fields?
      
      for (int m=0; m< a.mf.size(); m++)
	cout << "<th>" << a.mf[m] << "</th>";

      if ( a.extended_search ) 
	cout << "<th>Locus</th>";
      
     
      // End of header row
      
      cout << "</tr>";
      
      // Rows

      map<int,string>::iterator k= a.table_row.begin();
      while ( k != a.table_row.end() )
	{
	  cout << k->second;
	  ++k;
	}
      cout << "</table>";

    }  


  //
  // Query for a single variant: no need to call any database iteraton functions
  //
  
  if ( q == Q_VARIANT )
    {
      
      cout << "Back to "
	   << "<a href=\"pbrowse.cgi?q=g&" << a.print_form_value("proj")
	   << "&gene="<< genename  
	   << "&meta=" << a.mf_print()
	   << "&masks=" << a.msk_print() 
	   << "&pheno=" << a.phenotype_name 
	   << "\">gene report</a>"
	   << " for <b>" << genename << "</b>";
 
       cout << "<hr>";
      


       // Expect here we have variant ID
       
       // note == first param is actually redundant      
       
       bool okay = true;
       Region r( var_value , okay  );
       if ( ! okay ) 
	 {
	   cout << "Problem processing region code [ " << var_value << " ]</p>";
	   cout << "</body></html>";
	   exit(0);
	 }



       //
       // Get variant, based on physical position
       //
       
       Variant var = g.vardb.fetch( r.chromosome() , r.start.position() );

       
       //
       // Display
       //


       cout << "<table width=90%><tr><td width=40% valign=top>";

       cout << "<h3><font color=\"blue\">Variant information</font></h3>";
         
       cout << "<table border=1><tr><th>Field</th><th>Value</th></tr>";
       cout << "<tr><td>Name</td><td>" << rs_link( var.name() ) << "</td></tr>";
       cout << "<tr><td>Chromosome</td><td>" << var.chromosome() << "</td></tr>";
       cout << "<tr><td>Position</td><td>" << var.position() << "</td></tr>";
       cout << "<tr><td>Reference allele</td><td>" << var.reference() << "</td></tr>";
       cout << "<tr><td>Alternate allele(s)</td><td>" << var.alternate() << "</td></tr>";
       cout << "<tr><td>Samples</td><td>" << var.n_samples() << "</td></tr>";

       vector<string> keys = var.meta.keys();
       for (int m=0; m<keys.size(); m++)
	 cout << "<tr><td>" << keys[m] << "</td><td>" << var.meta.print( keys[m] ) << "</td></tr>";	  
       cout << "</table>";

       cout << "</p>";
       
       // Consensus meta-information
       
       SampleVariant & sample = var.consensus;
       
       cout << "<p>Consensus</p>";	   
       cout << "<table border=1><tr><th>Field</th><th>Value</th></tr>";	   
       cout << "<tr><td>Reference allele</td><td>" << sample.reference() << "</td></tr>";
       cout << "<tr><td>Alternate allele(s)</td><td>" << sample.alternate() << "</td></tr>";

       if ( sample.quality() < 0 ) 
	 cout << "<tr><td>Quality</td><td>" << "NA" << "</td></tr>";
       else
	 cout << "<tr><td>Quality</td><td>" << sample.quality() << "</td></tr>";

       cout << "<tr><td>Filter</td><td>" << sample.filter() << "</td></tr>";
//     cout << "<tr><td>Information</td><td>" << sample.info() << "</td></tr>";       
       vector<string> keys1 = sample.meta.keys();
       for (int m=0; m< keys1.size(); m++)
	 cout << "<tr><td>" << keys1[m] 
	      << "</td><td>" 
	      << sample.meta.print( keys1[m] ) 
	      << "</td></tr>";
       
	   cout << "</table>";

       // SampleVariant Tables
       
       var.set_first_sample();
       
       while ( 1 ) 
	 {

	   SampleVariant & sample = var.sample();

	   cout << "<p>Sample " << g.vardb.file_tag( sample.fileset() ) << "</p>";
	   
	   cout << "<table border=1><tr><th>Field</th><th>Value</th></tr>";
	   
	   cout << "<tr><td>Reference allele</td><td>" << sample.reference() << "</td></tr>";
	   cout << "<tr><td>Alternate allele(s)</td><td>" << sample.alternate() << "</td></tr>";

	   if ( sample.quality() < 0 ) 
	     cout << "<tr><td>Quality</td><td>" << "NA" << "</td></tr>";
	   else
	     cout << "<tr><td>Quality</td><td>" << sample.quality() << "</td></tr>";
	   
	   cout << "<tr><td>Filter</td><td>" << sample.filter() << "</td></tr>";
//	   cout << "<tr><td>Information</td><td>" << sample.info() << "</td></tr>";
	   
	   vector<string> keys = sample.meta.keys();
	   for (int m=0; m< keys.size(); m++)
	     cout << "<tr><td>" << keys[m] << "</td><td>" << sample.meta.print( keys[m] ) << "</td></tr>";	  
	   
	   cout << "</table>";
		   
	   if ( ! var.next_sample() ) break;
	 }
       

       //
       // Genotype information
       //
       
       cout << "</td><td valign=top>";
       
       cout << "<h3><font color=\"blue\">Individual genotypes</font></h3>";
       
       
       // Get list of IDs for consensus set
       vector<string> id_list = var.ind_id();
       
       set<string> gmeta;
       set<string> gRef;
       set<string> gHet;
       set<string> gHom;
       set<string> gMis;

       SampleVariant & con = var.consensus;
       
       
       if ( id_list.size() != con.calls.size() ) 
	 Helper::halt("internal problem in indmap/var-con sizes");
       
       for (int i=0; i < id_list.size(); i++)
	 {
	   std::vector<std::string> k = con.calls.genotype(i).meta.keys();
	   for (int j=0; j<k.size(); j++) gmeta.insert(k[j]);
	   Genotype & g = con.calls.genotype(i);
	   if ( g.null() ) gMis.insert( id_list[i] ) ;
	   else 
	     {
	       int ac = g.allele_count( &var );
	       if ( ac == 0 ) gRef.insert( id_list[i] );
	       else if ( ac == 1 ) gHet.insert( id_list[i] );
	       else if ( ac == 2 ) gHom.insert( id_list[i] );
	     }	  
	 }
       
       cout << "<table border=1><tr><th>Individual ID</th>";
      
       cout << "<th nowrap>Sample #</th>";
       
      if ( a.show_phenotype ) 
	cout << "<th>Phenotype</th>";
            
      std::cout << "<th>Genotype</th>";
      
      if ( var.flat() )
	{
	  std::set<std::string>::iterator i = gmeta.begin();
	  while ( i != gmeta.end() )
	    {
	      cout << "<th>" << *i << "</th>";
	      ++i;
	    }
	}
      else
	cout << "<th>Sample/genotype meta-information</th>";
      
      cout << "</tr>";


      // Done header, now list individuals, in genotypic order
      
      for (int gt=0; gt<4; gt++)
	{
	  for (int i=0; i<id_list.size(); i++)
	    {
	      
	      if ( gt == 3 && ! con.calls.genotype(i).null() ) continue;
	      else if ( con.calls.genotype(i).allele_count( &var ) != 2-gt ) continue;
	      
	
	      // Include link to individual-report

	      cout << "<tr><td>" 
		   << "<a href=\"pbrowse.cgi?q=i&" << a.print_form_value("proj")
		   << "&gene="<< genename  
		   << "&meta=" << a.mf_print()
		   << "&masks=" << a.msk_print() 
		   << "&ind=" << id_list[i] 
		   << "&pheno=" << a.phenotype_name 
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
		      if ( person->affected() == CASE )
			cout << "<td>CASE</td>";
		      else if ( person->affected() == CONTROL )
			cout << "<td>CONTROL</td>";
		      else
			cout << "<td>MISSING</td>";
		    }
		  else
		    cout << "<td>n/a</td>";
		}

	      if ( gt == 0 ) 
		cout << "<td bgcolor=\"red\">";
	      else if ( gt == 1 ) 
		cout << "<td bgcolor=\"yellow\">";
	      else if ( gt == 2 )
		cout << "<td bgcolor=\"lightgreen\">";
	      else 
		cout << "<td bgcolor=\"gray\">";
	      
	      if ( var.flat() )
		cout << var.consensus.label( con.calls.genotype(i) ) ;
	      else
		std::cout << var.label(i) ;
	      
	      std::cout << "</td>";
	      
	      if ( var.flat() )
		{
		  set<string>::iterator k = gmeta.begin();
		  while ( k != gmeta.end() )
		    {
		      if ( con.calls.genotype(i).meta.has_field( *k ) )
			{
			  meta_name_t t = *k;
			  cout << "<td>" << con.calls.genotype(i).meta.print( t ) << "</td>";
			}
		      ++k;
		    }
		}
	      else
		std::cout << "<td>" << var.gmeta_label(i) << "</td>";
	      
	      cout << "</tr>";	  
	    }
	}
      
      cout << "</table>";
             

      //
      // End of table
      //
      

      cout << "</td></tr></table>";
      
    }


  //
  //
  //

  if ( q == Q_INDIV )
    {

      a.genename = genename;

      cout << "Back to "
	   << "<a href=\"pbrowse.cgi?q=g&" << a.print_form_value("proj")
	   << "&gene="<< genename  
	   << "&meta=" << a.mf_print()
	   << "&masks=" << a.msk_print() 
	   << "&pheno=" << a.phenotype_name 
	   << "\">gene report</a>"
	   << " for <b>" << genename << "</b>";
      
      cout << "<hr>";
      

      // 1) Report all phenotypic information for this individual
      
      if ( a.show_phenotype ) 
	{

	  // Expect here we have a individual ID
	  
	  Individual * person = g.indmap.ind( ind_value );
	  
	  if ( ! person ) 
	    {
	      cout << "Individual [ " << ind_value << " ] not found...</p>";
	      cout << "</BODY></HTML>";
	      exit(0);	  
	    }
	  
	  cout << "Phenotype for variable <em>" << pheno << "</em> is [ " ;
	  if ( person->affected() == CASE )
	    cout << "CASE";
	  else if ( person->affected() == CONTROL )
	    cout << "CONTROL";
	  else
	    cout << "MISSING";
	  cout << " ]</p>";
	}

      
      // 2) Report all variants for gene for ths individual
      
      a.indiv_id = ind_value;


      // Crate Mask, ( use existing mask, no?)

//       std::string mstr;
//       for (int i = 0 ; i < a.msk.size(); i++)
// 	mstr += ( i ? " " : "" ) + a.msk[i];
//       Mask m( mstr );


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
      if ( a.regions.size() > 0 )
	{
	  for (int r = 0 ; r < a.regions.size() ; r++) 
	    m.require_reg( a.regions[r] );
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

  cout << "</body></html>";
  

  exit(0);
  
}


void ExomeBrowser::f_display(Variant & var, void *p)
{

  // Aux. data structure
  
  Aux * a = (Aux*)p;

  stringstream o1;
  
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
  
  o1 << "<td>" << var.reference() << "/" << var.alternate() << "</td>";

  o1 << "<td>" << var.n_samples() << "</td>";
  
  o1 << "<td>" << var.print_meta_filter( "<br>" ) << "</td>";

  SampleVariant & con = var.consensus;


  // Allele count?
  
  if ( a->show_phenotype )
    {
      int case_n = 0 , control_n = 0;
      for (int j=0; j<con.calls.size(); j++)
	{
	  affType aff = var.ind( j )->affected();
	  if ( con(j).nonreference() )
	    {
	      if ( aff == CASE ) ++case_n;
	      else if ( aff == CONTROL ) ++control_n;
	    }
	}
      o1 << "<td>" << case_n << "/" << control_n << "</td>";
    }
  
  // Optional meta-information?
  
  for (int m=0; m < a->mf.size(); m++)
    o1 << "<td nowrap>" << var.print_meta( a->mf[m] , "<br>" ) << "</td>";
  
  
  // Appends? 
  
  if ( a->extended_search ) 
    {
      if ( var.meta.has_field( PLINKSeq::META_LSET() ) )
	{
	  string s = PLINKSeq::META_LSET();
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
  
  cout << "Individual report for individual <b>" << a->indiv_id << "</b> for gene " << vars.name() << "</b></p>";
  
  cout << vars.size() << " variants in total<br>";
  
  if ( vars.size() == 0 ) return;
  
  //vars.prepare();

  // LINKS CHR POS NAME ALT/REF  GENO  GENE_META...
  
  cout << "<table border=1>"
       << "<tr><th>#</th><th>Indiv</th><th>Chr</th><th>Pos</th><th>ID</th><th>Ref/Alt</th>";

  // Optional meta-fields?
  
  for (int m=0; m< a->mf.size(); m++)
    cout << "<th>" << a->mf[m] << "</th>";
  
  cout << "<th>Genotype</th>";
  
  //
  // Get genotype meta-fields
  //

  set<string> gmeta;
  
  for (int i=0; i < vars.size(); i++)
    {            
      SampleVariant & con = vars(i).consensus;            
      vector<string> k = con.calls.genotype( ni ).meta.keys();
      for (int z=0; z<k.size(); z++) 
	gmeta.insert(k[z]);      
    }
  

  // We now have a list of all relevant genotypic meta fields, and we
  // know which position in each file this individual is sitting in
  
  set<string>::iterator i = gmeta.begin();
  while ( i != gmeta.end() )
    {
      cout << "<th>" << *i << "</th>";
      ++i;
    }
  
  // End of header row
  
  cout << "</tr>";
  


  //
  // Per-variant summary
  //

  for (int i=0; i<vars.size(); i++)
    {
      
      SampleVariant & con = vars(i).consensus;

      cout << "<tr>";

      // counter
      cout << "<td>" << i+1 << "</td>";

      // Links per variant
      cout << "<td>";

      // call pbrowse.cgi, but with all arguments explicitly formed.      

      cout << "<a href=\"pbrowse.cgi?q=v&" << a->print_form_value("proj")
	   << "&gene=" << a->genename
	   << "&meta=" << a->mf_print()
	   << "&masks=" << a->msk_print() 
	   << "&pheno=" << a->phenotype_name 
       	   << "&val=" << vars.var(i).coordinate()	
	   << "\">view</a>";

      cout << "</td>";

      cout << "<td>" << vars.var(i).chromosome() << "</td>" 
	   << "<td>" << vars.var(i).position() << "</td>";

      // Print name, with link to dbSNP if appropriate
      cout << "<td>" << rs_link( vars.var(i).name() ) << "</td>";
           
      cout << "<td>" << vars.var(i).reference() << "/" << vars.var(i).alternate() << "</td>" ;


      // Optional meta-information?
      
      for (int m=0; m< a->mf.size(); m++)
	cout << "<td>" << vars.var(i).print_meta( a->mf[m] , "<br>" ) << "</td>";
      

      // Genotype; genotype meta-fields;

      int gt = con(ni).null() ? 3 : 
	2 - con(ni).allele_count( &vars.var(i) );
      
      if ( gt == 0 ) 
	cout << "<td bgcolor=\"red\">";
      else if ( gt == 1 ) 
	cout << "<td bgcolor=\"yellow\">";
      else if ( gt == 2 )
	cout << "<td bgcolor=\"lightgreen\">";
      else 
	cout << "<td bgcolor=\"gray\">";
            
      cout << vars.var(i).consensus.label( con( ni ) )
	   << "</td>";
      
      set<string>::iterator k = gmeta.begin();
      while ( k != gmeta.end() )
	{
	  if ( vars.var(i).consensus.calls.genotype(ni).meta.hasField( *k ) )
	    {
	      meta_name_t t = *k;
	      cout << "<td>" << vars.var(i).consensus.calls.genotype(ni).meta.print( t ) << "</td>";
	    }
	  ++k;
	}
      cout << "</tr>";	  

    }
  
  cout << "</table>";

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



string ExomeBrowser::rs_link(const string & label )
{  
  if ( label == "0" ) return "n/a";  // fix old missing = 0 code      
  if ( label.size() > 2 && label.substr(0,2) == "rs" ) 
    return "<a href=\"http://www.ncbi.nlm.nih.gov/sites/entrez?db=snp&cmd=search&term=" 
      + label + "\" target=\"_blank\">" + label + "</a>"; 
  return label;
}



void ExomeBrowser::make_gene_list(Aux * a)
{

  // true param indicates to fetch alternate names
  vector<string> genes = a->g->locdb.fetch_names( a->loc_set , true );
  
  // these are sorted, so only display unique ones
  string lastgene = "";
  for (int g=0; g<genes.size(); g++)
    {      
      if ( genes[g] != lastgene )
	{
	  cout << "<a href=\"pbrowse.cgi?q=g&" << a->print_form_value("proj")
	       << "&gene="<< genes[g]
	       << "&meta=" << a->mf_print()
	       << "&masks=" << a->msk_print() 
	       << "&pheno=" << a->phenotype_name 
	       << "\">" << genes[g] << "</a><br>";
	  lastgene = genes[g];
	}
    }
  cout << "</body></html>";
  exit(0);
}


void ExomeBrowser::make_phe_list(Aux * a)
{
  map<string,vector<string> > m = a->g->inddb.fetch_phenotype_info();
  map<string,vector<string> >::iterator i = m.begin();
  cout << "<h3>Available phenotypes</h3>";
  cout << "<table border=1><tr>"
       << "<th align=left>Phenotype</th>"
       << "<th align=left>Type</th>"
       << "<th align=left>Description</th></tr>";
  while ( i != m.end() )
    {
      cout << "<tr>"
	   << "<td>" << i->first << "</td>"
	   << "<td>" << i->second[0] << "</td>"
	   << "<td>" << i->second[1] << "</td>"
	   << "</tr>";
      ++i; 
    }
  cout << "</table></body></html>";
  exit(0);
}



void ExomeBrowser::make_mf_list(Aux * a)
{
  
  map<int,string> f = a->g->vardb.fetch_files();
  map<int,string>::iterator i = f.begin();
  while ( i != f.end() )
    {
      cout << "<h3><font color=blue>" << i->second << "(" << i->first << ")</font></h3>";
      // Meta-types
      cout << "<h4>Variant meta-information fields</h4>";
      cout << "<table border=1><th align=left>Meta-field</th>"
 	   << "<th align=left>Type</th>"
 	   << "<th align=left>Length</th>"
 	   << "<th align=left>Group</th>"
 	   << "<th align=left>Description</th>"
 	   << "</tr>";
      
      vector<map<string,string> > m = a->g->vardb.fetch_metatypes( i->first );
      for (int j=0; j<m.size(); j++)
	{
	  cout << "<tr><td>";
	  
	  if ( MetaMeta::static_variant( m[j]["NAME"] ) )
	    cout << "<font color=\"green\"> " << m[j]["NAME"] << "</font></td><td>";
	  else
	    cout << m[j]["NAME"] << "</td><td>";
	  
	  cout << m[j]["TYPE"] << "</td><td>"
	       << m[j]["NUM"] << "</td><td>"
	       << m[j]["GRP"] << "</td><td>"
	       << m[j]["DESC"] << "</td></tr>";
	}
      cout << "</table>";

      // headers
      cout << "<h4>Header rows in VCF</h4>";
      cout << "<table border=1>"
 	   << "<th align=left>Key</th>"
 	   << "<th align=left>Value</th></tr>";
      
     vector<map<string,string> > h = a->g->vardb.fetch_headers( i->first );
      for (int j=0; j<h.size(); j++)
	{
	  // do not print out first line of VCF #CHROM...
	  string ky =  h[j]["KEY"];
	  if ( ky.size() > 4 && ky.substr(0,4) == "#CHR" )
	    continue;
	  cout << "<tr><td>" 
	       << ky << "</td><td>"
	       << h[j]["VALUE"] << "</td></tr>";
	}
      cout << "</table>";      
      cout <<"<hr>";
      ++i;
    }
  cout << "</table></body></html>";
  exit(0);
}


void ExomeBrowser::make_locset_list(Aux * a)
{
  
  set<GroupInfo> g = a->g->locdb.group_information();
  set<GroupInfo>::iterator i = g.begin();

  cout << "<h3>Available locus sets</h3>";
  cout << "<table border=1><tr>"
       << "<th align=left>Name</th>"
       << "<th align=left>Number of entries</th>"
       << "<th align=left>Decscription</th>"
       << "</tr>";
  
  while ( i != g.end() )
    {
      cout <<"<tr>"
	   << "<td>" << i->name << "</td>"
	   << "<td>" << a->g->locdb.count( i->idx )<< "</td>"
	   << "<td>" << i->description << "</td>"
	   << "</tr>";
      ++i; 
    }
  cout << "</table></body></html>";
  exit(0);
}

void ExomeBrowser::make_proj_summary(Aux * a)
{
  std::cout << "<h3>Project Summary</h3>"
	    << "<p><pre>"
	    << a->g->summary() 
	    << "</pre></p>";
  cout << "</table></body></html>";
  exit(0);
}
