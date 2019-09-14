<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html lang="en">
  <head>
    <meta http-equiv="X-UA-Compatible" content="IE=edge"/>
      
      <link rel="stylesheet" type="text/css" href="/scripts/shadowbox/shadowbox.css"/>
      	<script type="text/javascript" src="/includes_content/nextgen/scripts/jquery/jquery-latest.js"></script>
      <!-- START OF GLOBAL NAV -->

  <link rel="stylesheet" href="/matlabcentral/css/sitewide.css" type="text/css">
  <link rel="stylesheet" href="/matlabcentral/css/mlc.css" type="text/css">
  <!--[if lt IE 7]>
  <link href="/matlabcentral/css/ie6down.css" type="text/css" rel="stylesheet">
  <![endif]-->

      
      
      
      <meta http-equiv="content-type" content="text/html; charset=UTF-8">
<meta name="keywords" content="file exchange, matlab answers, newsgroup access, link exchange, matlab blog, matlab central, simulink blog, matlab community, matlab and simulink community">
<meta name="description" content="File exchange, MATLAB Answers, newsgroup access, Links, and Blogs for the MATLAB &amp; Simulink user community">
<link rel="stylesheet" href="/matlabcentral/css/fileexchange.css" type="text/css">
<link rel="stylesheet" type="text/css" media="print" href="/matlabcentral/css/print.css" />
<title> File Exchange - MATLAB Central</title>
<link rel="canonical" href="https://uk.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation/content/DERIVESTsuite/derivest.m" />
<script src="/matlabcentral/fileexchange/assets/application.js" type="text/javascript"></script>
<link href="/matlabcentral/fileexchange/assets/application.css" media="screen" rel="stylesheet" type="text/css" />
<link href="/matlabcentral/fileexchange/assets/profile_link/application.css" media="all" rel="stylesheet" type="text/css" />
<link href="/includes_content/responsive/fonts/mw_font.css" media="all" rel="stylesheet" type="text/css" />
<script src="/matlabcentral/fileexchange/assets/profile_link/application.js" type="text/javascript"></script>


<link rel="search" type="application/opensearchdescription+xml" title="Search File Exchange" href="/matlabcentral/fileexchange/search.xml" />


<!-- BEGIN Adobe DTM -->
<script src="/scripts/dtm/d0cc0600946eb3957f703b9fe43c3590597a8c2c/satelliteLib-e8d23c2e444abadc572df06537e2def59c01db09.js"></script>
<!-- END Adobe DTM -->  </head>
    <body>
      <div id="header" class="site6-header">
  <div class="wrapper">
  <!--put nothing in left div - only 11px wide shadow -->
    <div class="main">
        <div id="headertools">


<script language="JavaScript1.3" type="text/javascript">

function submitForm(query){

        choice = document.forms['searchForm'].elements['search_submit'].value;
        
        if (choice == "entire1" || choice == "contest" || choice == "matlabcentral" || choice == "blogs"){
        
           var newElem = document.createElement("input");
           newElem.type = "hidden";
           newElem.name = "q";
           newElem.value = query.value;
           document.forms['searchForm'].appendChild(newElem);
              
           submit_action = '/search/site_search.html';
        }
        
        switch(choice){
           case "matlabcentral":
              var newElem = document.createElement("input");
              newElem.type = "hidden";
              newElem.name = "c[]";
              newElem.value = "matlabcentral";
              document.forms['searchForm'].appendChild(newElem);
        
              selected_index = 0;
              break
           case "fileexchange":
              var newElem = document.createElement("input");
              newElem.type = "hidden";
              newElem.name = "term";
              newElem.value = query.value;
              newElem.classname = "formelem";
              document.forms['searchForm'].appendChild(newElem);
        
              submit_action = "/matlabcentral/fileexchange/";
              selected_index = 1;
              break
           case "answers":
              var newElem = document.createElement("input");
              newElem.type = "hidden";
              newElem.name = "term";
              newElem.value = query.value;
              newElem.classname = "formelem";
              document.forms['searchForm'].appendChild(newElem);
        
              submit_action = "/matlabcentral/answers/";
              selected_index = 2;
              break
           case "cssm":
              var newElem = document.createElement("input");
              newElem.type = "hidden";
              newElem.name = "search_string";
              newElem.value = query.value;
              newElem.classname = "formelem";
              document.forms['searchForm'].appendChild(newElem);
        
                  submit_action = "/matlabcentral/newsreader/search_results";
              selected_index = 3;
              break
           case "linkexchange":
              var newElem = document.createElement("input");
              newElem.type = "hidden";
              newElem.name = "term";
              newElem.value = query.value;
              newElem.classname = "formelem";
              document.forms['searchForm'].appendChild(newElem);
        
              submit_action = "/matlabcentral/linkexchange/";
              selected_index = 4;
              break
           case "blogs":
              var newElem = document.createElement("input");
              newElem.type = "hidden";
              newElem.name = "c[]";
              newElem.value = "blogs";
              document.forms['searchForm'].appendChild(newElem);
        
              selected_index = 5;
              break
           case "cody":
              var newElem = document.createElement("input");
              newElem.type = "hidden";
              newElem.name = "term";
              newElem.value = query.value;
              newElem.classname = "formelem";
              document.forms['searchForm'].appendChild(newElem);
        
              submit_action = "/matlabcentral/cody/";
              selected_index = 6;
              break
           case "contest":
              var newElem = document.createElement("input");
              newElem.type = "hidden";
              newElem.name = "c[]";
              newElem.value = "contest";
              document.forms['searchForm'].appendChild(newElem);
        
              selected_index = 7;
              break
           case "entire1":
              var newElem = document.createElement("input");
              newElem.type = "hidden";
              newElem.name = "c[]";
                  newElem.value = "entire_site";
              document.forms['searchForm'].appendChild(newElem);
              
              selected_index = 8;
              break
           default:
              var newElem = document.createElement("input");
              newElem.type = "hidden";
              newElem.name = "c[]";
              newElem.value = "entire_site";
              document.forms['searchForm'].appendChild(newElem);
           
              selected_index = 8;
              break
        }

        document.forms['searchForm'].elements['search_submit'].selectedIndex = selected_index;
        document.forms['searchForm'].elements['query'].value = query.value;
        document.forms['searchForm'].action = submit_action;
}

</SCRIPT>


  <form name="searchForm" method="GET" action="" style="margin:0px; margin-top:5px; font-size:90%" onSubmit="submitForm(query)">
          <label for="search">Search: </label>
        <select name="search_submit" style="font-size:9px ">
                 <option value = "matlabcentral">MATLAB Central</option>
                <option value = "fileexchange" selected>&nbsp;&nbsp;&nbsp;File Exchange</option>
                <option value = "answers">&nbsp;&nbsp;&nbsp;Answers</option>
            <option value = "cssm">&nbsp;&nbsp;&nbsp;Newsgroup</option>
                <option value = "linkexchange">&nbsp;&nbsp;&nbsp;Link Exchange</option>
                <option value = "blogs">&nbsp;&nbsp;&nbsp;Blogs</option>
                <option value = "cody">&nbsp;&nbsp;&nbsp;Cody</option>
                <option value = "contest">&nbsp;&nbsp;&nbsp;Contest</option>
          <option value = "entire1">MathWorks.com</option>
        </select>
<input type="text" name="query" size="10" class="formelem" value="">
<input type="submit" value="Go" class="formelem gobutton" >
</form>

      <ol id="access2">
  <li class="first">
    Jan Sieber
  </li>
  <li>
      <a href="https://uk.mathworks.com/mwaccount/profiles/edit?t=community" id="view_community_profile_link">My Community Profile</a>
  </li>
  <li>
    <a href="https://uk.mathworks.com/login/logout?uri=https%3A%2F%2Fuk.mathworks.com%2Fmatlabcentral%2Ffileexchange%2F13490-adaptive-robust-numerical-differentiation%2Fcontent%2FDERIVESTsuite%2Fderivest.m" id="logout_link">Log Out</a>
  </li>
</ol>

      </div>
        <div class="logo_container hidden-xs hidden-sm">
          <a href="/index.html?s_tid=gn_logo" class="svg_link pull-left">
            <img src="/images/responsive/global/pic-header-mathworks-logo.svg" class="mw_logo" alt="MathWorks">
          </a>
        </div>
        <div id="globalnav">
        <div class="navbar-header">
        <div class="navbar-collapse collapse hidden-xs hidden-sm">
          <ul class="nav navbar-nav" id="topnav">
            <li class="topnav_products "><a href="/products/?s_tid=gn_ps">Products</a></li>
            <li class="topnav_solutions "><a href="/solutions/?s_tid=gn_sol">Solutions</a></li>
            <li class="topnav_academia "><a href="/academia/?s_tid=gn_acad">Academia</a></li>
            <li class="topnav_support "><a href="/support/?s_tid=gn_supp">Support</a></li>
            <li class="topnav_community active"><a href="/matlabcentral/?s_tid=gn_mlc" class="dropdown-toggle"  role="button" aria-haspopup="true" aria-expanded="false">Community</a>
            </li>
            <li class="topnav_events "><a href="/company/events/?s_tid=gn_ev">Events</a></li>
          </ul>
        </div>
      </div>
        <!-- from includes/global_nav.html -->

      </div>
    </div>
  </div>
</div>

      <div id="middle">
  <div class="wrapper">
  	<div class="fileexchange-header">
  		<a href="/matlabcentral/fileexchange/?s_tid=gn_mlc_fx">File Exchange Home</a>
  	</div>

    <div id="mainbody" class="columns2">
  
  

  <div class="manifest">

      <div id="download_toolbox_button">
            <div class="btnCont ctaBtn ctaBlueBtn btnSmall">
              <div class="btn download"><a href="/matlabcentral/mlc-downloads/downloads/submissions/13490/versions/7/download/mltbx" class="link--download" data-file-format="mltbx" data-logintodownload="false" title="Download Now">Download Toolbox</a></div>
            </div>
          </div>
      <div id="download_zip_button">
            <div class="btnCont ctaBtn ctaBlueBtn btnSmall">
              <div class="btn download"><a href="/matlabcentral/mlc-downloads/downloads/submissions/13490/versions/7/download/zip" class="link--download" data-file-format="zip" data-logintodownload="false" title="Download Now">Download Zip</a></div>
            </div>
          </div>

      <p class="license">
      <a href="/matlabcentral/fileexchange/view_license?file_info_id=13490" popup="new_window height=500,width=640,scrollbars=yes">View License</a>
  </p>


  <style>
	#addons-callout {
		
	}
	#addons-callout .callout_content {
		border-top: 5px solid #009915;
		padding: 10px;
	}

	#addons-callout p.cta {
		margin-bottom: 0;
	}
	#addons-callout p.cta a {
		color: #005FCE;
	}
	#addons-callout p.cta a:hover {
		text-decoration: underline;
	}
</style>

<div class="callout_container solid margined_20" id="addons-callout">
	<div class="card_image card_medium" style="background-image:url(/images/nextgen/callouts/fileexchange-addons-spotlight.png); height: 110px; background-position: 42% 50%;"></div>
	<div class="callout_content card_color_tertiary "> <p>Download apps, toolboxes, and other File Exchange content using Add-On Explorer in MATLAB.</p>
		<p class="cta">
			<a href="/videos/add-on-explorer-106745.html">&raquo; Watch video</a>
		</p>
	</div>
</div>

  
  <h3 class="highlights_title">Highlights from <br/>
    <a href="https://uk.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation" class="manifest_title">Adaptive Robust Numerical Differentiation</a>
  </h3>
  <ul class='manifest'>
      <li class='manifest'><a href="https://uk.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation/content/DERIVESTsuite/demo/html/derivest_demo.html" class="example" title="Example">derivest_demo</a></li>
      <li class='manifest'><a href="https://uk.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation/content/DERIVESTsuite/demo/html/multivariable_calc_demo.html" class="example" title="Example">multivariable_calc_demo</a></li>
      <li class='manifest'><a href="https://uk.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation/content/DERIVESTsuite/derivest.m" class="function" title="Function">derivest(fun,x0,varargin)</a><span class="description">DERIVEST: estimate the n'th derivative of fun at x0, provide an error estimate</span></li>
      <li class='manifest'><a href="https://uk.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation/content/DERIVESTsuite/directionaldiff.m" class="function" title="Function">directionaldiff(fun,x0,vec)</a><span class="description">directionaldiff: estimate of the directional derivative of a function of n variables</span></li>
      <li class='manifest'><a href="https://uk.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation/content/DERIVESTsuite/gradest.m" class="function" title="Function">gradest(fun,x0)</a><span class="description">gradest: estimate of the gradient vector of an analytical function of n variables</span></li>
      <li class='manifest'><a href="https://uk.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation/content/DERIVESTsuite/hessdiag.m" class="function" title="Function">hessdiag(fun,x0)</a><span class="description">HESSDIAG: diagonal elements of the Hessian matrix (vector of second partials)</span></li>
      <li class='manifest'><a href="https://uk.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation/content/DERIVESTsuite/hessian.m" class="function" title="Function">hessian(fun,x0)</a><span class="description">hessian: estimate elements of the Hessian matrix (array of 2nd partials)</span></li>
      <li class='manifest'><a href="https://uk.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation/content/DERIVESTsuite/jacobianest.m" class="function" title="Function">jacobianest(fun,x0)</a><span class="description">gradest: estimate of the Jacobian matrix of a vector valued function of n variables</span></li>
    <li class="manifest_allfiles">
      <a href="https://uk.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation/all_files" id="view_all_files">View all files</a>
    </li>
  </ul>


</div>


  <table cellpadding="0" cellspacing="0" class="details file contents">
    <tr>
      <th class="maininfo">
        


<div id="details">
  <h1 itemprop="name">Adaptive Robust Numerical Differentiation</h1>
  <p id="author">
    by
    <span itemprop="author" itemscope itemtype="http://schema.org/Person">
      <span itemprop="name"><a href="/matlabcentral/profile/authors/869215-john-d-errico" class="author_inline results_author" data-cp-link-id="1">John D&#x27;Errico</a>

<div class="profile mlc_author_popover" data-cp-popup-id="1" vocab="http://schema.org/" typeof="Person">
  <h3 class="profile__name" property="name"><a href="/matlabcentral/profile/authors/869215-john-d-errico">John D&#x27;Errico <span>(view profile)</span></a></h3>

    <div class="profile__image">
      <a href="/matlabcentral/profile/authors/869215-john-d-errico"><img src="/responsive_image/100/0/0/0/0/cache/matlabcentral/profiles/869215.jpg" /></a>
    </div>
  <div class="profile__stats" typeof="CreativeWork">
    
  <ul>
    <li property="fileExchange">
      <span class="icon-fileexchange"></span> <a href="/matlabcentral/fileexchange/?term=authorid%3A679"><span property="fileExchange_count">56</span> files</a>
    </li>
    <li property="downloads">
      <span class="icon-download"></span> <span property="downloads_count">2832</span> downloads
    </li>
      <li property="rating">

        <div class="rating" title="4.80937">
          <div class="rate_scale">
            <div class="rated" property="aggregateRating" style="width: 96.18740000000001%">4.80937</div>
          </div>
        </div>

      </li>
  </ul>

  </div>
</div>
</span>
    </span>
  </p>
  <p>&nbsp;</p>
  <p>
    <span id="submissiondate"
>
      27 Dec 2006
    </span>
      <span id="date_updated">(Updated
        <span itemprop="datePublished" content="2014-12-03">03 Dec 2014</span>)
      </span>
  </p>

  <p id="summary">Numerical derivative of an analytically supplied function, also gradient, Jacobian &amp; Hessian</p>


  
</div>

        </div>
      </th>
    </tr>
    <tr>
      <td class="file">
        <table cellpadding="0" cellspacing="0" border="0" class="fileview section">
          <tr class="title">
            <th><span class="heading">derivest(fun,x0,varargin)</span></th>
          </tr>
          <tr>
            <td>
              <iframe id="content_iframe" style="width: 645px;min-height: 600px; border: none" frameborder="0" scrolling="no" sandbox="allow-popups allow-same-origin " src="/matlabcentral/mlc-downloads/downloads/submissions/13490/versions/7/previews/DERIVESTsuite/derivest.m/index.html"></iframe>
              <script>
                $('#content_iframe').load(function () { initIframe('content_iframe') })
                function getDocHeight(doc) {
                    doc = doc || document;
                    // stackoverflow.com/questions/1145850/
                    var body = doc.body, html = doc.documentElement;
                    var height = Math.max( body.scrollHeight, body.offsetHeight,
                        html.clientHeight, html.scrollHeight, html.offsetHeight );
                    return height;
                }
                function initIframe(id) {
                  var ifrm = document.getElementById(id);
                  var doc = ifrm.contentDocument ? ifrm.contentDocument :
                      ifrm.contentWindow.document;
                  disableAbsoluteLinks(doc);
                  setIframeHeight(ifrm, doc);
                }
                function disableAbsoluteLinks(doc) {
                  $(doc).find('a[href^="http"],a[href^="https"]').attr('target', '_TOP');
                }
                function setIframeHeight(ifrm, doc) {
                    ifrm.style.visibility = 'hidden';
                    ifrm.style.height = "10px"; // reset to minimal height ...
                    // IE opt. for bing/msn needs a bit added or scrollbar appears
                    ifrm.style.height = getDocHeight( doc ) + 25 + "px";
                    ifrm.style.visibility = 'visible';
                }
              </script>
            </td>
          </tr>
        </table>
      </td>
    </tr>

  </table>
  <script src="/matlabcentral/fileexchange/assets/file_infos/show.js" type="text/javascript"></script>

<p id="contactus"><a href="/company/feedback/">Contact us</a></p>

      	
      
</div>
<div class="clearboth">&nbsp;</div>
</div>
</div>
<!-- footer.html -->
<!-- START OF FOOTER -->
<style>
#mlc_footer_top { padding: 6px 5px 0; }
#mlc_footer_top .footer_column { width: 14%; float: left; }
#mlc_footer_top .footer_column:first-child { width: 30%; }
#mlc_footer_top .footer_column ul { text-align: left; list-style-type: none; }
#mlc_footer_top .footer_column ul li { display: block; padding: 0; border: 0; }
#mlc_footer_top .footer_column h3 { font-size: 16px; margin-bottom: 10px; }
#mlc_footer_bottom { padding: 6px 0 0 0; border-top: 1px #ccc solid; clear: both; }
</style>


<div id="mlc-footer">
  <script type="text/javascript">
function clickDynamic(obj, target_url, tracking_code) {
	var pos=target_url.indexOf("?");
	if (pos<=0) {
		var linkComponents = target_url + tracking_code;
		obj.href=linkComponents;
	}
}
</script>
  <div class="wrapper">
    <div id="mlc_footer_top">
      <div class="footer_column">
        <h3>MathWorks</h3>
        <p><em>Accelerating the pace of engineering and science</em></p>
        <p>MathWorks is the leading developer of mathematical computing software for engineers and scientists.</p>
        <p><a href="/discovery/">Discover...</a></p>
      </div>
      <div class="footer_column">
        <ul>
          <li><strong>Explore Products</strong></li>
          <li><a href="/products/matlab/?s_tid=mlc_fx_ff_p_matlab">MATLAB</a></li>
          <li><a href="/products/simulink/?s_tid=mlc_fx_ff_p_simulink">Simulink</a></li>
          <li><a href="/academia/student_version/?s_tid=mlc_fx_ff_p_student">Student Software</a></li>
          <li><a href="/hardware-support/home.html?s_tid=mlc_fx_ff_p_hwsupport">Hardware Support</a></li>
          <li><a href="/matlabcentral/fileexchange/?s_tid=mlc_fx_ff_p_fx">File Exchange</a></li>
        </ul>
      </div>
      <div class="footer_column">
        <ul>
          <li><strong>Try or Buy</strong></li>
          <li><a href="/downloads/web_downloads/?s_iid=mlc_fx_ff_t_downloads">Downloads</a></li>
          <li><a href="/programs/trials/trial_request.html?s_iid=mlc_fx_ff_p_trial">Trial Software</a></li>
          <li><a href="/company/aboutus/contact_us/contact_sales.html?s_iid=mlc_fx_ff_t_sales">Contact Sales</a></li>
          <li><a href="/pricing-licensing/?s_iid=mlc_fx_ff_t_pricing">Pricing and Licensing</a></li>
        </ul>
      </div>
      <div class="footer_column">
        <ul>
          <li><strong>Learn to Use</strong></li>
          <li><a href="/help/?s_tid=mlc_fx_ff_l_doc">Documentation</a></li>
          <li><a href="/support/learn-with-matlab-tutorials.html?s_tid=mlc_fx_ff_l_tutorials">Tutorials</a></li>
          <li><a href="/help/matlab/examples.html?s_tid=mlc_fx_ff_l_examples">Examples</a></li>
          <li><a href="/videos/?s_tid=mlc_fx_ff_l_videos">Videos and Webinars</a></li>
          <li><a href="/services/training/?s_tid=mlc_fx_ff_l_training">Training</a></li>
        </ul>
      </div>
      <div class="footer_column">
        <ul>
          <li><strong>Get Support</strong></li>
          <li><a href="/support/install-matlab.html?s_tid=mlc_fx_ff_s_install">Installation Help</a></li>
          <li><a href="/matlabcentral/answers/?s_tid=mlc_fx_ff_s_answers">Answers</a></li>
          <li><a href="/services/consulting/?s_tid=mlc_fx_ff_s_consulting">Consulting</a></li>
          <li><a href="/licensecenter/?s_tid=mlc_fx_ff_s_license">License Center</a></li>
        </ul>
      </div>
      <div class="footer_column">
        <ul>
          <li><strong>About MathWorks</strong></li>
          <li><a href="/company/jobs/opportunities/?s_tid=mlc_fx_ff_a_careers">Careers</a></li>
          <li><a href="/company/?s_tid=mlc_fx_ff_a_company">Company Overview</a></li>
          <li><a href="/company/newsroom/?s_tid=mlc_fx_ff_a_newsroom">Newsroom</a></li>
          <li><a href="/company/aboutus/soc_mission/?s_tid=mlc_fx_ff_a_socialmission">Social Mission</a></li>
        </ul>
      </div>
    </div>
    <div id="mlc_footer_bottom">
      <ul id="matlabcentral">
        <li class="copyright first">&copy; 1994-2017 The MathWorks, Inc.</li>
        <li class="first"><a href="/company/aboutus/policies_statements/patents.html?s_tid=gf_pat" title="patents" rel="nofollow">Patents</a></li>
        <li><a href="/company/aboutus/policies_statements/trademarks.html?s_tid=gf_trd" title="trademarks" rel="nofollow">Trademarks</a></li>
        <li><a href="/company/aboutus/policies_statements/?s_tid=gf_priv" title="privacy policy" rel="nofollow">Privacy Policy</a></li>
        <li><a href="/company/aboutus/policies_statements/piracy.html?s_tid=gf_pir" title="preventing piracy" rel="nofollow">Preventing Piracy</a></li>
        <li class="last"><a href="/matlabcentral/termsofuse.html?s_tid=gf_com_trm" title="Terms of Use" rel="nofollow">Terms of Use</a></li>
        <li class="icon"><a href="/company/rss/" title="RSS" class="rssfeed" rel="nofollow"><span class="text">RSS</span></a></li>
        <li class="icon"><a href="/programs/bounce_hub_generic.html?s_tid=mlc_lkd&url=http://www.linkedin.com/company/the-mathworks_2" title="LinkedIn" class="linkedin" rel="nofollow" target="_blank"></a></li>
        <li class="icon"><a href="/programs/bounce_hub_generic.html?s_tid=mlc_fbk&url=https://plus.google.com/117177960465154322866?prsrc=3" title="Google+" class="google" rel="nofollow" target="_blank"><span class="text">Google+</span></a></li>
        <li class="icon"><a href="/programs/bounce_hub_generic.html?s_tid=mlc_fbk&url=http://www.facebook.com/MATLAB" title="Facebook" class="facebook" rel="nofollow" target="_blank"><span class="text">Facebook</span></a></li>
        		<li class="last icon"><a href="/programs/bounce_hub_generic.html?s_tid=mlc_twt&url=http://www.twitter.com/MATLAB" title="Twitter" class="twitter" rel="nofollow" target="_blank"><span class="text">Twitter</span></a></li>

      </ul>

    </div>
  </div>
</div>
<!-- END OF FOOTER -->


      
      
<!-- BEGIN Adobe DTM -->
<script type="text/javascript">
try {
_satellite.pageBottom();
} catch (e) {
//something went wrong
}
</script>
<!-- END Adobe DTM -->
  

     

    </body>
</html>
