#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
##############
Text Templates
##############

*Created on Fri, Sep 1, 2017 by A. Pahl*

Text templates, mainly for the viewers.
"""

from string import Template


class PreTemplate(Template):
    delimiter = "§"


TABLE_OPTIONS = {"cellspacing": "1", "cellpadding": "1", "border": "1",
                 "align": "", "height": "60px", "summary": "Table", }  # "width": "800px",

# PAGE_OPTIONS = {"icon": "icons/chart_bar.png", "css": ["css/style.css", "css/collapsible_list.css"],
#                 "scripts": ["lib/jquery.js", "lib/highcharts.js", "script/folding.js"]}
PAGE_OPTIONS = {"icon": "icons/benzene.png"}
JSME = "lib/jsme/jsme.nocache.js"

HTML_FILE_NAME = "mol_table.html"

STYLESHEETS_LOCAL = """<link rel="stylesheet" href="lib/css/jquery.dataTables.min.css">
    <link rel="stylesheet" href="lib/css/bootstrap1.min.css">
    <link rel="stylesheet" href="lib/css/bootstrap2.min.css">
    <link rel="stylesheet" href="lib/css/fixedHeader.dataTables.min.css">
    """
# <link rel="stylesheet" href="lib/css/bootstrap-theme.min.css">
JS_LIBS_LOCAL = """    <!-- jQuery -->
    <script src="lib/jquery.min.js"></script>
    <!-- Popper -->
    <script src="lib/popper.js"></script>
    <!-- Bootstrap -->
    <script src="lib/bootstrap.min.js"></script>
    <!-- Data Table -->
    <script src="lib/jquery.dataTables.min.js"></script>
    <script src="lib/dataTables.fixedHeader.min.js"></script>"""

STYLESHEETS_NET = """<link rel="stylesheet" href="https://cdn.datatables.net/1.10.16/css/jquery.dataTables.min.css">
    <link rel="stylesheet" href="https://bootswatch.com/cosmo/bootstrap.min.css">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-beta/css/bootstrap.min.css" integrity="sha384-/Y6pD6FV/Vv2HJnA6t+vslU6fwYXjCFtcEpHbNJ0lyAFsXTsjBbfaDjzALeQsN6M" crossorigin="anonymous">
    <link rel="stylesheet" href="https://cdn.datatables.net/fixedheader/3.1.3/css/fixedHeader.dataTables.min.css" crossorigin="anonymous">
"""

JS_LIBS_NET = """<!-- jQuery -->
    <script
    src="https://code.jquery.com/jquery-3.2.1.min.js"
    integrity="sha256-hwg4gsxgFZhOsEEamdOYGBf13FyQuiTwlAQgxVSNgt4="
    crossorigin="anonymous"></script>
    <!-- Popper -->
    <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.5/esm/popper.min.js"></script>
    <!-- Bootstrap -->
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-beta/js/bootstrap.min.js" integrity="sha384-h0AbiXch4ZDo7tp9hKZ4TsHbi047NrKGLO3SEJAg45jXxnGIfYzk4Si90RDIqNm1" crossorigin="anonymous"></script>
    <!-- Data Table -->
    <script src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/fixedheader/3.1.3/js/dataTables.fixedHeader.min.js"></script>"""


SELECTION_BTN = """<button id="show_sel">Selected Ids</button>"""

SELECTION_JS = """
    <script>
    $(document).ready(function() {
      var table = $('#maintable').DataTable();

      $('#maintable tbody').on( 'click', 'tr', function () {
          $(this).toggleClass('selected');
      } );

      $('#show_sel').click( function () {
          var id_arr = [];
          var rows = table.rows('.selected').data()
          $.each($(rows), function(idx, row){
              id_arr.push(row[1]);
          });
          alert( id_arr );
      } );
    } );
    </script>"""

PANDAS_TABLE = """
<!DOCTYPE html>
<html>
  <head>
    <meta charset="utf-8">
    <title>$title</title>
    §style_sheets
    <style>
      .custom {
        font-size: 14px;
      }
      .bottomcustom, .topcustom {
        font-size: 10px;
      }
    </style>
  </head>
  <body>
    <div class=custom>$intro</div>
    <div class="container">
      <div class="page-header">
          <h3>$title</h3>
      </div>
      $selection_btn
      <table id="maintable" class="display custom" cellspacing="0" width="100%">
        $table
      </table>
    </div>

    §js_libs
    <script>
        $$(document).ready(function() {
            $$('#maintable').DataTable({
                fixedHeader: true,
                "pageLength": 25,
                "dom": '<"topcustom"lfr>t<"bottomcustom"ip>'
            });
        });
    </script>
    $selection_js
  </body>
</html>"""
t = PreTemplate(PANDAS_TABLE)
d = {"style_sheets": STYLESHEETS_LOCAL, "js_libs": JS_LIBS_LOCAL}
PANDAS_TABLE_LOCAL = t.substitute(d)
d = {"style_sheets": STYLESHEETS_NET, "js_libs": JS_LIBS_NET}
PANDAS_TABLE_NET = t.substitute(d)

TBL_JAVASCRIPT = '''<script type="text/javascript">

function toggleCpd(cpdIdent)
{{
  listPos = document.id_list{ts}.data.value.indexOf(cpdIdent);
  cpdIdentCell = document.getElementById(cpdIdent+"_{ts}");
  if (listPos == -1)
  {{
    if (document.id_list{ts}.remark.checked == true)
    {{
      rem = "\\t" + prompt("Remark (Enter for none):", "");
    }}
    else
    {{
      rem = "";
    }}
    document.id_list{ts}.data.value = document.id_list{ts}.data.value + cpdIdent + rem + "\\n";
    cpdIdentCell.style.backgroundColor = "yellow";
  }}
  else
  {{
    removeStr = cpdIdent;
    tempStr2 = document.id_list{ts}.data.value;
    if (listPos > 0) {{
      tempStr1 = tempStr2.substring(0, listPos);
      tempStr2 = tempStr2.substring(listPos, tempStr2.length);
    }} else {{
      tempStr1 = "";
    }}
    listPos = tempStr2.indexOf("\\n");
    if (listPos < tempStr2.length - 1) {{
      tempStr1 = tempStr1 + tempStr2.substring(listPos+1, tempStr2.length)
    }}
    document.id_list{ts}.data.value = tempStr1;
    cpdIdentCell.style.backgroundColor = "{bgcolor}";
  }}
  show_number_selected();
}}

function show_number_selected() {{
  // display the number of selected compounds:
  var count = (document.id_list{ts}.data.value.match(/\\n/g) || []).length;
  document.getElementById("selection_title{ts}").innerHTML = "Selection (" + count + "):";
}}

function highlight_cpds() {{
  // highlights compounds that were pasted into the selection list
  // and keeps those that could be found
  var lines = document.id_list{ts}.data.value.split("\\n");
  var found = "";
  for (var idx = 0; idx < lines.length; idx++) {{
    var cpd = lines[idx];
    var cpdIdentCell = document.getElementById(cpd+"_{ts}");
    if (cpdIdentCell != null) {{
      cpdIdentCell.style.backgroundColor = "yellow";
      found = found + cpd +  "\\n";
    }}
  }}
  // set the value of the selection list to the found compound Ids
  document.id_list{ts}.data.value = found;
  show_number_selected();
}}

function myShowSelection() {{
  document.location.hash = "#SelectionList";
}}
</script>
'''

ID_LIST = """<br><b><a name="SelectionList" id="selection_title{ts}">Selection (0):</a></b>
<form name="id_list{ts}">
<input type="checkbox" name="remark" value="prompt" > Prompt for Remarks<br>
<textarea name="data" cols="70" rows="10"></textarea>
<input type="button" name="highlight" value="highlight compounds" onclick="highlight_cpds()"
title="Paste a list of Compound_Ids here and press this button. The compounds will be highlighted in the report above. Compounds which were not found are removed from the list.">
</form>
"""

JSME_FORM = '''<script type="text/javascript" src="{jsme_loc}/jsme/jsme.nocache.js"></script>
<script type="text/javascript">

function jsmeOnLoad() {{
    //arguments: HTML id, width, height (must be string not number!)
    jsmeApplet{ts} = new JSApplet.JSME("appletContainer{ts}", "380px", "340px", {{
                     //optional parameters
                     "options" : "query,hydrogens"
    }});
}}

function onSubmit() {{
    var drawing = jsmeApplet{ts}.molFile();
    // document.getElementById('jsme_smiles{ts}').value = drawing;
    var command = '{var_name} = Chem.MolFromMolBlock("""' + drawing + '""")';
    console.log("Executing Command: " + command);

    var kernel = IPython.notebook.kernel;
    kernel.execute(command);
}}
</script>

<table align="left" style="border: none;">
<tr style="border: none;">
<td id="appletContainer{ts}" style="border: none;"></td>
<td style="vertical-align: bottom; border: none;">
<button onclick="onSubmit()">done !</button>
</td>
</tr>
</table>
'''


def tag(name, content, options=None, lf_open=False, lf_close=False):
    """creates a HTML stub with closed tags of type <name> around <content>
    with additional <options::dict> in the opening tag
    when lf_(open|close)==True, the respective tag will be appended with a line feed.
    returns: html stub as list"""

    if lf_open:
        lf_open_str = "\n"
    else:
        lf_open_str = ""
    if lf_close:
        lf_close_str = "\n"
    else:
        lf_close_str = ""

    option_str = ""
    if options:
        option_list = [" "]
        for option in options:
            option_list.extend([option, '="', str(options[option]), '" '])

        option_str = "".join(option_list)

    stub = ["<{}{}>{}".format(name, option_str, lf_open_str)]
    if type(content) == list:
        stub.extend(content)
    else:
        stub.append(content)

    stub.append("</{}>{}".format(name, lf_close_str))

    return stub


def page(content, title="Results", header=None, summary=None, options=PAGE_OPTIONS):
    """create a full HTML page from a list of stubs below
    options dict:
    css:     list of CSS style file paths to include.
    scripts: list of javascript library file paths to include.
    icon: path to icon image
    returns HTML page as STRING !!!"""

    # override the title if there is a title in <options>
    if "title" in options and len(options["title"]) > 2:
        title = options["title"]

    if "icon" in options and len(options["icon"]) > 2:
        icon_str = '<link rel="shortcut icon" href="{}" />'.format(options["icon"])
    else:
        icon_str = ""

    if "css" in options and options["css"]:
        css = options["css"]
        if not isinstance(css, list):
            css = [css]

        css_str = "".join(['  <link rel="stylesheet" type="text/css" href="{}">\n'.format(file_name) for file_name in css])

    else:
        # minimal inline CSS
        css_str = """<style>
  body{
  background-color: #FFFFFF;
  font-family: freesans, arial, verdana, sans-serif;
}
th {
  border-collapse: collapse;
  border-width: thin;
  border-style: solid;
  border-color: black;
  text-align: left;
  font-weight: bold;
}
td {
  border-collapse:collapse;
  border-width:thin;
  border-style:solid;
  border-color:black;
  padding: 5px;
}
table {
  border-collapse:collapse;
  border-width:thin;
  border-style:solid;
  //border-color:black;
  border: none;
  background-color: #FFFFFF;
  text-align: left;
}
</style>"""

    if "scripts" in options and options["scripts"]:
        scripts = options["scripts"]
        if type(scripts) != list:
            scripts = [scripts]

        js_str = "".join(['  <script src="{}"></script>\n'.format(file_name) for file_name in scripts])

    else:
        js_str = ""

    if header:
        if not isinstance(header, list):
            header = [header]

        header_str = "".join(h2(header))

    else:
        header_str = ""

    if summary:
        if not isinstance(summary, list):
            summary = [summary]

        summary_str = "".join(p(summary))

    else:
        summary_str = ""



    if isinstance(content, list):
        content_str = "".join(content)
    else:
        content_str = content


    html_page = """<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>{title}</title>
  {icon_str}
{css_str}
{js_str}
</head>
<body>
{header_str}
{summary_str}
{content_str}
</body>
</html>
""".format(title=title, icon_str=icon_str, css_str=css_str, js_str=js_str,
           header_str=header_str, summary_str=summary_str, content_str=content_str)

    return html_page


def write(text, fn=HTML_FILE_NAME):
    with open(fn, "w") as f:
        f.write(text)


def script(content):
    return tag("script", content, lf_open=True, lf_close=True)


def img(src, options=None):
    """takes a src, returns an img tag"""

    option_str = ""
    if options:
        option_list = [" "]
        for option in options:
            option_list.extend([option, '="', str(options[option]), '" '])

        option_str = "".join(option_list)

    stub = ['<img {}src="data:image/png;base64,{}" alt="icon" />'.format(option_str, src)]

    return stub


def table(content, options=TABLE_OPTIONS):
    tbody = tag("tbody", content, lf_open=True, lf_close=True)
    return tag("table", tbody, options, lf_open=True, lf_close=True)


def tr(content, options=None):
    return tag("tr", content, options, lf_close=True)


def td(content, options=None):
    return tag("td", content, options, lf_close=False)


def p(content):
    return tag("p", content, lf_open=True, lf_close=True)


def h1(content):
    return tag("h1", content, lf_open=True, lf_close=True)


def h2(content):
    return tag("h2", content, lf_open=True, lf_close=True)


def div(content, options=None):
    return tag("div", content, options, lf_close=False)


def ul(content):
    return tag("ul", content, lf_open=True, lf_close=True)


def li(content):
    return tag("li", content, lf_open=False, lf_close=True)


def li_lf(content):  # list item with opening line feed
    return tag("li", content, lf_open=True, lf_close=True)


def b(content, options=None):
    return tag("b", content, options, lf_close=False)


def a(content, options):
    """the anchor tag requires an "href" in options,
    therefore the "options" parameter is not optional in this case (klingt komisch, ist aber so)"""
    return tag("a", content, options, lf_close=False)
