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


STYLESHEETS_LOCAL = """<link rel="stylesheet" href="lib/css/jquery.dataTables.min.css">
    <link rel="stylesheet" href="lib/css/bootstrap1.min.css">
    <link rel="stylesheet" href="lib/css/bootstrap2.min.css">
    """
# <link rel="stylesheet" href="lib/css/bootstrap-theme.min.css">
JS_LIBS_LOCAL = """    <!-- jQuery -->
    <script src="lib/jquery.min.js"></script>
    <!-- Popper -->
    <script src="lib/popper.js"></script>
    <!-- Bootstrap -->
    <script src="lib/bootstrap.min.js"></script>
    <!-- Data Table -->
    <script src="lib/jquery.dataTables.min.js"></script>"""

STYLESHEETS_NET = """<link rel="stylesheet" href="https://cdn.datatables.net/1.10.16/css/jquery.dataTables.min.css">
    <link rel="stylesheet" href="https://bootswatch.com/cosmo/bootstrap.min.css">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-beta/css/bootstrap.min.css" integrity="sha384-/Y6pD6FV/Vv2HJnA6t+vslU6fwYXjCFtcEpHbNJ0lyAFsXTsjBbfaDjzALeQsN6M" crossorigin="anonymous">
"""

JS_LIBS_NET = """<!-- jQuery -->
    <script
    src="https://code.jquery.com/jquery-3.2.1.min.js"
    integrity="sha256-hwg4gsxgFZhOsEEamdOYGBf13FyQuiTwlAQgxVSNgt4="
    crossorigin="anonymous"></script>
    <!-- Popper -->
    <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.5/esm/popper.min.js">
    <!-- Bootstrap -->
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0-beta/js/bootstrap.min.js" integrity="sha384-h0AbiXch4ZDo7tp9hKZ4TsHbi047NrKGLO3SEJAg45jXxnGIfYzk4Si90RDIqNm1" crossorigin="anonymous"></script>
    <!-- Data Table -->
    <script src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"></script>"""


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
              id_arr.push(row[2]);
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
        font-size: 12px;
      }
      .bottomcustom, .topcustom {
        font-size: 10px;
      }
    </style>
  </head>
  <body>

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
