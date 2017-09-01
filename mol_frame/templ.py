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
    <link rel="stylesheet" href="lib/css/bootstrap-theme.min.css">"""

JS_LIBS_LOCAL = """    <!-- jQuery -->
    <script src="lib/jquery.min.js"></script>
    <!-- Bootstrap -->
    <script src="lib/bootstrap.min.js"></script>
    <!-- Data Table -->
    <script src="lib/jquery.dataTables.min.js"></script>"""

STYLESHEETS_NET = """<link rel="stylesheet" href="https://cdn.datatables.net/1.10.13/css/jquery.dataTables.min.css">
    <link rel="stylesheet" href="https://bootswatch.com/cosmo/bootstrap.min.css">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap-theme.min.css" integrity="sha384-rHyoN1iRsVXV4nD0JutlnGaslCJuC7uwjduW9SVrLvRYooPp2bWYgmgJQIXwl/Sp" crossorigin="anonymous">"""

JS_LIBS_NET = """<!-- jQuery -->
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js"></script>
    <!-- Bootstrap -->
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>
    <!-- Data Table -->
    <script src="https://cdn.datatables.net/1.10.13/js/jquery.dataTables.min.js"></script>"""

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
  </body>
</html>"""
t = PreTemplate(PANDAS_TABLE)
d = {"style_sheets": STYLESHEETS_LOCAL, "js_libs": JS_LIBS_LOCAL}
PANDAS_TABLE_LOCAL = t.substitute(d)
d = {"style_sheets": STYLESHEETS_NET, "js_libs": JS_LIBS_NET}
PANDAS_TABLE_NET = t.substitute(d)
