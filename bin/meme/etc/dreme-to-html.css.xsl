<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:template name="dreme-to-html.css">
    <![CDATA[
    /* START INCLUDED FILE "dreme-to-html.css" */
      table.dreme_motifs tr th, table.dreme_motifs tr td {
        padding: 0px 10px;
      }

      div.popup_wrapper {
        position:fixed; 
        z-index:2;
        width:100%; 
        height:0; 
        top:50%; 
        left:0;
      }

      div.popup {
        width: 400px; 
        z-index:2;
        margin-left: auto;
        margin-right: auto;
        padding: 5px;
        background: #FFF;
        border-style: double;
        border-width: 5px;
        border-color: #00666a;
        position:relative; 
      }

      div.grey_background {
        position:fixed; 
        z-index: 1;
        background-color: #000;
        -ms-filter: "progid:DXImageTransform.Microsoft.Alpha(Opacity=50)";
        filter: alpha(opacity=50);
        opacity: 0.5;
        left: 0;
        top: 0;
        width: 100%;
        height: 100%;

      }
      td.symaction {
        text-align: center;
      }
      *.symaction {
        font-size: 20px;
      }

      div.close {
        cursor: pointer;
        border: 1px solid black; 
        width:15px; 
        height:15px; 
        line-height:15px; /* this causes vertical centering */
        text-align:center; 
        background-color:#FFF; 
        color:#000; 
        font-size:15px;
        font-family:monospace;
      }

      div.close:hover {
        color:#FFF;
        background-color:#000; 
      }

      div.navnum {
        width:100%; 
        height:20px; 
        line-height:20px; 
        text-align:center; 
        font-size:medium;
      }

      a.navarrow {
        font-size: 30px;
        text-decoration:none;
      }

      a.inactive {
        color:#CCC;
      }

      div.actionbutton { 
        cursor: pointer;
        font-size: 18px;
        line-height:20px; 
        padding: 5px; 
        margin: 10px 0; 
        border: 1px solid black;
      }

      div.actionbutton:hover {
        color:#FFF;
        background-color:#000;
      }

      div.pop_content {
        position:absolute;
        z-index:1;
        width:300px;
        padding: 5px;
        background: #E4ECEC;
        font-size: 12px;
        font-family: Arial;
        border-style: double;
        border-width: 3px;
        border-color: #AA2244;
        display:none;
      }
      span.sort_dir {
        text-decoration: none;
      }

      div.section_title {
        font-weight: bold;
        cursor: pointer;
      }

      div.section_title.inactive {
        color: #000;
      }

      div.section_title.inactive:hover {
        color: #000;
        text-decoration:underline;
      }

      div.section_title label {
        cursor: pointer;
      }

      span.ellipsis {
        display: inline-block;
        border: 1px solid black;
        padding: 0 2px;
        margin: 0 2px;
      }

      div.section_title.inactive:hover span.ellipsis {
        color: #FFF;
        background-color: #000;
      }

      div.section_title.inactive span.toggle {
        color: #000;
      }

      div.section_data {
        margin-left: 20px;
      }
      tr.rule td, tr.rule th {
        border-bottom: 1px solid #CCC;
      }

      h1.compact, h2.compact, h3.compact, h4.compact, h5.compact, h6.compact {
        margin:0; 
        padding:0;
      }

      ul.programs {
        margin-top: 0;
        padding-top: 0;
        margin-bottom: 0;
        padding-bottom: 0;
        margin-left: 0;
        padding-left: 0;
        list-style: none;
        border-bottom: 1px solid black;
      }

      ul.programs li {
        border: 1px solid black;
        border-bottom-width: 0;
        background-color: #EFE;
        cursor: default;
      }

      ul.programs li.active {
        background-color: #CFC;
      }

      ul.programs li.selected {
        background-color: #262;
        color: #FFF;
      }

      div.programs_scroll {
        width: 100%; 
        height: 90px; 
        overflow-y: auto; 
        overflow-x: hidden;
        margin: 0 auto; 
      }
    /* END INCLUDED FILE "dreme-to-html.css" */
    ]]>
  </xsl:template>
</xsl:stylesheet>

