<?xml version="1.0"?>
<xsl:stylesheet version="1.0"
 xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns:cis="http://zlab.bu.edu/schema/cisml"
 xmlns:mem="http://noble.gs.washington.edu/meme"
>
  <xsl:output method="text" />

  <!-- A stylesheet to create a tab delimited file from MCAST CisML.-->

  <xsl:template match="/*">
    <xsl:text>#pattern name&#9;sequence name&#9;start&#9;stop&#9;score&#9;p-value&#9;q-value&#9;matched sequence&#10;</xsl:text>
    <xsl:apply-templates select="cis:multi-pattern-scan/mem:match">
      <xsl:sort order="descending" data-type="number" select="../@score"/>
    </xsl:apply-templates>
  </xsl:template>

  <xsl:template match="mem:match">
   <!-- pattern name -->
   <xsl:value-of select="@cluster-id"/><xsl:text>&#9;</xsl:text>
   <!-- sequence name -->
   <xsl:value-of select="@seq-name"/>
   <xsl:text>&#9;</xsl:text>
   <!-- start -->
   <xsl:value-of select="@start"/>
   <xsl:text>&#9;</xsl:text>
   <!-- end -->
   <xsl:value-of select="@stop"/>
   <xsl:text>&#9;</xsl:text>
   <!-- score -->
   <xsl:value-of select="../@score"/>
   <xsl:text>&#9;</xsl:text>
   <!-- p-value -->
   <xsl:value-of select="../@pvalue"/>
   <xsl:text>&#9;</xsl:text>
   <!-- e-value -->
   <xsl:value-of select="@evalue"/>
   <xsl:text>&#9;</xsl:text>
   <!-- q-value -->
   <xsl:value-of select="@qvalue"/>
   <xsl:text>&#9;</xsl:text>
   <!-- Matched sequence -->
   <xsl:value-of select="substring(text(), 2)"/>
  </xsl:template>

</xsl:stylesheet>
