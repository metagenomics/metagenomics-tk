# Rows have a varying number of columns (some end with a trailing tab). Cut them to the 11 columns
# SAMPLE, BIN_ID and the 9 Metabuli columns, then remove SAMPLE and BIN_ID by name.
# Example: one classifications file per bin, e.g. test1_binned_bin.1.gtdb_classifications.tsv
for cl in !{classifications}; do
  # Example row with trailing tab (12 fields): "test1 bin.1 1 contig_1 618910136 3054 0.16 0 no rank d_Bacteria;...;-_GCF_900129225.1 618910136:51 <empty>"
  # -> NF>11{NF=11} cuts it to 11 fields, the empty trailing field is gone. Rows with 11 fields stay as they are.
  # -> csvtk then removes the columns SAMPLE and BIN_ID by name, so "test1 bin.1 1 contig_1 ..." becomes "1 contig_1 ..."
  awk -F'\t' 'BEGIN{OFS="\t"} NF>11{NF=11}1' ${cl} | csvtk -C "&" -t -T cut -f -SAMPLE,-BIN_ID >${cl}.mod
done

# Create contig classification to gene mapping
# Example: contig "contig_1" has the genes with locus tags "test1_bin.1_94129_afe209_BJFFBIFK_03107" and "..._03108"
# -> the join (contig_1 == name) gives one row per gene with the taxonomy of its contig.
# The mapping files of all bins are concatenated first (csvtk concat), as are all *.mod classification files.
csvtk -C "&" -t -T join -f "CONTIG;name" <(csvtk -C "&" -t -T concat !{contig2GeneMapping} |
  # Example: from SAMPLE BIN_ID CONTIG locus_tag and many more columns only these four are kept: "test1 bin.1 contig_1 test1_bin.1_94129_afe209_BJFFBIFK_03107"
  csvtk -C "&" -t -T cut -f SAMPLE,BIN_ID,CONTIG,locus_tag) <(csvtk -C "&" -t -T concat *.mod) |
  # Select needed columns for EMGB compatibple output
  # Example: "test1_bin.1_94129_afe209_BJFFBIFK_03107  618910136  no rank  d_Bacteria;p_Bacillota;c_Negativicutes;o_Selenomonadales;f_Selenomonadaceae;g_Schwartzia;s_Schwartzia succinivorans;-_GCF_900129225.1"
  csvtk -C "&" -t -T cut -f locus_tag,taxID,rank,lineage |
  # Header row must be removed for EMGB
  # Example: the line "locus_tag taxID rank lineage" is dropped, the data rows start directly.
  tail -n +2 |
  # Contigs without a hit (taxID 0) are set to "no rank" and "unclassified".
  # For hits, only the standard ranks (d,p,c,o,f,g,s) of the lineage are kept, which drops the genome accession
  # and other rankless nodes. Name is the last kept element, a rank of "no rank" is replaced by the rank of that element
  # and missing ranks are filled up to 8 levels with uc_<last known name> as in the MMseqs2 taxonomy output.
  # Example 1 (genome hit):  gene_1  618910136  no rank  d_Bacteria;...;s_Schwartzia succinivorans;-_GCF_900129225.1
  #                       -> gene_1  618910136  species  Schwartzia succinivorans  Bacteria;...;Schwartzia succinivorans;uc_Schwartzia succinivorans
  # Example 2 (order hit):   gene_2  80  order  d_Bacteria;p_Bacteroidota;c_Bacteroidia;o_Bacteroidales
  #                       -> gene_2  80  order  Bacteroidales  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;uc_Bacteroidales;uc_Bacteroidales;uc_Bacteroidales;uc_Bacteroidales
  # Example 3 (no hit):      gene_3  0  -  -
  #                       -> gene_3  0  no rank  unclassified  (empty lineage)
  awk -F'\t' 'BEGIN{OFS="\t"; LEVELS=8; split("superkingdom phylum class order family genus species", RANKS, " ")}
    {
      # Example: $4 = "d_Bacteria;p_Bacillota;s_Schwartzia succinivorans;-_GCF_900129225.1" -> n = 4 elements
      n = split($4, elements, ";"); kept = 0; lastRank = 0
      for (i = 1; i <= n; i++) {
        # Example: element "s_Schwartzia succinivorans" -> first letter "s" is at position 7 of "dpcofgs", so pos = 7
        # Example: element "-_GCF_900129225.1" -> "-" is not found, pos = 0
        pos = index("dpcofgs", substr(elements[i], 1, 1))
        # Example: "s_Schwartzia succinivorans" is kept as "Schwartzia succinivorans" (prefix removed), lastRank = 7.
        # "-_GCF_900129225.1" and "sp_Craniata" are skipped (no standard rank, or second character is not an underscore).
        if (substr(elements[i], 2, 1) == "_" && pos > 0) { kept++; lineage[kept] = substr(elements[i], 3); lastRank = pos }
      }
      # Example: taxID 0 (no hit) or a lineage without any standard rank like "-_cellular organisms" -> unclassified row
      if ($2 == 0 || kept == 0) { print $1, 0, "no rank", "unclassified", ""; next }
      # Example: kept = 7 -> name = "Schwartzia succinivorans" (the last kept element)
      name = lineage[kept]
      # Example: kept = 4 (order hit) -> levels 5 to 8 become "uc_Bacteroidales"
      for (i = kept + 1; i <= LEVELS; i++) lineage[i] = "uc_" name
      # Example: result starts with the first element, e.g. "Bacteria"
      result = lineage[1]
      # Example: the other elements are appended with ";" -> "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;uc_Bacteroidales;..."
      for (i = 2; i <= LEVELS; i++) result = result ";" lineage[i]
      # Example: rank "no rank" with lastRank = 7 -> "species". A given rank like "order" is kept as it is.
      # Output columns: locus_tag, taxID, rank, name, lineage
      print $1, $2, ($3 == "no rank" ? RANKS[lastRank] : $3), name, result
    }' >"!{sample}_!{binType}.!{dbType}.taxonomy.tsv"
