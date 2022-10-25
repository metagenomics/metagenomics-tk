smetana !{xmls} !{params.steps.metabolomics.smetana.additionalParams.global} -o !{sample}
sed -i -e '1s/^/SAMPLE\t/' -e "2s/^/!{sample}\t/" !{sample}_global.tsv

