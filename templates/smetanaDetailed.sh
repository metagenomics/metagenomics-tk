smetana !{xmls} !{params.steps.metabolomics.smetana.additionalParams.detailed} -d -o !{sample}
sed -e "1s/^/SAMPLE\t/" -e "2,$ s/^/f\!{sample}/" !{sample}_detailed.tsv

