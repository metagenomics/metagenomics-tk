memote run !{params.steps.metabolomics.memote.additionalParams.run} --filename !{sample}_!{id}_report.json.gz !{model}
memote report snapshot !{params.steps.metabolomics.memote.additionalParams.report} --filename !{sample}_!{id}_report.html !{model}
values=$(zcat !{sample}_!{id}_report.json.gz  | jq -r ' [ .tests.test_stoichiometric_consistency.duration, 
               .tests.test_reaction_mass_balance.metric, 
               .tests.test_reaction_charge_balance.metric, 
               .tests.test_find_disconnected.metric, 
               .tests.test_find_reactions_unbounded_flux_default_condition.metric ] | @tsv ')
title=$(zcat !{sample}_!{id}_report.json.gz  | jq -r ' [ .tests.test_stoichiometric_consistency.title, 
               .tests.test_reaction_mass_balance.title, 
               .tests.test_reaction_charge_balance.title, 
               .tests.test_find_disconnected.title, 
               .tests.test_find_reactions_unbounded_flux_default_condition.title ] | @tsv ')

TOTAL_SCORE=$(grep -oP "total_score\":\K[^}]*"  *_report.html)
echo -e "BIN_ID\TOTAL_SCORE\t$title" > !{sample}_!{id}_metrics.tsv
echo -e "!{id}\t${TOTAL_SCORE}\t${values}" >> !{sample}_!{id}_metrics.tsv
