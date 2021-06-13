make run_small_full_test WORK_DIR="work" OPTIONS=" --steps.fragmentRecruitment.frhit.genomes=test/bins/small/bin.*.fa --steps.fragmentRecruitment.frhit.samples=test/reads/small/reads.tsv " PROFILE="local"  ENTRY="wFragmentRecruitment" PARAMS_FILE=example_params/fragmentRecruitment.yml
make check
