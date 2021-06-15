OPTIONS=$@
make run_small_full_test WORK_DIR="work" OPTIONS=" $OPTIONS  " PROFILE="local"  ENTRY="wFragmentRecruitment" PARAMS_FILE=example_params/fragmentRecruitment.yml
make check
