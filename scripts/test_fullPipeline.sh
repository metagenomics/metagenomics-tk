OPTIONS=$@
make run_small_full_test WORK_DIR="work"  PROFILE="local" OPTIONS=" $OPTIONS "
make check
