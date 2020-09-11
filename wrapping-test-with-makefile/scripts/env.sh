FILEPATH=`realpath ${BASH_SOURCE[0]}`
DIR=`dirname $FILEPATH`
export PATH=$DIR:$PATH
export PYTHONPATH=/home/nnthuynh/Bureau/Stage-Wrapping-Python/wrapping-test/wrapping-test-with-makefile/build/wrappers

alias run_twinexp='python $DIR/run_twinexp.py $*'
