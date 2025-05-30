code_id=$1
n_frames=$2
qber=$3
REST_OF_ARGS=$4
print_update_every_n_frames=20

seconds_since_1970=$(date +%s)
out_dirname="simout_$seconds_since_1970"
mkdir $out_dirname

git rev-parse HEAD >> status.txt
git diff >> status.txt

cp status.txt "$out_dirname/results0.txt"
cp status.txt "$out_dirname/results1.txt"
cp status.txt "$out_dirname/results2.txt"
cp status.txt "$out_dirname/results3.txt"
cp status.txt "$out_dirname/results4.txt"
cp status.txt "$out_dirname/results5.txt"
cp status.txt "$out_dirname/results6.txt"
cp status.txt "$out_dirname/results7.txt"
cp status.txt "$out_dirname/results8.txt"
cp status.txt "$out_dirname/results9.txt"
cp status.txt "$out_dirname/results10.txt"
cp status.txt "$out_dirname/results11.txt"
cp status.txt "$out_dirname/results12.txt"
cp status.txt "$out_dirname/results13.txt"
cp status.txt "$out_dirname/results14.txt"
cp status.txt "$out_dirname/results15.txt"

rm status.txt

#prefix="cmake-build-ofast"
prefix="build"


echo $REST_OF_ARGS

./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 0  >> "$out_dirname/results0.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 1  >> "$out_dirname/results1.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 3  >> "$out_dirname/results3.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 4  >> "$out_dirname/results4.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 5  >> "$out_dirname/results5.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 6  >> "$out_dirname/results6.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 2  >> "$out_dirname/results2.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 7  >> "$out_dirname/results7.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 8  >> "$out_dirname/results8.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 9  >> "$out_dirname/results9.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 10 >> "$out_dirname/results10.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 11 >> "$out_dirname/results11.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 12 >> "$out_dirname/results12.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 13 >> "$out_dirname/results13.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 14 >> "$out_dirname/results14.txt" 2>&1 &
./$prefix/benchmarks_error_rate/critical_rate_simulation -cp $code_id -nf $n_frames -upn $print_update_every_n_frames -p $qber $REST_OF_ARGS --seed 15 >> "$out_dirname/results15.txt" 2>&1 &

echo "Started 16 simulations at $out_dirname"
