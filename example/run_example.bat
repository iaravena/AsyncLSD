echo ""
echo "Solving stochastic unit commitment for WECC instance of a typical weekday in Autumn with 10 scenarios"
mpiexec -n 3 ..\AsyncLSD.exe "WECC_SUC_AutumnWD_0010" ".\" "WECC_SUC_AutumnWD_0010" "example_options_file.txt"
