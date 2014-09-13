#!/bin/bash
#Usage : ./test.sh test_number fixed_parameters number_of_tests start_n
#\"${NAMEFILE}.txt\" using 1:4 title 'temps_ordm' with linespoints,\

N=1
SAGE=${SAGE:-sage}
GNUPLOT=${GNUPLOT:-gnuplot}
CMD=batch.py


if [[ "$1" == 1 ]];then
    while [ -f "data_test1_p_fixed${N}.txt"  ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test1_p_fixed${N}"

    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt

    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"extension degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'temps_m' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'temps_E' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:4 title 'temps_ordm' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:5 title 'temps_period' with linespoints"| $GNUPLOT

elif [[ "$1" == 2 ]];then
    while [ -f "data_test2_n_fixed${N}.txt"  ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test2_n_fixed${N}"

    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt

    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"caractéristique\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'temps_m' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'temps_E' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:4 title 'temps_ordm' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:5 title 'temps_period' with linespoints"| $GNUPLOT
elif [[ "$1" == 3 ]];then
    while [ -f "data_test3_p_fixed${N}.txt"  ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test3_p_fixed${N}.txt"

    echo 'TODO'
elif [[ "$1" == 4 ]];then
    while [ -f "data_test4_n_fixed${N}.txt"  ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test4_n_fixed${N}"

    echo 'TODO'

elif [[ "$1" == 5 ]];then
    while [ -f "data_test5_p_fixed${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test5${N}.text"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_cycl' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_ell' with linespoints"| $GNUPLOT
elif [[ "$1" == 6 ]];then
    while [ -f "data_test6_p_fixed${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test6_ext${N}.text"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_cycl' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_ell' with linespoints"| $GNUPLOT
elif [[ "$1" == 7 ]];then
    while [ -f "data_test7_p_fixed${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test7_flint${N}.text"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_cycl' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_ell' with linespoints"| $GNUPLOT
elif [[ "$1" == 8 ]];then
    while [ -f "data_test8_ext_flint${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test8_ext_flint${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_cycl' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_ell' with linespoints"| $GNUPLOT
elif [[ "$1" == 9 ]];then
    while [ -f "data_test9_cmpflint${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test9_cmpflint${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
         
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_ell' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_ell_flint' with linespoints"| $GNUPLOT
elif [[ "$1" == 10 ]];then
    while [ -f "data_test10_cmpflint${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test10_cmpflint${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_ell' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_ell_flint' with linespoints"| $GNUPLOT
elif [[ "$1" == 11 ]];then
    while [ -f "data_test11_cmpellFFH${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test11_cmpellFFH${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_ell' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_naive' with linespoints"| $GNUPLOT
elif [[ "$1" == 12 ]];then
    while [ -f "data_test12_cmpellFFH${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test12_cmpellFFH${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_ell' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_naive' with linespoints"| $GNUPLOT
elif [[ "$1" == 13 ]];then
    while [ -f "data_test13_cmpellFFH${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test13_cmpellFFH${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_ell' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_all' with linespoints"| $GNUPLOT
elif [[ "$1" == 14 ]];then
    while [ -f "data_test14_cmpellFFH${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test14_cmpellFFH${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_ell_flint' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_all' with linespoints"| $GNUPLOT
elif [[ "$1" == 15 ]];then
    while [ -f "data_test15_testcycl${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test15_testcycl${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_cycl' with linespoints"| $GNUPLOT
elif [[ "$1" == 16 ]];then
    while [ -f "data_test16_testcycl_fixed_degree${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test16_testcycl_fixed_degree${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"caractéristique\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_cycl' with linespoints"| $GNUPLOT
elif [[ "$1" == 17 ]];then
    while [ -f "data_test17_testell${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test17_testell${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_ell' with linespoints"| $GNUPLOT
elif [[ "$1" == 18 ]];then
    while [ -f "data_test18_testell_fixed_degree${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test18_testell_fixed_degree${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"caractéristique\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_ell' with linespoints"| $GNUPLOT
elif [[ "$1" == 19 ]];then
    while [ -f "data_test19_testnorm${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test19_testnorm${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_coeff' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_isom' with linespoints"| $GNUPLOT
elif [[ "$1" == 20 ]];then
    while [ -f "data_test20_cmpellcycl_nocase${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test20_cmpellcycl_nocase${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_cycl' with linespoints,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_ell' with linespoints"| $GNUPLOT
elif [[ "$1" == 21 ]];then
    while [ -f "data_test21_cmpcyclFFH_nocase${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test21_cmpcyclFFH_nocase${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_cycl' with points,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_naive' with linespoints"| $GNUPLOT
elif [[ "$1" == 22 ]];then
    while [ -f "data_test22_cmpcyclall_nocase${N}.txt" ]
    do
        N=$((N+1))
    done

    NAMEFILE="data_test22_cmpcyclall_nocase${N}"
    
    echo 'Generating the data file...'
    $SAGE $CMD $1 $2 $3 $4> ${NAMEFILE}.txt
            
    echo 'Generating the graph...'

    echo "set term jpeg
    set output \"${NAMEFILE}.jpg\"
    set xlabel \"degré n\"
    set ylabel \"temps(sec)\"
    plot \"${NAMEFILE}.txt\" using 1:2 title 'T_cycl' with points,\
     \"${NAMEFILE}.txt\" using 1:3 title 'T_All' with linespoints"| $GNUPLOT
fi

