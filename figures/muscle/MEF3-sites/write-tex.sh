#!/bin/bash

preambule(){
echo '\begin{table}[!ht]'
echo '\centering'
echo '{\singlespacing'
echo '\begin{tabular}{|l|l|}'
echo '\hline'
}

writelines(){

    file=$1
    cat $file |\
        while read line
        do
            read name site <<< $line
            echo ''$name' & '\\texttt{$site}' \\'

        done
}

end(){
    echo '\hline'
    echo '\end{tabular}'
    echo '}'
    echo '\end{table}'
}


preambule
writelines mef3-pos-coordsavailable.dat
echo '\hline'
echo '\hline'
writelines mef3-pos-outoflearning.dat
end


echo
echo
preambule
writelines mef3-neg.dat
end
