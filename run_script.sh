
# Spectrum configuration 

# Punctures to calculate
# It begins at COUNTER+1
COUNTER=5
TOTAL=10

# The lattice type will have the following convention

# For "linear" lattices letters 'a', 'b', 'c', ...
#   this letter depends on the label of the top right puncture of 1
#   'a': for top right puncture 1, i.e. lattice with adjacencies
#   'b', 'c', ...: for top right 2, 3, ...

# For squared matrices use 'sq'
LTTYPE='b'
TR=2

# Hilbert space and Hamiltonian
ONLYINT=0
KVALUE=1
KTILD=0

# Range alpha (0 means default)
DEFALP=1

POINTS=225

#Save plots(0), or show the plots(1)
SVPLSW=0


echo ------- BEGIN SPECTRUM CALCULATION --------
# Loops over Number of punctures

while [  $COUNTER -lt $TOTAL ]; do
        let COUNTER=COUNTER+1 # see, it begins at COUNTER+1

        echo -------------------------------
        echo ----- Changing config --------
        python change_config.py $COUNTER $LTTYPE $TR $ONLYINT $KVALUE $KTILD $DEFALP $POINTS $SVPLSW

        # echo -------------------------------
        # echo ----- Loop $COUNTER out of $TOTAL --------
        # echo -------------------------------
        # echo --- Evaluating 'torus_arbsize.py' ---
        # python torus_arbsize.py 

        # echo -------------------------------
        # echo --- Evaluating 'eigencalc.py' ---
        # python eigencalc.py

        echo -------------------------------
        echo --- Evaluating 'process.py' ---
        python process.py 

        echo -------------------------------
        echo --- Evaluating 'plotting.py' ---
        #python plotting.py

        echo ---------------
        echo --- Save to CSV files
        python ordergeteigendata.py
        
done
