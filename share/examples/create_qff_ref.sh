#script to generate quickff force fields that are used as reference for tests
for subdir in $(ls)
do
    if [ -d "$subdir" ]
    then
        cd $subdir
        echo "$subdir/high_0_he"
        qff-est.py --atypes-level=high --ei-scheme=he --ei-scales=1.0,1.0,1.0 --ei-model=harm gaussian.fchk gaussian.fchk.h5 > high_0_he.qff
        echo "$subdir/high_1_he"
        qff-est.py --atypes-level=high --ei-scheme=he --ei-scales=0.0,1.0,1.0 --ei-model=harm gaussian.fchk gaussian.fchk.h5 > high_1_he.qff
        echo "$subdir/high_2_he"
        qff-est.py --atypes-level=high --ei-scheme=he --ei-scales=0.0,0.0,1.0 --ei-model=harm gaussian.fchk gaussian.fchk.h5 > high_2_he.qff
        echo "$subdir/high_3_he"
        qff-est.py --atypes-level=high --ei-scheme=he --ei-scales=0.0,0.0,0.0 --ei-model=harm gaussian.fchk gaussian.fchk.h5 > high_3_he.qff
        echo "$subdir/high_-1_he"
        qff-est.py --atypes-level=high --ei-scheme=he --ei-scales=0.0,0.0,0.0 --ei-model=zero gaussian.fchk gaussian.fchk.h5 > high_-1_he.qff
        rm pars_* system.chk
        cd ..
    fi
done
