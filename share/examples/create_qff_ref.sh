#script to generate quickff force fields that are used as reference for tests
for subdir in $(ls)
do
    if [ -d "$subdir" ]
    then
        cd $subdir
        rm -v *.qff
        echo "$subdir/high_he-1_uff-1"
        qff-est.py --atypes-level=high --ei-scheme=he --ei-scales=0.0,0.0,0.0  --ei-model=zero \
                                       --vdw-from=uff --vdw-scales=0.0,0.0,0.0 --vdw-model=zero gaussian.fchk gaussian.fchk.h5 > high_he-1_uff-1.qff
        echo "$subdir/high_he0_uff-1"
        qff-est.py --atypes-level=high --ei-scheme=he --ei-scales=1.0,1.0,1.0  --ei-model=harm \
                                       --vdw-from=uff --vdw-scales=0.0,0.0,0.0 --vdw-model=zero gaussian.fchk gaussian.fchk.h5 > high_he0_uff-1.qff
        echo "$subdir/high_he2_uff2"
        qff-est.py --atypes-level=high --ei-scheme=he --ei-scales=0.0,0.0,1.0  --ei-model=harm \
                                       --vdw-from=uff --vdw-scales=0.0,0.0,1.0 --vdw-model=harm gaussian.fchk gaussian.fchk.h5 > high_he2_uff2.qff
        rm pars_* system.chk
        cd ..
    fi
done
