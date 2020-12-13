#!/bin/bash
echo "Enter PODAAC username:"
read user
echo "Enter destination directory for data:"
read dd
echo $dd

python ~/read-GRACE-harmonics/scripts/podaac_grace_sync.py --user=$user --release=RL06 --directory=$dd
python ~/read-GRACE-harmonics/scripts/run_grace_date.py --directory=$dd --release=RL06 
