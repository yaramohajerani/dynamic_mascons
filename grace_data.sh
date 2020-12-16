#!/bin/bash
echo "Enter PODAAC username:"
read user
echo "Enter destination directory for data:"
read dd
echo $dd

podaac_grace_sync.py --user=$user --release=RL06 --directory=$dd
run_grace_date.py --directory=$dd --release=RL06 
