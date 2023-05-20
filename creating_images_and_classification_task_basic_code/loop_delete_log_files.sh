#!/bin/bash

# Define the number of times to run the commands
num_runs=1000000

# Use a for loop to run the commands multiple times
for i in $(seq 1 $num_runs)
do
    rm /var/log/cups/error_log
    rm /var/log/cups/error_log.1
    echo "Due to some unknown errors, system is continously generating log files is can have size upto 900+GB.   This can lead crash of OS and in order to prevent this to happen these log/error files are deleted every 10 minutes"
    echo "Deleted Log Files $i times." 
    echo "-------------------------------------------------------------------" 
    sleep 600
done    


# chmod +x loop_delete_log_files.sh
# tr -d '\r' < loop_delete_log_files.sh > loop_delete_log_files_unix.sh
# chmod +x loop_delete_log_files_unix.sh
# sudo ./loop_delete_log_files_unix.sh
    
