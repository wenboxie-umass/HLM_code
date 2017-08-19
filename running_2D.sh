#!/bin/bash

#compile the c++ file
g++-7 -O3 -Wall -std=c++14 -fopenmp -Itrng-4.19 -Ltrng-4.19/src/.libs -ltrng4 HLM_KMP_2D_YL_Parallel.cpp -o 2D_Parallel

N=5
N_Max=30
N_Normal_Value=5
M=5
Step=100000
line_number=1
n=0
number_of_attempts=10
TL=1.0
TL_Max=3.0
TL_Value=1.0
TL_Step=0.5
TR=2.0
TR_Max=3.0
TR_Value=2.0
TR_Step=0.5
#Rate_Func is indicating the rate function we are using. 1 = sqrt(xy/(x+y)), 2 = sqrt(x+y)
Rate_Func=1
N_thread=8
current_attempt=1
file_name="data_2d.csv"
starting_line="N,M,TL,TR,current_attempt,"

#Check the existence of data file
if [ -f $file_name ];
	then
		rm -v $file_name
		touch $file_name
	else
		touch $file_name
fi

#Set title for each column

thread=1

while [ $thread -le $N_thread ]
do
	starting_line=$starting_line"Thread_$thread,"
	thread=$(( thread + 1 ))
done

starting_line=$starting_line"Mean,std_div"

echo $starting_line >> $file_name

unset thread

line_number=$(( line_number + 1 ))

while [ $Rate_Func -le 2 ]
do
	if [ $Rate_Func = 1 ];
		then
			echo "Rate Function: sqrt(x*y/(x+y))"
			echo "Rate Function: sqrt(x*y/(x+y)),,,,,,,,,,,,,," >> $file_name
			line_number=$(( line_number + 1 ))
		else
			echo "Rate Function: sqrt(x+y)"
			echo "Rate Function: sqrt(x+y),,,,,,,,,,,,,," >> $file_name
			line_number=$(( line_number + 1 ))
	fi
	
	#Test_1: Mean energy flux for changing length of the chain. Let N change from 2 to 100.
	echo "Test_1: Mean energy flux for changing length of the chain. Let N change from 2 to 100."
	echo "Test_1: Mean energy flux for changing length of the chain. Let N change from 2 to 100.,,,,,,,,,,,,,," >> $file_name
	line_number=$(( line_number + 1 ))
	
	while [ $N -le $N_Max ]
	do
		
		while [ $current_attempt -le $number_of_attempts ]
		do
			echo "N = $N, M = $M Current Attempt: $current_attempt"
			./2D_Parallel $N $M $N_thread $Step $Rate_Func $TL $TR $file_name $line_number $current_attempt $number_of_attempts
			#if SegFault was caught, then make another attempt
			if [ $? -ne 0 ];
				then
					echo "Seg Fault"
					continue
			fi
			current_attempt=$(( current_attempt + 1 ))
			line_number=$(( line_number + 1 ))
		done
		
		echo "-----------------------------"
		echo ",,,,,,,,,,,,,," >> $file_name
		line_number=$(( line_number + 2 ))
		
		current_attempt=1
		
		if [ $N -le 9 ];
			then
				N=$(( N + 1 ))
			else
				N=$(( N + 5 ))
		fi
	done
	
	N=$N_Normal_Value
	TR=1.5
	
	#Test_2: Mean energy flux for changing difference of boundary temperatures. Fix TL = 1 and change TR from 1.5 to 10.0.
	echo "Test_2: Mean energy flux for changing difference of boundary temperatures. Fix TL = 1 and change TR from 1.5 to 10.0."
	echo "Test_2: Mean energy flux for changing difference of boundary temperatures. Fix TL = 1 and change TR from 1.5 to 10.0.,,,,,,,,,,,,,," >> $file_name
	line_number=$(( line_number + 1 ))
	
	while (( $(bc <<< "$TR <= $TR_Max") ))
	do
		
		while [ $current_attempt -le $number_of_attempts ]
		do
			echo "N = $N, M = $M, TR = $TR Current Attempt: $current_attempt"
			./2D_Parallel $N $M $N_thread $Step $Rate_Func $TL $TR $file_name $line_number $current_attempt $number_of_attempts
			#if SegFault was caught, then make another attempt
			if [ $? -ne 0 ];
				then
					echo "Seg Fault"
					continue
			fi
			current_attempt=$(( current_attempt + 1 ))
			line_number=$(( line_number + 1 ))
		done
		
		echo "-----------------------------"
		echo ",,,,,,,,,,,,,," >> $file_name
		line_number=$(( line_number + 2 ))
		
		current_attempt=1
		
		TR=$(bc <<< "scale=3; $TR + $TR_Step")
	done
	
	TR=$(bc <<< "scale=3; $TL * 1.1")
	
	#Test_3: Mean energy flux for changing temperatures. Change TL from 1 to 100 and let TR be 1.10 times TL.
	echo "Test_3: Mean energy flux for changing temperatures. Change TL from 1 to 100 and let TR be 1.10 times TL."
	echo "Test_3: Mean energy flux for changing temperatures. Change TL from 1 to 100 and let TR be 1.10 times TL.,,,,,,,,,,,,,," >> $file_name
	line_number=$(( line_number + 1 ))
	
	while (( $(bc <<< "$TL <= $TL_Max") ))
	do
		
		while [ $current_attempt -le $number_of_attempts ]
		do
			echo "N = $N, M = $M, TL = $TL, TR = $TR Current Attempt: $current_attempt"
			./2D_Parallel $N $M $N_thread $Step $Rate_Func $TL $TR $file_name $line_number $current_attempt $number_of_attempts
			#if SegFault was caught, then make another attempt
			if [ $? -ne 0 ];
				then
					echo "Seg Fault"
					continue
			fi
			current_attempt=$(( current_attempt + 1 ))
			line_number=$(( line_number + 1 ))
		done
		
		echo "-----------------------------"
		echo ",,,,,,,,,,,,,," >> $file_name
		line_number=$(( line_number + 2 ))
		
		current_attempt=1
		
		TL=$(bc <<< "scale=3; $TL + $TL_Step")
		TR=$(bc <<< "scale=3; $TL * 1.1")
	done
	
	N=$N_Normal_Value
	TL=$TL_Value
	TR=$TR_Value
	
	Rate_Func=$(( Rate_Func + 1 ))
done 

unset N
unset M
unset Step
unset number_of_attempts
unset TL
unset TR
unset Rate_Func
unset n
unset file_name
unset starting_line
unset N_thread
unset line_number
unset current_attempt
unset N_Normal_Value
unset TL_Max
unset TL_Value
unset TR_Max
unset TR_Value