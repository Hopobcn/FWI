COMPILATION WORKFLOW:
	- Source (or create if it is not available for your system) the environment file. Ej. "source environment_mn.sh" for the MareNostrum case.
	- Make. This should create 3 binaries.


EXECUTION WORKFLOW:

	1*Edit these two configuration files:
		- fwi_frequencies.txt is a list of wavelef frequencies. One freq. value per line. Be careful, the number of cells in the domain is multiplied by 8 everytime the frequency is incremented by 1 Hz.
		- fwi_params.txt: simulation definition parameters.
											- z dimension in meters
											- x dimension in meters
											- y dimension in meters (in real world, this is the largest)
											- vmin: min velocity propogation (it does not make sense to change this)
											- srclen: unused
											- rcvlen: unused
											- nshots: number of shots to be processed. Each shot will be managed by one slave, so, in practice this is the number of slave nodes.
											- ngradients: number of gradient iterations. Between 2 and 25 in real cases.
											- ntests: number of test iterations: Between 2 and 6 in real cases.
											- slavemem: slave's node memory in GB
											- workermem: worker's memory in GB

			** Notice that in real world FWI the number of shots will be in the order of thousands to tens of thousands.
			** The domain is divided according to the worker's memory, less memory equals more workers.

	2*Run the schedule. This will generate a ../SetupParams/fwi_schedule.txt

	3*The last column of the ../SetupParams/fwi_schedule.txt file is the number of workers needed for each case. The number of shots to process is the number of slave nodes. Edit the job file accordingly.

	4*Submit the job. Have fun :)
