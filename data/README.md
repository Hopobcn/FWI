
- fwi_params.txt:

    ```
    lenz: physical domain size in z dimension (meters)
    lenx: physical domain size in x dimension (meters)
    leny: physical domain size in y dimension (meters)
    vmin: minimum wavelet velocity (expected) (m/s)
    srclen: (not in use)
    rcvlen: (not in use)
    nshots: number shots
    ngrads: number of gradient iterations (typically, under 20)
    ntests: number of test iterations (typically, between 7 and 10)
    slavemem: memory (in GB) of slave nodes
    workermem: memory (in GB) of worker nodes
    shotdirectory: path for storing local IO output
    ```

- fwi_frequency.txt: contains a list of frequencies to simulate:

    ```
    1.00
    2.00
    etc...
    ```
- fwi_shcedule.txt: generated file using `fwi-sched-generator` and `fwi_params.txt` + `fwi_frequencies.txt` files

    ```
    nfreqs
    nshots
    ngrads
    ntests
    outputfolder
    freq forws backs stacki dt dz dx dy dimmz dimmx dimmy ppd nworkers
    ```

    
When `PERFORM_IO` is enabled, data is initialized from files present in `data/inputmodels`. 
Those files are generated using the `fwi-data-generator`.
