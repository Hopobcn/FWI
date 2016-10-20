FWI has two input files

- fwi_params.txt with contents:

    ```
    lenz (meters)
    lenx (meters)
    leny (meters)
    vmin (m/s)
    srclen (units - not used)
    rcvlen (units - not used)
    nshots (units - not used)
    ngrads (units - not used)
    ntests (units - not used)
    slavemem (GB - not used)
    workermem (GB - not used)
    shotdirectory (char - not used)
    ```

- fwi_frequency.txt that contains a list of frequencies to simulate:

    ```
    1.00
    2.00
    etc...
    ```

When PERFORM_IO is enabled, data is initialized from files present in `data/inputmodels`. Those files are generated using the `fwi-data-generator`.
