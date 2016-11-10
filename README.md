# cdmspy

Query the CDMS database and read ASCII tables from CDMS.

The Cologne Database for Molecular Spectroscopy
http://www.astro.uni-koeln.de/cdms

# Why?

The CDMS web interface is allways updated first. Any changes to their database goes directly into the web search. The VAMDC database takes a couple of weeks to update.


# Install

Open terminal in directory after cloning it.
Then run

    python install setup.py


# How?

## Searching in frequency interval and for certain species

Once in python import the library with

    import cdmspy

Then to search for a specific molecule you can use the built in function to get the IDs etc.

    cdmspy.find_molecules("HDCO")

will return 

    array(['031501 HDCO', '032507 HDC-13-O', '033510 HDCO-18', '043507 HDC2O'], 
      dtype='<U33')

so lets say we want to search for all HDCO isotopologues in a certain frequency range. We first grab the names of them with the search function. [:-1] here is to not include HDC2O as seen above.

    mols = cdmspy.find_molecules("HDCO")[:-1]
    table = cdmspy.query(freqs=[200, 230], molecules=mols)

this will return an Astropy Table with all the information on lines of isotopologues of HDCO between 200 and 230 GHz.

    table.show_in_browser()

will pop-up a window with the results. The table also contains the Eu and El in Kelvin, instead of just cm-1. 



## Get full line list of specific species

Example usage

    import cdmspy
    results = cdmspy.get_entry('032507 HDC-13-O')
    results.show_in_browser()


## Get partition function for specific species

Example usage


## TODO

TODO: implement parsing of Astropy Units as intput.
TODO: use molecular ID to get the partition function values and return a function to estimate it at temperature T.

