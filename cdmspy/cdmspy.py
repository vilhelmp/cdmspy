import requests
import bs4
import urllib.parse
from difflib import SequenceMatcher as SM
import numpy as np
from astropy.table import Table
import astropy.units as u
import astropy.constants as c

# define the base URLs to use
BASEURL= "http://www.astro.uni-koeln.de"
FORMURL = urllib.parse.urljoin(BASEURL, "/cgi-bin/cdmssearch")

# GET THE molecular line list
# do this on import
# Temporary workaround: added "_" infront
# of names to make hidden in namespace.
_searchpage = requests.get(FORMURL)
_searchpage.close()
_searchform = bs4.BeautifulSoup(_searchpage.content, "lxml")
_mollist_html = _searchform.find_all("select")[0]
molecules = [i for i in _mollist_html.stripped_strings]
mol_lists = np.array([[i[:6], i[7:]] for i in molecules])
mol_ids = mol_lists.T[0]
mol_names = mol_lists.T[1]
# now we can search

def find_molecules(tofind, lim=0.8):
    # remove "-" from each species name string
    # this is to remove unecessary symbols before
    # comparing the strings.
    # TODO: What about species strings with
    #"v=1" and similar in their name.
    ignore = lambda x: x == "-"
    scores = [SM(ignore,
                 tofind,
                 i.replace("-", "")).ratio() for i in mol_names]
    # score should be above 80% (i.e. 0.8) to qualify
    return np.array(molecules)[np.array(scores)>=lim]


def query(freqs=None,
          molecules=None):
    payload = dict(MinNu=freqs[0],
                   MaxNu=freqs[1],
                   UnitNu="GHz",
                   StrLim=-10,
                   Molecules=molecules,
                   temp=0, output="text",
                   sort="frequency",
                   mol_sort_query="tag",
                   logscale="yes",
                   but_action="Submit")
    postrequest = requests.post(FORMURL, data=payload)
    postrequest.close()
    soup = bs4.BeautifulSoup(postrequest.content, "lxml")
    link = soup.find_all('a')[0].get('href')
    newurl = urllib.parse.urljoin(BASEURL, link)
    resultstable = requests.get(newurl)
    resultstable.close()
    newsoup = bs4.BeautifulSoup(resultstable.content, "lxml")
    asciitable = newsoup.find_all('pre')[0].get_text()
    lines = parse_results_table(asciitable)
    return lines




def parse_results_table(asciitable):
    cdms_colnames = (
    'freq_rest', 'freqerr', 'aij', 'dofrot', 'elow_cm', 'gup', 'tag', 'qenq', 'qnum1', 'qnum2', 'species')
    # (F13.4,F8.4, F8.4,  I2,F10.4,  I3,  I7,    I4,  6I2,  6I2)
    cdms_colstarts = (0, 13, 24, 35, 37, 47, 50, 57, 61, 72, 89)
    lines = Table.read(asciitable,
                       format='ascii.fixed_width_no_header',
                       names=cdms_colnames,
                       col_starts=cdms_colstarts,
                       )
    # give out some units in one go, so that we don't run it twice on the same table
    lines['freq_rest'] = lines['freq_rest'] * 1e-3
    lines['freq_rest'].unit = u.GHz
    lines['freqerr'] = lines['freqerr'] * 1e-3
    lines['freqerr'].unit = u.GHz

    lines['elow_cm'].unit = u.cm ** -1

    # if elow in cm**-1 (CDMS table), then transfer into K
    lines['elow'] = (lines['elow_cm'].quantity * c.c * c.h / c.k_B).decompose().to(u.K)

    # calculate the E_up (in Kelvin) from the E_low (in Kelvin)
    lines['eup'] = lines['elow'] + ((c.h * lines['freq_rest'].quantity) / c.k_B).decompose()
    lines['eup_cm'] = (lines['eup'].quantity * c.k_B / (c.c * c.h)).decompose().to(1 / u.cm)

    return lines