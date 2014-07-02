"""Helper for various methods that put values in the proper format
    e.g. chebi ids should have CHEBI: stripped
"""

import util

def fix_chebi_id(id_value):
    """The Chebi ID should just be the number, we do not prepend the CHEBI:"""
    if not util.isnan(id_value):
        id_value = util.make_string(id_value)
        if id_value.startswith("CHEBI:"):
            id_value = id_value[6:]
    return id_value


def fix_inchi_parent_key(id_value):
    """The database requires the InChIKey= prepended to id"""
    if not util.isnan(id_value):

        id_value = util.make_string(id_value)

        if not id_value.startswith("InChIKey="):
            id_value = "InChIKey=" + id_value
    return id_value


def fix_inchi_parent(id_value):
    """The database requires the InChI= prepended to id"""
    if not util.isnan(id_value):

        id_value = util.make_string(id_value)

        if not id_value.startswith("InChI="):
            id_value = "InChI=" + id_value
    return id_value


def fix_all(df):
    """Loops through the dataframe fixing what we can"""
    for key in df:
        if "ChEBI_ID" == key:
            df[key] = map(fix_chebi_id, df[key])
        elif "InChi_Key_Parent" == key:
            df[key] = map(fix_inchi_parent_key, df[key])
        elif "InChi_Parent" == key:
            df[key] = map(fix_inchi_parent, df[key])

