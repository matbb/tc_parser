#!/usr/bin/env python3
# coding: utf-8

# Thermocalc files parsing utility
# ==================
# 
# Provides:
# * parse a file into a pandas dataframe
# * report on composition of the alloy based on the parsed dataframe
# 
# Report contains:
# * First all present phases
#     TC phase name | chem. comp (avg) | chem. comp (with ranges)
# * By region: all phases in the region with chem. comp (avg and ranges)

import regex as re
import pandas as pd
import numpy as np
import io

try:
    set_trace
except NameError:
    from ipdb import set_trace

molar_col = "NP({phase:s})"
mass_col = "BP({phase:s})"
fraction_col = "X({phase:s},{element:s})"

def element_to_symbol(el):
    return el[0] + el[1:].lower()

def has_molar_composition( columns ):
    """ Returns true if compositin columns in molar fractions exist """
    if type(columns) == type(pd.DataFrame):
        columns = columns.columns
    return True in [ col[0:3] == "NP(" for col in  columns ]

def has_mass_composition( columns ):
    """ Returns true if compositin columns in mass g/100g exist """
    if type(columns) == type(pd.DataFrame):
        columns = columns.columns
    return True in [ col[0:3] == "BP(" for col in  columns ]

def get_phases( columns ):
    if type(columns) == type(pd.DataFrame):
        columns = columns.columns
    if has_molar_composition(columns):
        return [ col[3:-1] for col in columns if col[0:3] == molar_col[0:3] ]
    else:
        return [ col[3:-1] for col in columns if col[0:3] == mass_col[0:3] ]

def get_elements( columns ):
    if type(columns) == type(pd.DataFrame):
        columns = columns.columns
    elements = []
    return sorted( set( [ col.split(",")[1][:-1] for col in columns if
        col[0:2] == fraction_col[0:2] ] ) )

def convert_composition( df, from_molar_to_mass_per_100g=True):
    import periodictable as pt
    phases = get_phases(df.columns)
    elements = get_elements(df.columns)
    if from_molar_to_mass_per_100g:
        from_col, to_col = molar_col, mass_col
        f_convert = lambda el, val : val * pt.elements.symbol(element_to_symbol(el)).mass
        total = 100.0
    else:
        from_col, to_col = mass_col, molar_col
        f_convert = lambda el, val : val / pt.elements.symbol(element_to_symbol(el)).mass
        total = 1.0

    sum_phases = np.zeros(len(df))
    for phase in phases:
        # 1. Calculate mass of 1-mol of the phase
        mass = np.zeros(len(df))
        for el in elements:
            col = fraction_col.format(phase=phase,element=el)
            mass += df[col] * pt.elements.symbol(element_to_symbol(el)).mass

        if from_molar_to_mass_per_100g:
            new_val = df[from_col.format(phase=phase)] * mass
        else:
            new_val = df[from_col.format(phase=phase)] / mass

        df[to_col.format(phase=phase)] = new_val

        sum_phases += np.where(np.isnan(new_val),0.0,new_val)

    # 2. normalize the total
    for phase in phases:
        df[to_col.format(phase=phase)] /= sum_phases * total

    return df

def convert_composition_from_molar_to_mass_per_100g(composition):
    import periodictable as pt
    mass_composition = {}
    for el,mol in composition.items():
        atomic_mass = pt.elements.symbol(element_to_symbol(el)).mass
        mass_composition[el] = mol * atomic_mass
    total_g = np.sum(list(mass_composition.values()))
    mass_composition_100g = { el : m/total_g*100.0 for el, m in mass_composition.items() }
    return mass_composition_100g

def convert_composition_from_mass_per_100g_to_molar(composition):
    import periodictable as pt
    mol_composition = {}
    for el, mass in composition.items():
        atomic_mass = pt.elements.symbol(element_to_symbol(el)).mass
        mol_composition[el] = mass / atomic_mass
    sum_mol = sum( [ val for val in mol_composition.values() ] )
    mol_composition_x = { el : mol/sum_mol for el, mol in mol_composition.items() }
    return mol_composition_x

def get_composition( df, return_range=False ):
    composition = {}
    ranges = {}
    elements = get_elements(df.columns)
    phases = get_phases(df.columns)
    if not has_molar_composition(df.columns):
        convert_composition(df,from_molar_to_mass_per_100g=False)

    for el in get_elements(df.columns):
        el_content = np.zeros(len(df))
        for phase in phases:
            el_phase_frac = df[molar_col.format(phase=phase)]
            el_phase_el = df[fraction_col.format(
                phase=phase,element=el)]
            el_this_phase = el_phase_frac * el_phase_el
            el_content += np.where( np.isnan(el_this_phase), 0.0,
                    el_this_phase )
        composition[el] = np.average(el_content)
        if return_range:
            ranges[el] = np.max(el_content) - np.min(el_content)
    if return_range:
        return composition,ranges
    return composition

def get_phase_composition( df, phase, return_range=False ):
    composition = {}
    ranges = {}
    for el in get_elements(df.columns):
        col = fraction_col.format(phase=phase,element=el)
        composition[el] = np.nanmean(df[col])
        if return_range:
            ranges[el] = np.nanmax( df[col] ) - np.nanmin( df[col] )
    if return_range:
        return composition,ranges
    return composition


def parse_tc_data(
        filename,
        keep_T=True,
        phase_region_keep="last", # both, first or last
        # add mass/molar composition columns if missing
        add_mass_and_molar_composition=True,
        ):
    """
    Parses a thermocalc txt output data file into a pandas table.

    Thermocalc file is of the format:

     Phase Region for:
         <phase1>
         <phase2>
         <phase3>
     col-1=T, col-2=NP(<phase1>, col-3=NP(<phase2>, ..., col-5=X(<phase1>,<el1>), col-6=X(<phase1>,<el2>),
     < numbers, space separated >

     <empty line>
     <another phase region (repeat of the above), for all phase regions>

    phase composition columns are named NP(... if data is in molar
        fractions and BP(... if data is in g/100g
    Element names are in capital letters.
    """
    with open(filename, "rt") as f:
        txt = f.read()
        regions = txt.split("\n Phase Region for:")[1:]
        df = None
        i_reg = 1
        for reg in regions:
            icolumns = reg.find("col-1")
            idata = reg.find("\n", icolumns)
            columns = reg[icolumns:idata]
            columns = [x.split("=")[1][:-1]
                       for x in re.findall("col-[0-9]*=[^ ]*,", columns)]
            data = reg[idata+1:]
            dff = pd.read_table(io.StringIO(data), names=columns, sep="\s+")
            dff["region"] = i_reg
            i_reg += 1
            dff["Tx100"] = np.rint(dff["T"] * 100).astype(int)
            dff.set_index("Tx100", inplace=True)
            dff.drop(["T"], inplace=True, axis=1)
            dff.sort_index(inplace=True, ascending=False)
            if df is None:
                df = dff
            else:
                df = pd.concat([df, dff])
        if phase_region_keep != "both":
            df = df.loc[~df.index.duplicated(keep=phase_region_keep)]
    df.sort_index(inplace=True)
    if keep_T:
        df["T"] = df.index / 100.0
    if add_mass_and_molar_composition:
        if not has_mass_composition(df.columns):
            convert_composition(df,from_molar_to_mass_per_100g=True)
        if not has_molar_composition(df.columns):
            convert_composition(df,from_molar_to_mass_per_100g=False)
    return df



def convert_composition_from_molar_to_mass_per_100g(composition):
    import periodictable as pt
    mass_composition = {}
    for el,mol in composition.items():
        atomic_mass = pt.elements.symbol(element_to_symbol(el)).mass
        mass_composition[el] = mol * atomic_mass
    total_g = np.sum(list(mass_composition.values()))
    mass_composition_100g = { el : m/total_g*100.0 for el, m in mass_composition.items() }
    return mass_composition_100g

def convert_composition_from_mass_per_100g_to_molar(composition):
    import periodictable as pt
    mol_composition = {}
    for el, mass in composition.items():
        atomic_mass = pt.elements.symbol(element_to_symbol(el)).mass
        mol_composition[el] = mass / atomic_mass
    sum_mol = sum( [ val for val in mol_composition.values() ] )
    mol_composition_x = { el : mol/sum_mol for el, mol in mol_composition.items() }
    return mol_composition_x

if __name__ == "__main__":
    import sys
    if sys.argv[1] == "-h" or sys.argv[1] == "--help": # for report
        print("Usage: {:s} [-r|-h|--help] <in file> <out file>".format(
            sys.argv[0]) )
        print("\t (default) : just convert")
        print("\t -s : print report and convert")
        print("\t -r : just print report")
        print("\t <out file> : if ommitted set to <in file>_out.csv")
        sys.exit(0)
    print_report = False
    write_to_disk = True
    options = sys.argv[1:]
    while options[0][0] == "-":
        if options[0] == "-s": # for report
            print_report = True
        if options[0] == "-r": # for report
            write_to_disk = False
            print_report = True
        options = options[1:]

    fname = options[0]
    if len(options) == 1:
        fname_out = fname.replace(".txt","") + "_out.csv"
    else:
        fname_out = options[1]

    df = parse_tc_data(fname,add_mass_and_molar_composition=True)
    if write_to_disk:
        df.to_csv(fname_out)

    if print_report == False:
        sys.exit(0)

    print("Elements: ", get_elements(df))
    print()
    print("Phases: ", get_phases(df))
    print()
    print("Composition:")
    comp, comp_ranges = get_composition(df,return_range=True)
    comp_mass = convert_composition_from_molar_to_mass_per_100g(comp)
    print("{:<4s} {:>12s} | {:>12s} [{:12s}]".format(
        "el", "mol", "g/100g", "mol(max-min)" ))
    for el in get_elements(df):
        print("{:<4s} {: 12.6f} | {: 12.6f} [{: 12.6f}]".format(
            el, comp[el], comp_mass[el], comp_ranges[el] ))

