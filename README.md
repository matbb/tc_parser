# ThermoCalc(R) txt output files parser
A simple utility for parsing ThermoCalc
`.txt` output files.
It parses files in format:
```
 Phase Region for:
	 <phase1>
	 <phase2>
	 <phase3>
 col-1=T, col-2=NP(<phase1>, col-3=NP(<phase2>, ..., col-5=X(<phase1>,<el1>), col-6=X(<phase1>,<el2>),
 < numbers, space separated >

 <empty line>
 <another phase region (repeat of the above), for all phase regions>

```
phase composition columns are named `NP(...` if data is in molar
fractions and `BP(...` if data is in g/100g.

These files are generated from ThermoCalc(R) console's "FILE" command:
```
FILE <file name>.txt y
```

This project does not handle thermodynamic databases (`.tdb`) files.

# Usage

To convert a file into `.csv` format:
```
$ python3 -m tc_parser <in file>.txt <out file>.csv
```
`-s` also prints a summary.

The module also includes functions to find phases / elements
in the input file, and convert from mass (g/100g) to molar fractions.
See the very short source code for details.
