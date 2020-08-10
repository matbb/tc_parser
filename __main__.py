from tc_parser import *
import sys

if __name__ == "__main__":
    if len(sys.argv) == 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help": # for report
        print("Usage: {:s} [-r|-t|-h|--help] <in file> <out file>".format(
            sys.argv[0]) )
        print("\t (default) : just convert")
        print("\t -s : print report and convert")
        print("\t -r : just print report")
        print("\t -t : use tabs as delimiters in output")
        print("\t -p : plot phase composition vs temperature")
        print("\t <out file> : if ommitted set to <in file>_out.csv")
        sys.exit(0)
    print_report = False
    write_to_disk = True
    plot_composition = False
    separator=","
    options = sys.argv[1:]
    while options[0][0] == "-":
        if options[0] == "-s":
            print_report = True
        if options[0] == "-r":
            write_to_disk = False
            print_report = True
        if options[0] == "-t":
            separator="\t"
        if options[0] == "-p":
            plot_composition = True
            from matplotlib import pyplot as plt
        options = options[1:]

    fname = options[0]
    if len(options) == 1:
        fname_out = fname.replace(".txt","") + "_out.csv"
    else:
        fname_out = options[1]

    df = parse_tc_data(fname,add_mass_and_molar_composition=True)
    if write_to_disk:
        print("Separator = " + separator)
        df.to_csv(fname_out,sep=separator)

    if print_report == False and plot_composition == False:
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

    if plot_composition:
        # Plot composition
        fig = plt.figure(figsize=(10,8))
        plt.yscale("log")
        print("Columns: ", df.columns)
        for phase in get_phases(df):
            print("Phase: ", phase)
            plt.plot( df["T"].values, df["NP("+phase+")"].values,
                     label=phase )
        plt.legend()
        plt.xlabel("T[C]")
        plt.ylabel("NP(phase)")
        xlim = plt.xlim()
        plt.xlim(xlim[1], xlim[0])
        plt.savefig(fname_out.replace(".csv","") + "_comp.png", dpi=200)
        plt.close()

        # Plot solubility of elements in each phase
        for phase in get_phases(df):
            # Plot element solubility
            fig = plt.figure(figsize=(10,8))
            plt.yscale("log")
            for el in get_elements(df):
                print("Phase: ", phase, ", element: ", el)
                plt.plot( df["T"].values, df["X("+phase+","+el+")"].values,
                        label=el )
            plt.legend()
            plt.xlabel("T[C]")
            plt.ylabel("X("+phase+",element)")
            xlim = plt.xlim()
            plt.xlim(xlim[1], xlim[0])
            plt.savefig(fname_out.replace(".csv","") +\
                        "_phasecomp_"+ phase + ".png", dpi=200)
            plt.close()


