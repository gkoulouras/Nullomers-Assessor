import sys

if __name__ == "__main__":
    try:
        genome_file = sys.argv[1]
        nullomers_files = sys.argv[2]
    except:
        print("\n**An error occured**\n")
        raise SystemExit() 

with open(genome_file) as gf:
    for i,lines in enumerate(gf):
        if (i != 0):
            gfline = str(lines)
    print("The file was read succesfully")
    with open(nullomers_files) as nf:
        for lns in nf:
            occur = str(gfline).count(str(lns).rstrip("\r\n"))
            print(str(lns).rstrip("\r\n") + '\t' + str(occur).rstrip("\r\n"))
