import glob
import sys



outlines = []
for tbl1file in glob.glob(sys.argv[1]):
    sys.stdout.write("Searching file {}\n".format(tbl1file))
    query_names = set()
    for i in sys.argv[2:]:
        query_names.add(i.upper())
    with open(tbl1file) as f:
        f.readline()
        for line in f:
            labnum, name, idqld, dob, collected, epi, lineage = line.rstrip().split(",")
            if name.startswith('"'):
                name = name[1:]
            if name.endswith('"'):
                name = name[:-1]
            names = set(name.upper().split())
            if labnum.startswith('"'):
                labnum = labnum[1:]
            if labnum.endswith('"'):
                labnum = labnum[:-1]
            if idqld.startswith('"'):
                idqld = idqld[1:]
            if idqld.endswith('"'):
                idqld = idqld[:-1]
            if query_names.issubset(names) or labnum in sys.argv[2:] or idqld in sys.argv[2:] or "QLD" + idqld in sys.argv[2:]:
                print(labnum, "Ding")
                outlines.append(line)

sys.stdout.write("########## ANALAYSIS DONE ###############\n")
for i in outlines:
    sys.stdout.write(i)