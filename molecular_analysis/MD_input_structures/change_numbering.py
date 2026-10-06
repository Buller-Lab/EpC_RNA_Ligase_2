import sys
import os

if len(sys.argv) < 3:
    print("Usage: python renumber_pdb.py <input.pdb> <offset>")
    print("Example: python renumber_pdb.py protein.pdb 100")
    print("This script adds the offset to all residue numbers and saves the result as <input>_renumbered.pdb")
    sys.exit(1)

input_pdb = sys.argv[1]
offset = int(sys.argv[2])

base, ext = os.path.splitext(input_pdb)
output_pdb = base + "_renumbered" + ext

with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
    for line in infile:
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                res = int(line[22:26])
                new_res = res + offset
                new_line = line[:22] + f"{new_res:4d}" + line[26:]
                outfile.write(new_line)
            except ValueError:
                outfile.write(line)
        else:
            outfile.write(line)

print(f"Renumbering completed successfully!")
print(f"Output saved as: {output_pdb}")