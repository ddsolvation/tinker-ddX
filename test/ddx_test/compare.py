#!/usr/bin/python3
import sys

threshold = 1e-5

def read_log(path):
    energy = 0.0
    with open(path, "r", encoding="utf-8") as log:
        for line in log:
            if "Solvation energy (Hartree):" in line:
                tokens = line.split()
                energy = float(tokens[3])
                break
    return energy

output_file = sys.argv[1]
ref_file = sys.argv[2]
if len(sys.argv) == 4:
    if sys.argv[3] == "loose":
        threshold = 5e-2

energy = read_log(output_file)
ref_energy = read_log(ref_file)

rel_error = abs(energy - ref_energy)/abs(ref_energy)

if rel_error > threshold:
    print(f"Energy:         {energy:20.10f}")
    print(f"Ref. energy:    {ref_energy:20.10f}")
    print(f"Rel. error:     {rel_error:20.10f}")
    quit(1)
