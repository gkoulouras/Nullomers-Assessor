################# Two-line FASTA file creator #########################
import sys

if __name__ == "__main__":
    try:
        file_path = sys.argv[1]
        seq_type = sys.argv[2]
    except:
        print("\n**The full path of the file including the filename and the filetype are required as a first argument.**\n")
        print("**Please choose either DNA or PROT to declare the sequence type.**\n")
        print("**Usage: >python fasta_formater.py {absolute_path} {sequence_type}**\n")
        raise SystemExit() 

print("******* Two-line FASTA file formatter  *******\n")


with open(file_path) as f1:
    f1lines = f1.read().splitlines()
    if (">" in f1lines[0]):
        header = f1lines[0]
    else:
        header = file_path
    if (seq_type == "DNA"):
        f1lines = [x for x in f1lines if x and (">" not in x)]
        f1lines = [item.replace("N", "") for item in f1lines]
        f1lines = [item.replace("n", "") for item in f1lines]
    elif (seq_type == "PROT"):
        f1lines = [x for x in f1lines if x and (">" not in x)]
        f1lines = [item.replace("B", "") for item in f1lines]
        f1lines = [item.replace("X", "") for item in f1lines]
        f1lines = [item.replace("Z", "") for item in f1lines]
        f1lines = [item.replace("J", "") for item in f1lines]
print("The file was read successfully")

# Open File
resultFile = open(file_path,'w')

# Write data to file
for ind, val in enumerate(f1lines):
    if (ind == 0):
        resultFile.write(header + "\n")
        resultFile.write(val.upper())
    else:
        resultFile.write(val.upper())

resultFile.close()
     
print("The file was updated successfully!\n*Success!*\n")
