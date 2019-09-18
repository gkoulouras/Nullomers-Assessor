################# oneline fasta file creator #########################
import sys

if __name__ == "__main__":
    try:
        file_path = sys.argv[1]
    except:
        print("\n**The full path of the file including the filename and the filetype are required to be given as an argument.**\n**Please correct and try again.**\n")
        raise SystemExit() 

print("******* Multiple lines to one line fasta formatter  *******\n")


with open(file_path) as f1:
    f1lines = f1.read().splitlines()
    if (">" in f1lines[0]):
        header = f1lines[0]
    else:
        header = file_path
    f1lines = [x for x in f1lines if x and (">" not in x)] #and ("X" not in x) and ("Z" not in x) and ("J" not in x)
print("The file was read successfully")

# Open File
resultFile = open(file_path,'w')

# Write data to file
for ind, val in enumerate(f1lines):
    if (ind == 0):
        resultFile.write(">" + header + "\n")
        resultFile.write(val.upper())
    else:
        resultFile.write(val.upper())

resultFile.close()
     
print("The file was read successfully!\n*Success!*\n")

