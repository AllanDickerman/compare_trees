import sys

if len(sys.argv) < 2:
    print("python selectNperCategory.py F N K")
    print("select N random items for each value of column K from file F")
    print("output to stdout")
    print("K is 0-based, default 1")
    sys.exit(0)

F = open(sys.argv[1])
num_per = int(sys.argv[2])
key_column = 0
if len(sys.argv) > 3:
    key_column = int(sys.argv[3])

header_line = F.readline()
headers = header_line.rstrip().split("\t")
#print("header fields = "+"\n".join(headers))
#print(f"selecting {num_per} from column {key_column} = {headers[key_column]}")
print(header_line, end='')
data = {}
for line in F:
    fields = line.rstrip().split("\t")
    key = fields[key_column]
    if not key in data:
        data[key] = set()
    data[key].add(line)

for key in data:
    x = list(data[key])[:num_per]
    print("".join(x), end='')
