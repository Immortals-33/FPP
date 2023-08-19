import csv

txt_filename = "test.txt"
csv_filename = "recovery.csv"

data = []

with open(txt_filename, 'r') as txt_file:
    for line in txt_file:
        values = line.strip().split()
        if len(values) == 2:
            data.append((int(values[0]), float(values[1])))

with open(csv_filename, 'w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(["Value", "Score"])
    csv_writer.writerows(data)

print(f"Data written to {csv_filename}")

