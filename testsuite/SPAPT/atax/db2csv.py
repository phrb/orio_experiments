#!/usr/bin/env python3

import dataset, csv

database = dataset.connect("sqlite:///results.db")

experiments = database["results"]

with open("results.csv", "w", newline = "") as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames = experiments.find_one(id = 1).keys())

    writer.writeheader()
    writer.writerows(experiments)
