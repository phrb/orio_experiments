#!/usr/bin/env python

import dataset, csv

database = dataset.connect("sqlite:///search_space.db")

experiments = database["experiments"]

with open("search_space.csv", "w", newline = "") as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames = experiments.find_one(id = 1).keys())

    writer.writeheader()
    writer.writerows(experiments)
