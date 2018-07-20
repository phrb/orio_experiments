#!/usr/bin/env python2

import dataset, csv

database = dataset.connect("sqlite:///results.db")

experiments = database["results"]

with open("results.csv", "w") as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames = experiments.find_one(id = 1).keys())

    writer.writeheader()
    writer.writerows(experiments)

database = dataset.connect("sqlite:///search_space.db")

experiments = database["experiments"]

with open("search_space.csv", "w") as csv_file:
    writer = csv.DictWriter(csv_file, fieldnames = experiments.find_one(id = 1).keys())

    writer.writeheader()
    writer.writerows(experiments)
