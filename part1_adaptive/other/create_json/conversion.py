import json
import csv
from datetime import datetime


now = datetime.now()
timestamp = now.strftime("%Y%m%d%H%M%S")

CSV_FILEPATH = "assets/trials.csv"
JSON_OUTPUTFILE = "assets/trials_" + timestamp + ".json"
EXPECTED_COLUMN_NAMES = ['Label','Name','TargetEmo','ComparisonEmo','Distance','Label1','Label2',
                         'Soundtrack','Level']


trials = {}

# read csv and convert into arranged python dictionary
with open(CSV_FILEPATH) as csv_file:
    dict_reader = csv.DictReader(csv_file)

    # read the rows of csv and create the expected dictionary
    for row in dict_reader:
        # assert that you got correct column names from csv, if wrong the json file won't be created
        csv_column_names = list(row.keys())
        assert set(csv_column_names) == set(EXPECTED_COLUMN_NAMES)

        level = row["Level"]
        emotion = row["Label1"]
        if level not in trials.keys():
            trials[level] = {}

        if emotion not in trials[level].keys():
            trials[level][emotion] = []

        trials[level][emotion].append(row)

# write created result into a file
json_trials = json.dumps(trials)
with open(JSON_OUTPUTFILE, 'w') as json_file:
    json.dump(trials, json_file)
