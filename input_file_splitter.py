"""Helper script to split MIBiG submission portal input files for enzymes in "Tailoring"""

import json
import os
import re
import sys


def generate_unique_filename(json_file_path, digits=7):
    directory = os.path.dirname(json_file_path)
    counter = 1
    while True:
        filename = f"MITE{counter:0{digits}d}.json"
        filepath = os.path.join(directory, filename)
        if not os.path.exists(filepath):
            return filepath
        counter += 1


def write_json(data, json_file_path):
    with open(generate_unique_filename(json_file_path), "w") as outfile:
        json.dump(data, outfile, indent=4)


def split_files(data, json_file_path):
    changelog_list = data.get("Changelog")
    tailoring_list = data.get("Tailoring")

    if tailoring_list is None:
        print(f"Error: File '{json_file_path}' has no tailoring info - pass.")
        sys.exit(1)

    tailoring_dict = {}
    for entry in tailoring_list:
        tailoring_dict[entry[0]] = entry[1]

    enzymes = set()
    for key, _value in tailoring_dict.items():
        match = re.match(r"^enzymes-(\d+)", key)
        if match:
            enzymes.add(match.group(1))

    for enzyme in enzymes:
        payload = []
        for key, value in tailoring_dict.items():
            match = re.match(rf"^enzymes-{enzyme}", key)
            if match:
                new_key = re.sub(rf"^enzymes-{enzyme}", r"enzymes-0", key, count=1)
                payload.append([new_key, value])

        write_json(
            data={"Changelog": changelog_list, "Tailoring": payload},
            json_file_path=json_file_path,
        )


def main():
    if len(sys.argv) != 2:
        print("Usage: python parse_json.py <path_to_json_file>")
        sys.exit(1)

    json_file_path = sys.argv[1]

    try:
        with open(json_file_path) as file:
            data = json.load(file)

        split_files(data, json_file_path)

    except FileNotFoundError:
        print(f"Error: File '{json_file_path}' not found.")
        sys.exit(1)
    except json.JSONDecodeError:
        print(f"Error: File '{json_file_path}' is not a valid JSON file.")
        sys.exit(1)


if __name__ == "__main__":
    main()
