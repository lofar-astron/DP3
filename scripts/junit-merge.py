#!/usr/bin/env python3
import sys
from xml.etree import ElementTree


def run(output_file, input_files):
    suites = ElementTree.Element("testsuites")
    for filename in input_files:
        data = ElementTree.parse(filename).getroot()
        suites.extend(data.iter("testsuite"))
    ElementTree.ElementTree(suites).write(output_file, xml_declaration=True)


if __name__ == "__main__":
    run(sys.argv[1], sys.argv[2:])
