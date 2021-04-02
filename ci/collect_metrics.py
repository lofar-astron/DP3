# Retrieved from https://gitlab.com/ska-telescope/ci-metrics-utilities/

"""Script to collect and parse CI metrics."""
import sys
import os
import urllib.request
import json
import logging
import xml.etree.ElementTree as etree
from datetime import datetime


LOGGER_NAME = "ci-metrics"


def parse_junit_tree(element, test_summary=None):
    """Recursively parse and extract a summary of a JUnit xml element tree.

    Args:
        element (xml.etree.Element): XML Element
        test_summary (dict): Summary of the JUnit XML report

    Returns:
        dict: Summary of Junit results.

    Notes:
        For JUnit XML schema see:
            <https://github.com/windyroad/JUnit-Schema>
    """
    # Create default entries in the stats dictionary if required.
    if not test_summary:
        test_summary = dict(
            testsuites=dict(tests=0, failures=0),
            testsuite=dict(tests=0, failures=0, skipped=0, errors=0),
            testcase=dict(tests=0, failures=0, skipped=0, errors=0),
        )

    # Parse top level <testsuites> or <testsuite> tags.
    if element.tag in ["testsuites", "testsuite"]:

        for attr in test_summary[element.tag]:
            test_summary[element.tag][attr] += int(element.get(attr, 0))

        for child in element:
            parse_junit_tree(child, test_summary)

    # Parse <testcase> tag. This is a child of testsuite.
    if element.tag == "testcase":
        key = "testcase"

        # Incrememnt test case counter.
        test_summary[key]["tests"] += 1

        # Parse child <error>, <skipped>, <failure> tags.
        for child in element:
            parse_junit_tree(child, test_summary)

    # Parse <error>, <skipped>, and <failure> tags. Children of testcase.
    if element.tag in ["error", "skipped", "failure"]:
        key = element.tag
        if element.tag == "error":
            key = "errors"
        elif element.tag == "failure":
            key = "failures"
        test_summary["testcase"][key] += 1

    # Parse <system-out>, <system-err>, and <properties> tags.
    # Children of testsuite.
    if element.tag in ["system-out", "system-err", "properties"]:
        pass

    return test_summary


def count_junit_metrics(filename):
    """Collect metrics from a JUnit XML file.

    Used to parse unit tests and linting results.

    Args:
        filename (str): Filename path of JUnit file

    Returns:
        dict: Summary of Unit test or linting results

    """
    log = logging.getLogger(LOGGER_NAME)
    try:
        root_elem = etree.parse(filename).getroot()
        if root_elem.tag not in ['testsuites', 'testsuite']:
            raise ValueError('Invalid JUnit XML file.')
        stats = parse_junit_tree(root_elem)
        result = dict(errors=0, failures=0, tests=0, skipped=0)
        for key in result:
            if key in ["tests", "failures"]:
                result[key] = max(
                    stats["testsuites"][key],
                    stats["testsuite"][key],
                    stats["testcase"][key],
                )
            else:
                result[key] = max(stats["testsuite"][key],
                                  stats["testcase"][key])
        result["total"] = result["tests"]
        del result["tests"]
    except Exception as expt:
        log.exception("Exception caught parsing '%s', returning 0 since the CI does not allow any linting errors/warnings", filename)
        result = dict(errors=0, failures=0, total=0, skipped=0)

    return result


def parse_coverage():
    """Parse the coverage report to return the percentage coverage.

    Returns:
        int or str: coverage percentage or 'unknown'
    """
    cov_tree = None
    cov_percent = "unknown"
    try:
        cov_tree = etree.parse("build/reports/code-coverage.xml")
    except FileNotFoundError:
        print("WARNING code-coverage.xml file not found")

    if cov_tree:
        try:
            cov_root = cov_tree.getroot()
            cov_percent = 100 * float(cov_root.get("line-rate"))
        except AttributeError as err:
            print(
                "WARNING: Attribute not found. Make sure that the file "
                "code-coverage.xml has the correct 'line-rate' attribute: ",
                err,
            )
    return cov_percent


def main():
    """Build ci-metric.json file from JUnit reports."""
    # pylint: disable=broad-except, too-many-locals
    log = logging.getLogger(LOGGER_NAME)

    # Exit if 'ci-metrics.json' already exists
    ci_metrics_fp = "build/reports/ci-metrics.json"
    if os.path.isfile(ci_metrics_fp):
        log.info(
            "ci-metrics.json file already exists, using the data "
            "present in that file."
        )
        sys.exit(0)

    # Read CI environment variables
    env_commit_sha = os.environ["CI_COMMIT_SHA"]
    env_project_id = os.environ["CI_PROJECT_ID"]
    env_default_branch = os.environ["CI_DEFAULT_BRANCH"]

    # Read pipeline status
    project_url = "https://gitlab.com/api/v4/projects/" + env_project_id
    project_pipelines_url = project_url + "/pipelines"

    # Setup defaults in case of problems getting the data
    latest_build_timestamp = "unknown"
    try:
        # Load data about last builds (pipelines) using GitLab API
        pipeline_data = urllib.request.urlopen(project_pipelines_url).read()
        project_pipelines = json.loads(pipeline_data)
    except Exception as err:
        log.warning("failed accessing pipeline data: %s", err)

    try:
        for pipeline in project_pipelines:
            if pipeline["ref"] == env_default_branch:
                # latest_pipeline_id = str(pipeline["id"])
                latest_build_date = pipeline["created_at"]
                latest_build_timestamp = datetime.timestamp(
                    datetime.strptime(latest_build_date,
                                      "%Y-%m-%dT%H:%M:%S.%fZ")
                )
                break
    except Exception as err:
        log.warning("failed to parse pipeline data: %s", err)

    # Parse coverage report
    coverage = parse_coverage()

    # Parse JUnit unit XML tests report
    test_result = count_junit_metrics("build/reports/unit-tests.xml")

    # Parse JUnit XML linting report
    lint_result = count_junit_metrics("build/reports/linting.xml")

    # Create and write data object with all the collected info
    ci_metrics_data = {
        "commit-sha": env_commit_sha,
        "build-status": {"last": {"timestamp": latest_build_timestamp}},
        "coverage": {"percentage": coverage},
        "tests": test_result,
        "lint": lint_result,
    }

    with open(ci_metrics_fp, "w") as write_file:
        json.dump(ci_metrics_data, write_file)


def init_logging(name="", level=logging.DEBUG):
    """Initialise python logger."""
    root = logging.getLogger(name)
    root.setLevel(level)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    root.addHandler(handler)


if __name__ == "__main__":
    init_logging()
    main()