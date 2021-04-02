# Copied from https://gitlab.com/ska-telescope/templates/cpp-template
# Generates the required metrics to generate the badges in a next step.
import json
import os
import subprocess
import sys
import time
import xml.dom.minidom


if sys.version_info[0] < 3:
    b2s = lambda s: s
else:
    b2s = lambda s: s.decode('utf-8')

def _has_subelement(element, subelement_name):
    return len(element.getElementsByTagName(subelement_name)) > 0

def _has_error(test_case):
    return _has_subelement(test_case, 'error')

def _has_failure(test_case):
    return _has_subelement(test_case, 'failure')

def get_git_sha():
    if 'CI_COMMIT_SHA' in os.environ:
        return os.environ['CI_COMMIT_SHA']
    return b2s(subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip())

def get_coverage_metrics(cov_xml_file):
    cov_dom = xml.dom.minidom.parse(cov_xml_file)
    coverage = cov_dom.getElementsByTagName('coverage')[0]
    coverage_line_rate = float(coverage.attributes['line-rate'].value)
    return {
        'percentage': coverage_line_rate * 100
    }

def get_tests_metrics(utests_xml_file):
    utests_xml = xml.dom.minidom.parse(utests_xml_file)
    test_cases = utests_xml.getElementsByTagName('testcase')
    errors = len(list(filter(_has_error, test_cases)))
    failures = len(list(filter(_has_failure, test_cases)))
    return {
        'errors': errors,
        'failed': failures,
        'total': len(test_cases)
    }

def get_lint_metrics(lint_xml_file):
    return {
        'errors': 0,
        'failures': 0,
        'tests': 0
    }

def get_build_status(ci_metrics):
    now = time.time()
    test_metrics = ci_metrics['tests']
    if test_metrics['errors'] > 0 or test_metrics['failed'] > 0:
        last_build_status = 'failed'
    else:
        last_build_status = 'passed'

    return {
        'last': {
            'status': last_build_status,
            'timestamp': now
        },
        'green': {
            'timestamp': now
        }
    }

def produce_ci_metrics(build_dir):
    cov_xml_file = os.path.join(build_dir, 'code-coverage.xml')
    utests_xml_file = os.path.join(build_dir, 'unit-tests.xml')
    lint_xml_file = os.path.join(build_dir, 'linting.xml')

    ci_metrics = {
        'commit_sha': get_git_sha(),
        'coverage': get_coverage_metrics(cov_xml_file),
        'tests': get_tests_metrics(utests_xml_file),
        'lint': get_lint_metrics(lint_xml_file)
    }
    ci_metrics['build-status'] = get_build_status(ci_metrics)
    print(json.dumps(ci_metrics, indent=2))

if __name__ == '__main__':
    produce_ci_metrics(sys.argv[1])