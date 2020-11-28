# Copyright 2016 UCSC Computational Genomics Lab
# Original contributor: Arjun Arkal Rao
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

define help

Supported targets: develop, sdist, clean, ci_test, test, pypi.

Please note that all build targets require a virtualenv to be active.

The 'develop' target creates an editable install of ProTECT and its runtime requirements in the
current virtualenv. The install is called 'editable' because changes to the source code
immediately affect the virtualenv.

The 'sdist' target creates a source distribution of ProTECT suitable for hot-deployment (not
implemented yet).

The 'clean' target undoes the effect of 'develop', and 'sdist'.

The 'ci_test' target runs the CI test for ProTECT. I.e. it runs the pipeline end-to-end using a test
dataset.

The 'test' target runs ProTECT's unit tests. These tests mandatorily require Toil to be installed in
the same virtualenv as ProTECT. Set the 'tests' variable to run a particular test, e.g.

	make test tests=src/protect/test/test_file_downloads.py::TestFileDownloads::test_file_downloads_from_s3

The 'pypi' target publishes the current commit of ProTECT to PyPI after enforcing that the working
copy and the index are clean, and tagging it as an unstable .dev build.

endef
export help
help:
	@echo "$$help"


python=python
pip=pip
tests=src/protect/test/unit
extras=
green=\033[0;32m
normal=\033[0m
red=\033[0;31m

# WIP 
special_install: check_venv
	git clone https://github.com/Dranion/bd2k-extras.git
	make -C bd2k-extras/bd2k-python-lib develop
	make -C bd2k-extras/s3am develop

prepare: check_venv
	@$(pip) install toil pytest  

develop: check_venv
	$(pip) install -e .$(extras)
clean_develop: check_venv
	- $(pip) uninstall -y protect
	- rm -rf src/*.egg-info


sdist: check_venv
	$(python) setup.py sdist
clean_sdist:
	- rm -rf dist

check_build_reqs:
	@$(python) -c 'import pytest; import toil' \
		|| ( echo "$(red)Build requirements (pytest or Toil) is missing. Run 'make prepare' to install them.$(normal)" ; false )
	@s3am --version \
		|| ( echo "$(red)Build requirement (s3am) is missing. Please install before running ProTECT.$(normal)" ; false )

check_toil_in_venv:
	@$(python) -c 'import toil; import os; assert toil.__file__.startswith(os.getcwd())' \
		|| ( echo "$(red)Build requirement (Toil) is not installed in the same venv as ProTECT. Install Toil in the venv before continuing.$(normal)" ; false )

test: check_venv check_toil_in_venv check_build_reqs
	$(python) -m pytest -vv -pyargs $(tests) --junitxml=test-report.xml

ci_test: check_venv check_build_reqs
	$(python) -m pytest -vv -pyargs src/protect/test/ci/test_protect.py -m 'all_fastq' --junitxml=test-report.xml

ci_mix_bam_fastq_test: check_venv check_build_reqs
	$(python) -m pytest -vv -pyargs src/protect/test/ci/test_protect.py -m 'mix_bam_fastq' --junitxml=test-report.xml

ci_vcf_fastq_test: check_venv check_build_reqs
	$(python) -m pytest -vv -pyargs src/protect/test/ci/test_protect.py -m 'vcf_fastq' --junitxml=test-report.xml

ci_test_all: ci_test ci_mix_bam_fastq_test ci_vcf_fastq_test

pypi: check_venv check_clean_working_copy check_running_on_jenkins
	set -x \
	&& tag_build=`$(python) -c 'pass;\
		from version import version as v;\
		from pkg_resources import parse_version as pv;\
		import os;\
		print "--tag-build=.dev" + os.getenv("BUILD_NUMBER") if pv(v).is_prerelease else ""'` \
	&& $(python) setup.py egg_info $$tag_build sdist bdist_egg upload
clean_pypi:
	- rm -rf build/


clean: clean_develop clean_sdist clean_pypi


check_venv:
	@$(python) -c 'import sys; sys.exit( int( not hasattr(sys, "base_prefix") ) )' \
		|| ( echo "$(red)A virtualenv must be active.$(normal)" ; false )


check_clean_working_copy:
	@echo "$(green)Checking if your working copy is clean ...$(normal)"
	@git diff --exit-code > /dev/null \
		|| ( echo "$(red)Your working copy looks dirty.$(normal)" ; false )
	@git diff --cached --exit-code > /dev/null \
		|| ( echo "$(red)Your index looks dirty.$(normal)" ; false )
	@test -z "$$(git ls-files --other --exclude-standard --directory)" \
		|| ( echo "$(red)You have are untracked files:$(normal)" \
			; git ls-files --other --exclude-standard --directory \
			; false )


check_running_on_jenkins:
	@echo "$(green)Checking if running on Jenkins ...$(normal)"
	@test -n "$$BUILD_NUMBER" \
		|| ( echo "$(red)This target should only be invoked on Jenkins.$(normal)" ; false )

docker:
	cd docker && make

clean_docker:
	cd docker && make clean

push_docker:
	cd docker && make push

.PHONY: help \
		develop clean_develop \
		sdist clean_sdist \
		test \
		pypi clean_pypi \
		clean \
		check_venv \
		check_clean_working_copy \
		check_running_on_jenkins \
		docker \
		clean_docker \
		push_docker
