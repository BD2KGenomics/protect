#!/usr/bin/env python2.7
# Copyright 2016 Arjun Arkal Rao
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

# A lot of this code was taken from toil/test/src/__init__.py

from __future__ import absolute_import
import logging
import os
import tempfile
import unittest
import shutil
import re

from bd2k.util.files import mkdir_p

log = logging.getLogger(__name__)


class ProtectTest(unittest.TestCase):
    """
    A common base class for Protect tests. Please have every test case directly or indirectly
    inherit this one.

    When running tests you may optionally set the PROTECT_TEST_TEMP environment variable to the path
    of a directory where you want temporary test files be placed. The directory will be created
    if it doesn't exist. The path may be relative in which case it will be assumed to be relative
    to the project root. If PROTECT_TEST_TEMP is not defined, temporary files and directories will
    be created in the system's default location for such files and any temporary files or
    directories left over from tests will be removed automatically removed during tear down.
    Otherwise, left-over files will not be removed.
    """

    _temp_base_dir = None

    @classmethod
    def setUpClass(cls):
        super(ProtectTest, cls).setUpClass()
        cls._tempDirs = []
        temp_base_dir = os.environ.get('PROTECT_TEST_TEMP', None)
        if temp_base_dir is not None and not os.path.isabs(temp_base_dir):
            temp_base_dir = os.path.abspath(os.path.join(cls._projectRootPath(), temp_base_dir))
            mkdir_p(temp_base_dir)
        cls._temp_base_dir = temp_base_dir

    @classmethod
    def tearDownClass(cls):
        if cls._temp_base_dir is None:
            while cls._tempDirs:
                temp_dir = cls._tempDirs.pop()
                if os.path.exists(temp_dir):
                    shutil.rmtree(temp_dir)
        else:
            cls._tempDirs = []
        super(ProtectTest, cls).tearDownClass()

    def setUp(self):
        log.info("Setting up %s ...", self.id())
        super(ProtectTest, self).setUp()

    def _createTempDir(self, purpose=None):
        prefix = ['toil', 'test', self.id()]
        if purpose:
            prefix.append(purpose)
        prefix.append('')
        temp_dir_path = tempfile.mkdtemp(dir=self._temp_base_dir, prefix='-'.join(prefix))
        self._tempDirs.append(temp_dir_path)
        return temp_dir_path

    def tearDown(self):
        super(ProtectTest, self).tearDown()
        log.info("Tore down %s", self.id())

    def _getTestJobStorePath(self):
        path = self._createTempDir(purpose='jobstore')
        # We only need a unique path, directory shouldn't actually exist. This of course is racy
        # and insecure because another thread could now allocate the same path as a temporary
        # directory. However, the built-in tempfile module randomizes the name temp dir suffixes
        # reasonably well (1 in 63 ^ 6 chance of collision), making this an unlikely scenario.
        os.rmdir(path)
        return path

    def _getTestUnivOptions(self):
        return {'patient': 'test',
                'output_folder': self._createTempDir(purpose='pipeline_outputs'),
                'storage_location': 'local',
                'dockerhub': 'aarjunrao',
                'java_Xmx': '20G',
                'max_cores': 2,
                'ref': 'hg19'}

    @classmethod
    def _projectRootPath(cls):
        """
        Returns the path to the project root, i.e. the directory that typically contains the .git
        and src subdirectories. This method has limited utility. It only works if in "develop"
        mode, since it assumes the existence of a src subdirectory which, in a regular install
        wouldn't exist. Then again, in that mode project root has no meaning anyways.
        """
        assert re.search(r'__init__\.pyc?$', __file__)
        project_root_path = os.path.dirname(os.path.abspath(__file__))
        package_components = __name__.split('.')
        expected_suffix = os.path.join('src', *package_components)
        assert project_root_path.endswith(expected_suffix)
        project_root_path = project_root_path[:-len(expected_suffix)]
        return project_root_path


try:
    # noinspection PyUnresolvedReferences
    from _pytest.mark import MarkDecorator
except ImportError:
    # noinspection PyUnusedLocal
    def _mark_test(name, test_item):
        return test_item
else:
    def _mark_test(name, test_item):
        return MarkDecorator(name)(test_item)
