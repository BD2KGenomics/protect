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
from __future__ import print_function

from collections import Counter

from protect.common import export_results, get_files_from_filestore, untargz
from protect.haplotyping.phlat import parse_phlat_file

import os
import pandas as pd
import textwrap

TCGAToGTEx = {'STAD': 'Stomach', 'KIRP': 'Kidney', 'THCA': 'Thyroid', 'PAAD': 'Pancreas',
              'KICH': 'Kidney', 'ESCA': 'Esophagus', 'CESC': 'Cervix Uteri',
              'PCPG': 'Adrenal Gland', 'BRCA': 'Breast', 'KIRC': 'Kidney', 'LUAD': 'Lung',
              'BLCA': 'Bladder', 'GBM': 'Brain', 'SKCM': 'Skin', 'LIHC': 'Liver',
              'COAD': 'Colon', 'LUSC': 'Lung', 'UCEC': 'Uterus', 'PRAD': 'Prostate'}

