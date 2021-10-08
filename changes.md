# Manual Changes to ProTECT 
:star: indicates changes to the algorithm that *could* potentially change results, though best efforts were made for 1:1 conversion  

:black_square_button: indicates changes i hope to reverse and are only 'temp fixes'

- Originally ran 2to3 in commit a5d062fab68f8bbbebc2bbe9f4192b47b451146e
	- removed explicit versioning in the [Makefile](https://github.com/BD2KGenomics/protect/commit/a5d062fab68f8bbbebc2bbe9f4192b47b451146e#diff-76ed074a9305c04054cdebb9e9aad2d818052b07091de1f20cad0bbac34ffb52) since py3 version still in dev
- :black_square_button: [removed version checks in setup.py](https://github.com/BD2KGenomics/protect/commit/f04f22fb9f50270e5c0307d4a64aca0f3f7022d3) and obsolete setuptools 
	- [along with setup version](https://github.com/BD2KGenomics/protect/commit/f70d3196198a2530406906b8af5a55b848aa0b14)
- [changed default references](https://github.com/BD2KGenomics/protect/commit/c2fe3a8b8223682e6d63cccb4fccf0787227c525) from s3://cgl-pipeline-inputs to s3://protect-data 
	- this s3 bucket is pay to access, however currently s3am is untested and only automatically converted 
- :star::star: [common.py chromosome sorting](https://github.com/BD2KGenomics/protect/commit/b5ca956f3dfe05bf6714be8135cc90fe48140d98)
-  [docker image decodes to utf-8](https://github.com/BD2KGenomics/protect/commit/1d2bdb941548bdf4703113140d1f0758791bf88a)
-  [IOBase check rather than file](https://github.com/BD2KGenomics/protect/commit/351e855184ae218242a42e1aaa5781d22aba0511) 
-  [some binary vs string adaptations](https://github.com/BD2KGenomics/protect/commit/5a4c50d1d2b8c71f3bc2f512f3679e80368044be#diff-e46b0e6e9cc33d9130334ab6994c9684b0972aaca58c889b6c1f4819751f1c79) 
-  :star: [change obselete ix to loc](https://github.com/BD2KGenomics/protect/commit/66bc12db0b815ab2099ee0174a06b923240322a4#diff-3347ae223ced4e929cf7f273bf839bdeb219d82681f8e66a951d85cbeb079685)
-  :black_square_button: quick fix: running into problem with '80.0' default for cores. can't figure out where the default is coming from so manual changes: 
	-  [phlat cores to 10](https://github.com/BD2KGenomics/protect/commit/66bc12db0b815ab2099ee0174a06b923240322a4#diff-7e85a3e4e9c911fded129ff48b2dd983d800c5190412f641eee85ff23ed9295c)
	-  [rsem cores to 20](https://github.com/BD2KGenomics/protect/commit/66bc12db0b815ab2099ee0174a06b923240322a4#diff-1615337ffdbffe39413f26e4ccbb5309ed10b61987559df23ce6fc57cb5dd86a) 
	-  [star fusion cores to 20](https://github.com/BD2KGenomics/protect/commit/66bc12db0b815ab2099ee0174a06b923240322a4#diff-60e2cfd2feabfe71442d69d1d0d44ff293f8fe1e12aa74c3fe52101d5b32e60e)
-  [string.maketrans is obselete, str.maketrans is better](https://github.com/BD2KGenomics/protect/commit/66bc12db0b815ab2099ee0174a06b923240322a4#diff-60e2cfd2feabfe71442d69d1d0d44ff293f8fe1e12aa74c3fe52101d5b32e60eR300) 
-  changed gunzip file write to use a library (faster?) 
