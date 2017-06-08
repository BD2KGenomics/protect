# Definitions
docker_versions = 1.13.0 1.12.3 1.11.2 1.10.3 1.9.1 1.8.3 1.7.1 1.6.2
gen_directories = $(docker_versions:=.docker)
build_tool = runtime-container.DONE

#Get versions of toil's docker container, protect, toil, and s3am
protectversion := $(shell grep ^version ../src/protect/version.py | cut -d"'" -f2)
toilversion := $(shell grep ^toil_version ../setup.py | cut -d"'" -f2)
s3amversion :=$(shell grep ^s3am_version ../setup.py | cut -d"'" -f2)
toilcontainerversion := $(shell git clone --bare https://github.com/BD2KGenomics/toil && cd toil.git && git show-ref --dereference releases/$(toilversion) | tail -1 | cut -d" " -f1)

#Any pre-existing Toil repository in directory must be removed or else installation fails.
$(info The Toil version is $(toilversion)-$(toilcontainerversion))

ifeq (,$(toilcontainerversion))
    $(info Git tag for Toil not obtained correctly. This usually happens because Toil's git is already present in the directory.)
    $(info Try removing toil.git from the directory then trying again.)
    $(error Fatal error: Git tag for Toil not obtained correctly)
endif

#To allow building a particular dev version, such as dev404, for testing purposes,
#add a nested "ifneq ($(findstring dev404,$(protectversion)),dev404)"
ifeq ($(findstring a,$(protectversion)),a)
    $(info This version of ProTECT is $(protectversion))
    $(error Fatal error: Cannot build container for dev version)
endif

build: $(build_tool)

%.docker: dockerver=$(@:.docker=)
%.docker: Dockerfile.template Makefile.template wrapper.py pipelineWrapper.py
	mkdir -p $@
	sed -e "s:DOCKERVER:$(dockerver):g" Dockerfile.template > ${@}/Dockerfile
	sed -e "s:DOCKERVER:$(dockerver):g" Makefile.template > ${@}/Makefile
	sed -i "s/TOIL_DOCKER_CONTAINER/$(toilversion)-$(toilcontainerversion)/g" ${@}/Dockerfile
	sed -i "s/PROTECT_VERSION/$(protectversion)/g" ${@}/Dockerfile
	sed -i "s/TOIL_VERSION/$(toilversion)/g" ${@}/Dockerfile
	sed -i "s/S3AM_VERSION/$(s3amversion)/g" ${@}/Dockerfile
	cp wrapper.py ${@}/
	cp pipelineWrapper.py ${@}/
	rm -rf toil.git

generate: $(gen_directories)

$(build_tool): generate
	$(foreach subdir, $(gen_directories), cd $(subdir) && $(MAKE) build && cd ..;)
	touch $(build_tool)

push: build
	# Requires ~/.dockercfg
	$(foreach subdir, $(gen_directories), cd $(subdir) && $(MAKE) push && cd ..;)

test: build test.py
	$(foreach subdir, $(gen_directories), cd $(subdir) && $(MAKE) test && cd ..;)

clean:
	-rm -rf $(gen_directories) && rm $(build_tool)
