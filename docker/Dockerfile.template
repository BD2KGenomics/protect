FROM quay.io/ucsc_cgl/toil:TOIL_DOCKER_CONTAINER

RUN apt-get update && apt-get install -y \
    git \
    python-dev \
    python-pip \
    wget \
    curl \
    apt-transport-https \
    ca-certificates \
    libcurl4-openssl-dev \
    libyaml-dev \
    zip

# Get the Docker binary
RUN curl https://get.docker.com/builds/Linux/x86_64/docker-DOCKERVER.tgz \
         | tar -xvzf - --transform='s,[^/]*/,,g' -C /usr/local/bin/ \
         && chmod u+x /usr/local/bin/docker

#Upgrade pip
RUN pip install --upgrade pip setuptools

# Install Toil
RUN pip install toil[aws]==TOIL_VERSION

# Install s3am
RUN pip install s3am==S3AM_VERSION

# Install GDC Client
RUN cd /usr/local/bin \
    && wget -c https://gdc.cancer.gov/files/public/file/gdc-client_v1.1.0_Ubuntu14.04_x64.zip \
    && unzip gdc-client_v1.1.0_Ubuntu14.04_x64.zip \
    && chmod +x gdc-client \
    && rm gdc-client_v1.1.0_Ubuntu14.04_x64.zip

# Install ProTECT
RUN pip install protect==PROTECT_VERSION


# s3am requires an older version of bd2k-python-lib than Toil requires.
# But s3am will work with the newer version of bd2k-python-lib.
# Related: BD2KGenomics/protect#132

# Install bd2k-python-lib
RUN pip install bd2k-python-lib==1.14a1.dev35


# Copy relevant files to image folder
COPY wrapper.py /opt/pipeline/
COPY pipelineWrapper.py /opt/pipeline/

ENTRYPOINT ["python", "/opt/pipeline/wrapper.py"]
CMD ["--help"]
