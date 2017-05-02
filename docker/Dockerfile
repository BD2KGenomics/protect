FROM quay.io/ucsc_cgl/toil:3.5.2-378fffa320ded1ed1ebade5ec7d01138699db3f6

RUN apt-get update && apt-get install -y \
    git \
    python-dev \
    python-pip \
    wget \
    curl \
    apt-transport-https \
    ca-certificates \
    libcurl4-openssl-dev \
    libyaml-dev

# Get the Docker binary
RUN curl https://get.docker.com/builds/Linux/x86_64/docker-1.12.3.tgz \
         | tar -xvzf - --transform='s,[^/]*/,,g' -C /usr/local/bin/ \
         && chmod u+x /usr/local/bin/docker

#Upgrade pip
RUN pip install --upgrade pip setuptools

# Install Toil
RUN pip install toil[aws]==3.5.2

# Install s3am
RUN pip install s3am

# Install ProTECT
RUN pip install protect==2.4.1

# Install bd2k-python-lib
RUN pip install bd2k-python-lib==1.14a1.dev35

# Copy relevant files to image folder
COPY wrapper.py /opt/pipeline/
COPY pipelineWrapper.py /opt/pipeline/

ENTRYPOINT ["python", "/opt/pipeline/wrapper.py"]
CMD ["--help"]
