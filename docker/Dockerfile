FROM quay.io/ucsc_cgl/toil:3.5.1-eb5032ad339a5114b97bd64c6dcd041a078884b5

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

# Install ProTECT
RUN virtualenv /root/protectvenv --system-site-packages
RUN /root/protectvenv/bin/pip install protect==2.3.0

# Install bd2k-python-lib
RUN /root/protectvenv/bin/pip install bd2k-python-lib==1.14a1.dev29

# Copy relevant files to image folder
COPY wrapper.py /opt/pipeline/
COPY pipelineWrapper.py /opt/pipeline/

ENTRYPOINT ["python", "/opt/pipeline/wrapper.py"]
CMD ["--help"]
