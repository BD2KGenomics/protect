FROM quay.io/ucsc_cgl/toil:3.5.0a1.dev321-6b22036e1bb4227c6d15f2aeda126dfb5cfab716

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

# Install Toil
RUN pip install toil[aws]==3.3.1

# Install s3am
RUN pip install s3am

# Install ProTECT
RUN pip install protect==2.3.0

# s3am requires an older version of bd2k-python-lib than ProTECT requires.
# But s3am will work with the newer version of bd2k-python-lib. Also, if
# this fix isn't applied after protect is installed, protect will get mad
# about something with s3am's version or whatever.
# Related: BD2KGenomics/protect#132

# Install bd2k-python-lib
RUN pip install bd2k-python-lib==1.14a1.dev29

# Copy relevant files to image folder
COPY wrapper.py /opt/pipeline/
COPY pipelineWrapper.py /opt/pipeline/

ENTRYPOINT ["python", "/opt/pipeline/wrapper.py"]
CMD ["--help"]
