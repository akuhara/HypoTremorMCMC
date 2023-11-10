FROM ubuntu:22.04

ARG USER=seismologist

RUN apt-get update && apt-get install -y \
  gfortran \
  liblapack-dev \
  libopenmpi-dev \
  libfftw3-dev \
  openmpi-bin \
  python3-pip \
  ssh

COPY . /usr/local/HypoTremorMCMC/

WORKDIR /usr/local/HypoTremorMCMC/src

RUN make FFTW="-I/usr/include -lfftw3"

RUN useradd -m ${USER}

USER ${USER}

WORKDIR /home/${USER}

ENV PATH=$PATH:/usr/local/HypoTremorMCMC/bin \
    OMPI_MCA_btl_vader_single_copy_mechanism=none

RUN cp -r /usr/local/HypoTremorMCMC/sample . \
  && mkdir /home/${USER}/wrk

CMD ["/bin/bash"] 

