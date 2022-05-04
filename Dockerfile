# based on existing Docker image
FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive

# dependencies, some are probably unnecessary
RUN apt-get update && apt-get install -y wget && apt-get install -y --no-install-recommends build-essential r-base python3.9 python3-pip python3-setuptools python3-dev
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    r-base-dev \
    libgsl0-dev \ 
    libxml2-dev \
    libcairo2-dev \
    libsqlite-dev \
    libpq-dev \
    libicu-dev \
    libbz2-dev \
    liblzma-dev \
    libfontconfig1-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libnetcdf-dev \
    udunits-bin \
    libopenblas-dev \ 
    libudunits2-dev \
    curl
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    autoconf \
    automake \
    g++ \
    gcc \
    gfortran \
    make \
    && apt-get clean all \
    && rm -rf /var/lib/apt/lists/*

# Python packages
WORKDIR /app
COPY requirements.txt /app/requirements.txt
RUN pip3 install -r requirements.txt

# R packages
RUN LC_ALL=C.UTF-8 Rscript -e "install.packages('data.table')"
RUN LC_ALL=C.UTF-8 Rscript -e "install.packages('logger')"
RUN LC_ALL=C.UTF-8 Rscript -e "install.packages('Rfast')"


COPY . /app


