FROM bioconductor/bioconductor_docker:RELEASE_3_22

LABEL maintainer="Vy Dang <dang.vy.t@gmail.com>" \
      description="Reproduction of Tsukui et al. 2020 Fig 1 and 5 (GSE132771)" \
      version="1.0.0"

RUN mkdir -p /root/.R && printf 'CFLAGS += -w\nCXXFLAGS += -w\nCXX11FLAGS += -w\nCXX14FLAGS += -w\nCXX17FLAGS += -w\n' > /root/.R/Makevars

# hdf5:read .h5 file, glpk:graph clustering, gsl:Azimuth, bz2 & lzma:Azimuth
RUN apt-get update && apt-get install -y --no-install-recommends \
        libhdf5-dev libglpk-dev libgsl-dev libbz2-dev liblzma-dev libmagick++-dev \
    && rm -rf /var/lib/apt/lists/*

# Freeze CRAN to one date (see note below on choosing it)
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2026-09-17"))' \
      >> /usr/local/lib/R/etc/Rprofile.site

RUN Rscript -e "install.packages('remotes')"

COPY 00_install_packages.R /tmp/00_install_packages.R
RUN Rscript /tmp/00_install_packages.R && rm /tmp/00_install_packages.R