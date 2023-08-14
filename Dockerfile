# start with linux/amd64 image that contains R, RStudio.
FROM rocker/rstudio:4.2.1

RUN apt-get update \ 
&& apt-get install -y libnetcdf-dev=1:4.7.3-1 \
&& apt-get install -y libxt-dev=1:1.1.5-1 \
&& apt-get install -y libcairo2-dev=1.16.0-4ubuntu1 \ 
&& apt-get install -y libxml2-dev=2.9.10+dfsg-5ubuntu0.20.04.6 \
&& apt-get install -y libharfbuzz-dev=2.6.4-1ubuntu4.2 \
&& apt-get install -y libfribidi-dev=1.0.8-2ubuntu0.1 \
&& apt-get install -y libfreetype6-dev=2.10.1-2ubuntu0.3  \
&& apt-get install -y libpng-dev=1.6.37-2  \
&& apt-get install -y libtiff5-dev=4.1.0+git191117-2ubuntu0.20.04.8 \ 
&& apt-get install -y libjpeg-dev=8c-2ubuntu8 \
&& apt-get install -y libglpk-dev=4.65-2 \
&& apt-get install -y --no-install-recommends build-essential libopenmpi-dev=4.0.3-0ubuntu1 mpi-default-dev=1.13

ENV RENV_VERSION=0.16.0
RUN Rscript -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" \
&& Rscript -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# brms backend (cmdstanr, cmdstan)
# specify versions (cmdstanr is commit of release 0.5.3)
ENV CMDSTAN_VERSION=2.32.2
ENV CMDSTANR_VERSION=22b391e68c9577bafcc0ae0721d8dc32a14e341b 
ENV CMDSTANR=stan-dev/cmdstanr@$CMDSTANR_VERSION
ENV CMDSTAN=/opt/cmdstan-$CMDSTAN_VERSION

RUN mkdir $CMDSTAN
WORKDIR /opt/

RUN Rscript -e "remotes::install_github('${CMDSTANR}', dependencies = TRUE)"
RUN Rscript -e "cmdstanr::install_cmdstan(version = '${CMDSTAN_VERSION}', dir = '${CMDSTAN}')"

# get renv in docker image
WORKDIR /home/rstudio
COPY renv.lock renv.lock

RUN mkdir -p renv
COPY renv/activate.R renv/activate.R
RUN Rscript -e "renv::restore()"

WORKDIR /home/rstudio
CMD ["/init"]