FROM condaforge/mambaforge
COPY pxnotebook_env.yml /
#RUN conda install -c conda-forge mamba
RUN mamba env create --file /pxnotebook_env.yml -p /opt/conda/envs/pxnotebook && \
    mamba clean --all --yes
RUN apt-get update -qq && \
    apt-get install -y zip procps ghostscript
# Add conda installation dir to PATH
ENV PATH /opt/conda/envs/pxnotebook/bin:$PATH
# Dump the details of the installed packates to a file for posterity
RUN mamba env export --name pxnotebook > pxnotebook.yml
# Instruct R processes to use these empty files instead of clashing with a local config
RUN touch .Rprofile
RUN touch .Renviron
