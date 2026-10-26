# Image for the docker and apptainer profiles: the tools of environment.yml,
# installed from the linux-64 lock (containers/update_lock.sh) so a rebuild
# gets the same packages.
FROM mambaorg/micromamba:2.3.0
USER root
RUN apt-get update && apt-get install -y --no-install-recommends procps && rm -rf /var/lib/apt/lists/*
COPY containers/linux-64.lock /tmp/linux-64.lock
RUN micromamba install -y -n base -f /tmp/linux-64.lock && micromamba clean -a -y
ENV PATH=/opt/conda/bin:$PATH
USER $MAMBA_USER
