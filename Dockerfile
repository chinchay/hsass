FROM julia:latest

RUN echo 'import Pkg; Pkg.add("PyCall"); Pkg.add("Conda")' >> xinstall.jl
RUN echo 'using Conda; Conda.add("pymatgen"; channel="conda-forge")' >> xinstall.jl
RUN julia xinstall.jl

# # this is for pyenv
# RUN apt-get update \
#     && apt-get upgrade -y \
#     && apt-get -yq --no-install-recommends install \
#     make build-essential libssl-dev zlib1g-dev \
#     libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev \
#     libncursesw5-dev xz-utils tk-dev libffi-dev liblzma-dev python-openssl \
#     git nano vim\
#     && apt-get autoremove \
#     && apt-get clean

# # https://github.com/jprjr/docker-pyenv/blob/master/Dockerfile
# RUN useradd -m python_user
# WORKDIR /home/python_user
# USER python_user
# RUN git clone git://github.com/yyuu/pyenv.git .pyenv
# ENV HOME  /home/python_user
# ENV PYENV_ROOT $HOME/.pyenv
# ENV PATH $PYENV_ROOT/shims:$PYENV_ROOT/bin:$PATH
# RUN pyenv install -v 3.8.0

# # # for global python environment:
# RUN pyenv global 3.8.0

# # # install python packages requirements:
# RUN pip install --upgrade pip
# RUN pip install pymatgen ase

# # RUN git clone https://username:password@gitlab.bcamath.org/cleon/hsass/tree/DGHS-algo3
# # RUN git clone https://gitlab.bcamath.org/cleon/hsass.git
# # RUN git checkout DGHS-algo3

COPY hsass-DGHS-algo3.tar .
RUN tar -xf hsass-DGHS-algo3.tar
RUN cd hsass-DGHS-algo3
WORKDIR hsass-DGHS-algo3




RUN echo 'import Pkg; Pkg.add("PyPlot"); Pkg.add("StatsBase")' >> xinstall.jl
RUN julia xinstall.jl

# RUN echo 'using Conda; Conda.add("jupyter")' >> xinstall.jl
# RUN echo 'import Pkg; Pkg.add("IJulia")' >> xinstall.jl

RUN echo 'import Pkg; Pkg.add("IJulia")' >> xinstall.jl
RUN julia xinstall.jl


# RUN pip3 install jupyter
# RUN echo 'import Pkg; Pkg.add("IJulia")' >> xinstall.jl
# RUN julia xinstall.jl


# # # Add Tini. Tini operates as a process subreaper for jupyter. This prevents kernel crashes.
# # ENV TINI_VERSION v0.6.0
# # ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
# # # RUN chmod +x /usr/bin/tini
# # CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]


# # RUN pip3 install ipython
# # RUN pip3 install qtconsole

# # CMD ["ipython", "notebook", "--profile", "julia", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]
# # CMD ["ipython", "notebook", "--profile", "julia", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]

# # CMD ["jupyter", "notebook", "--profile", "julia", "--port=8888", "--ip=0.0.0.0"]

# # CMD ["jupyter", "qtconsole", "--profile", "julia", "--port=8888", "--ip=0.0.0.0"]

# COPY testing.ipynb .
# CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]

# RUN echo "using IJulia; notebook()" >> oo.jl

# CMD ["julia", oo.jl, "--port=8888", "--ip=0.0.0.0", "--allow-root"]


ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    # && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
# RUN conda --version
RUN conda install -y jupyter
#COPY tutorial.ipynb .
CMD ["jupyter", "notebook", "--port=8888", "--ip=0.0.0.0", "--allow-root"]


# docker build -t dghsalgo3 .
# docker container run -p 8888:8888 dghsalgo3