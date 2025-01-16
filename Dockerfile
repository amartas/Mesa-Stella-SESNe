FROM ubuntu:24.04

LABEL maintainer="aidmart@berkeley.edu"

# Change shell to bash such that source works
SHELL ["/bin/bash", "-l", "-c"]

# Remove the interactive requirement for building everything
RUN sed -i '/\[ -z "\$PS1" \] && return/d' $HOME/.bashrc

# Set up the instance
RUN apt-get update -qq && apt-get install git nano unzip curl wget -y && \
    cd ~/ && \
    mkdir Setup && \
    cd Setup && \
    curl -O https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh && \
    bash Anaconda3-2024.10-1-Linux-x86_64.sh -b -p $HOME/anaconda3 && \
    source ~/anaconda3/bin/activate && \
    echo "source ~/anaconda3/bin/activate" >> ~/.bashrc && \
    conda create -n mesaenv python=3.12 -y && \
    conda activate mesaenv && \
    pip install numpy pandas -q -q -q

# Install the MESA SDK and its dependencies
RUN cd ~/Setup && \
    apt install binutils make perl libx11-6 libx11-dev zlib1g zlib1g-dev tcsh -y && \
    wget http://user.astro.wisc.edu/~townsend/resource/download/mesasdk/mesasdk-x86_64-linux-24.7.1.tar.gz --user-agent="" && \
    tar xvfz mesasdk-x86_64-linux-24.7.1.tar.gz -C ~/ && \
    export MESASDK_ROOT=~/mesasdk && \
    echo "export MESASDK_ROOT=~/mesasdk" >> ~/.bashrc && \
    echo "source $MESASDK_ROOT/bin/mesasdk_init.sh" >> ~/.bashrc && \
    source ~/.bashrc

# Get MESA
RUN cd ~/ && \
    wget https://zenodo.org/records/13353788/files/mesa-24.08.1.zip?download=1 && \
    unzip -qq 'mesa-24.08.1.zip?download=1'

# Get Mesa-Stella-SESNe
RUN cd ~/ && \
    git clone https://github.com/amartas/Mesa-Stella-SESNe.git -q && \
    mv Mesa-Stella-SESNe/ MESA && \
    mv mesa-24.08.1 MESA/mesa-24.08.1

# Install MESA and Mesa-Stella-SESNe
RUN cd ~/MESA && \
    ls && \
    pwd && \
    bash allow_root && \
    echo "export MESA_DIR=/root/MESA/mesa-24.08.1" >> ~/.bashrc && \
    echo "export OMP_NUM_THREADS=2" >> ~/.bashrc && \
    echo "export PATH=$PATH:$MESA_DIR/scripts/shmesa" >> ~/.bashrc && \
    echo "echo 'bashrc sourced'" >> ~/.bashrc && \
    source ~/.bashrc && \
    cd $MESA_DIR && \
    pwd && \
    gfortran --version && \
    bash install && \
    cd .. && \
    cp -r ModelGrids $MESA_DIR/ModelGrids && \
    mkdir Logs && \
    mkdir DataExports && \
    mkdir ProgOptimize && \
    chmod -R 777 $MESA_DIR/ModelGrids/000*

# Add the interactive requirement back
RUN sed -i '1i [ -z "$PS1" ] && return' $HOME/.bashrc

CMD ["bash"]
