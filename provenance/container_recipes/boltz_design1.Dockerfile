# Use Miniconda as base image
FROM continuumio/miniconda3:latest

# set environment variables
env http_proxy="http://proxy.ethz.ch:3128"
env https_proxy="http://proxy.ethz.ch:3128"

# configure apt to use proxy
run echo 'acquire::http::proxy "http://proxy.ethz.ch:3128";' > /etc/apt/apt.conf.d/00proxy && \
    echo 'acquire::https::proxy "http://proxy.ethz.ch:3128";' >> /etc/apt/apt.conf.d/00proxy
# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# Set working directory
WORKDIR /app

# Set environment variables
ENV CONDA_DEFAULT_ENV=boltz_design
ENV PATH="/opt/conda/envs/boltz_design/bin:$PATH"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    wget \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Copy project files (assumes boltz directory and other files are in build context)
COPY . /app/

# Create conda environment and install dependencies
RUN echo "🚀 Setting up Boltz Design Environment..." && \
    conda create -n boltz_design python=3.10 -y && \
    echo "conda activate boltz_design" >> ~/.bashrc

# Activate environment and install packages
SHELL ["conda", "run", "-n", "boltz_design", "/bin/bash", "-c"]

# Install boltz
RUN if [ -d "boltz" ]; then \
        echo "📂 Installing Boltz..." && \
        cd boltz && \
        pip install -e . && \
        cd ..; \
    else \
        echo "❌ boltz directory not found. Please ensure boltz directory is in build context." && \
        exit 1; \
    fi

# Install conda dependencies
RUN echo "🔧 Installing conda dependencies..." && \
    conda install -c anaconda ipykernel -y

# Install Python dependencies
RUN echo "🔧 Installing Python dependencies..." && \
    pip install matplotlib seaborn prody tqdm PyYAML requests pypdb py3Dmol logmd==0.1.45

# Install PyRosetta
RUN echo "⏳ Installing PyRosetta (this may take a while)..." && \
    pip install pyrosettacolabsetup pyrosetta-installer && \
    python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'

# Download Boltz weights and dependencies
RUN echo "⬇️  Downloading Boltz weights and dependencies..." && \
    python -c "from boltz.main import download; from pathlib import Path; cache = Path('~/.boltz').expanduser(); cache.mkdir(parents=True, exist_ok=True); download(cache); print('✅ Boltz weights downloaded successfully!')"

# Setup LigandMPNN if directory exists
RUN if [ -d "LigandMPNN" ]; then \
        echo "🧬 Setting up LigandMPNN..." && \
        cd LigandMPNN && \
        bash get_model_params.sh "./model_params" && \
        cd ..; \
    fi

# Make DAlphaBall.gcc executable
RUN chmod +x "boltzdesign/DAlphaBall.gcc" || { echo "Error: Failed to chmod DAlphaBall.gcc"; exit 1; }

# Setup Jupyter kernel for the environment
RUN echo "📓 Setting up Jupyter kernel..." && \
    python -m ipykernel install --user --name=boltz_design --display-name="Boltz Design"

# Expose Jupyter port
EXPOSE 8888

# Set default command to activate environment and start bash
CMD ["conda", "run", "-n", "boltz_design", "/bin/bash"]
