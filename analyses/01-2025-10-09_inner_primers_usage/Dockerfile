FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy the analysis script
COPY check_inner_primer_usage.py .

# Create a non-root user
RUN useradd -m -u 1000 analyzer
USER analyzer

# Set entrypoint
ENTRYPOINT ["python", "check_inner_primer_usage.py"]

# Add labels
LABEL maintainer="Eugenio Mattei <emattei@broadinstitute.org>"
LABEL description="Inner primer usage analysis tool for DC-TAP readout"
LABEL version="1.0"
