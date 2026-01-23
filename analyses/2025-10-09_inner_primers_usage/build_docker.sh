#!/bin/bash

# Build Docker image for primer usage analysis
docker build -t primer_usage_analyzer:latest .

echo "Docker image built successfully!"
echo "To run the container:"
echo "docker run --rm -v \$(pwd):/data primer_usage_analyzer:latest /data/input.bam /data/metadata.csv -o /data/output"
