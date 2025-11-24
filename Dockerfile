# Use Python 3.10 as base image
FROM python:3.10-slim

# Set working directory
WORKDIR /app

# Install system dependencies (needed for compiling some Python packages)
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    libc6-dev \
    libffi-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements.txt and install Python dependencies
COPY requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt

# Copy entire project to container
COPY . /app/

# Set Python path to ensure lotus package can be imported
ENV PYTHONPATH=/app:$PYTHONPATH

# Set environment variables for production mode
ENV GUNICORN=1
ENV DOCKER=1

# Expose port 8080
EXPOSE 8080

# Keep working directory at /app (project root) for Gunicorn to import lotus.api.app
# Gunicorn needs to run from project root to properly import the lotus package

# Run the application with Gunicorn production server
CMD ["gunicorn", "lotus.api.app:app", "--bind", "0.0.0.0:8080", "--workers", "4", "--timeout", "120", "--access-logfile", "-", "--error-logfile", "-"]

