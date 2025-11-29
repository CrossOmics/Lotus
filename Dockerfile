# Use Python 3.10 as base image
FROM python:3.10

# Set working directory
WORKDIR /app

# Create app user
RUN useradd -m -u 1000 appuser

# Copy requirements and install dependencies
COPY requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt

# Copy project files
COPY . /app/
RUN chown -R appuser:appuser /app

# Ensure lotus package is importable
ENV PYTHONPATH=/app:$PYTHONPATH

# Let Hugging Face provide the PORT
ENV PORT=7860

# Expose the port
EXPOSE ${PORT}

# Run the application with Gunicorn
CMD ["sh", "-c", "gunicorn lotus.api.app:app --bind 0.0.0.0:${PORT} --workers 1 --worker-class sync --timeout 600 --graceful-timeout 30 --keep-alive 5 --max-requests 50 --max-requests-jitter 5 --access-logfile - --error-logfile - --log-level info"]
