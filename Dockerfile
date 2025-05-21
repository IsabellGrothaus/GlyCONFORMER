# syntax=docker/dockerfile:1

ARG PYTHON_VERSION=3.10
FROM python:${PYTHON_VERSION}-slim as base

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Create a non-privileged user that the app will run under.
#ARG UID=10001
#RUN adduser \
#    --disabled-password \
#    --gecos "" \
#    --home "/nonexistent" \
#    --shell "/sbin/nologin" \
#    --no-create-home \
#    --uid "${UID}" \
#    appuser

# Install system dependencies for Plumed and notebook rendering
RUN apt-get update && apt-get install -y \
    build-essential \
    libblas-dev \
    liblapack-dev \
    libgfortran5 \
    libxrender1 \
    libglib2.0-0 \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install pip, setuptools, and wheel
RUN pip install --upgrade pip setuptools wheel

# Switch to the non-privileged user to run the application.
#USER appuser

# Set working directory
WORKDIR /app

# Copy the source code into the container.
COPY . /app

# Install dependencies using pyproject.toml
RUN pip install .

# Install Jupyter Notebook
RUN pip install notebook

# Expose the port that the application listens on.
EXPOSE 8888

# Run the application.
CMD ["jupyter", "notebook", "--ip=0.0.0.0", "--port=8888", "--NotebookApp.token=''", "--NotebookApp.password=''", "--no-browser", "--allow-root", "--NotebookApp.notebook_dir=/app"]
