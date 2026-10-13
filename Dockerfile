# MitelmanProgression (Shubhodeep Chatterjee): reproducible environment for the full pipeline.
#   docker build -t mitelman-progression .
#   docker run --rm -v "$PWD:/work" mitelman-progression bash run_all.sh
FROM python:3.13-slim
RUN apt-get update && apt-get install -y --no-install-recommends curl make \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /work
COPY requirements-lock.txt /tmp/requirements-lock.txt
RUN pip install --no-cache-dir -r /tmp/requirements-lock.txt
ENV OMP_NUM_THREADS=1
CMD ["bash", "run_all.sh"]
