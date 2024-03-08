FROM python:3.10

COPY ./scanpy/requirements.txt /requirements.txt
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    cmake gcc gfortran libopenblas-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip install --upgrade pip setuptools distlib --no-cache-dir && \
    pip install --no-cache-dir -r /requirements.txt

ENTRYPOINT []