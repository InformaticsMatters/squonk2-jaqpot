# FROM python:3.11.5-slim-bullseye
FROM python:3.11.6-slim-bullseye

ENV PYTHONUNBUFFERED=1
ENV HOME=/code

COPY requirements.txt /tmp/
RUN pip install --upgrade pip && \
    pip install --no-cache-dir --requirement /tmp/requirements.txt


WORKDIR ${HOME}
COPY src/ ./
